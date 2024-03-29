"""=================================================================================================
orthogroups:og_interpro.py

send selected orthogroups to interpro to find domains and GO terms

04 March 2024     gribskov
================================================================================================="""
import argparse
import time
import datetime
import sys
import pickle
from os.path import basename, exists
from os import mkdir
# import os

from sequence.fasta import Fasta
from interpro.interpro import Interpro


def process_command_line():
    """---------------------------------------------------------------------------------------------
    Read command line options
    :return:
    ---------------------------------------------------------------------------------------------"""
    cl = argparse.ArgumentParser(
        description='Run interproscan on selected orthogroups',
        formatter_class=lambda prog: argparse.HelpFormatter(prog, width=120, max_help_position=40)
        )
    cl.add_argument('-g', '--orthogroup',
                    help='Orthogroups file',
                    type=str,
                    default='<not specified>')

    cl.add_argument('-o', '--out',
                    help='output directory for results',
                    type=str,
                    default='./')

    cl.add_argument('-s', '--skip',
                    help='skip (default: %(default)s)',
                    action='store_true',
                    default=False)

    # cl.add_argument('-n', '--ntop',
    #                 help='number of top/bottom groups to report',
    #                 type=int,
    #                 default=20)

    return cl.parse_args()


def read_list(opt, ogs):
    """---------------------------------------------------------------------------------------------
    Read the list of selected orthogroups. Use paths relative to the working directory or absolute
    paths. Blank lines or lines beginning with # are skipped.
    Exit with status = 0 if opt.orthogroup file cannot be read

    :param opt: namespace       command line options from process_command_line()
    :param ogs: dict            keys are orthogroup filepaths from input list
    :return: int                number of orthogroups read
    ---------------------------------------------------------------------------------------------"""
    oglist = opensafe(opt.orthogroup, 'r')

    og_n = 0
    for line in oglist:
        if line.startswith('#') or not line.strip():
            # skip blank lines and comments
            continue

        ogs.append(line.rstrip())
        og_n += 1

    oglist.close()

    return og_n


def opensafe(filename, mode):
    """---------------------------------------------------------------------------------------------
    open a file with error check. Failure results in exit(1)

    :param filename: string     path to file
    :param mode: string         open mode, typically 'r' or 'w'
    :return: filehandle         open filehandle
    ---------------------------------------------------------------------------------------------"""
    try:
        fh = open(filename, mode)
    except OSError:
        sys.stderr.write(f'\nUnable to open orthogroup list ({opt.orthogroup}\n\n')
        exit(1)

    return fh


def interpro_setup(fasta, og, pickledir):
    """---------------------------------------------------------------------------------------------
    creates a new Interpro object and loads the parameters for one search. The orthogroup is
    prefixed to the title so all results from the same og begin the same. The prefixed ID is stored
    in the Interpro object in the title attribute.

    :param fasta: Fasta         query sequence object
    :param og: string           orthogroup for this sequence
    :param pickledir: string    directory path for pickle output
    :return: Interpro           Interpro search object
    ---------------------------------------------------------------------------------------------"""
    ips = Interpro()
    ips.email = 'mgribsko@purdue.edu'
    ips.application_select(['PfamA', 'SMART', 'PrositePatterns', 'CDD', 'NCBIfam', 'PIRSF', 'SuperFamily'])
    ips.output_select('json')
    ips.parameter_select({'goterms': True, 'pathways': False})
    ips.title = f'{og}_{fasta.id}'
    fasta.seq = fasta.seq.rstrip('*')  # interproscan doesn't like * at the end
    ips.sequence = fasta.format()
    # get the output file name at the beginning in order to test if it already exists
    ips.outputfile = pickledir + ips.title.replace('|', '_') + '.pkl'
    # ips.pickle = os.path.join(pickledir, ips.title.replace('|', '_') + '.pkl')

    return ips


def next_og_sequence(ogs, fasta):
    """---------------------------------------------------------------------------------------------
    Generator for the main loop, iterates over all OGs and all sequences in each OG.
    Yields a fasta object with the next query sequence, and its OG.
    If all orthogroups have been processed, returns None.

    :param ogs: list        strings with OG sequence file names
    :param fasta: Fasta     Fasta object with query sequences from current OG
    :return: Fasta, og      next query sequence, or None if all have been processed
                            og is a string with the name of the current Orthogroup
                            or None when finished
    ---------------------------------------------------------------------------------------------"""
    og = ''
    query = None
    fasta_iter = None

    while True:

        # if an iterator is set try to get a sequence
        if fasta_iter:
            # print(f'iter exists {fasta_iter}')
            try:
                # success
                query = next(fasta_iter)
            except StopIteration:
                # failed to get a sequence
                # keeps StopIteration from crashing program
                query = None
                fasta_iter = None

            if query:
                # print(f'query found ({query})')
                yield fasta, og
                continue

        # no query could be obtained from current iterator, the current file is finished, look for the
        # next one in ogs

        if ogs:
            # OG list has more files
            ogfilename = ogs.pop()
            # print(f'new og {ogfilename}')
            fasta = Fasta(filename=ogfilename, mode='r')
            og = basename(ogfilename)
            og = og[:og.rindex('.')]
            fasta_iter = iter(fasta)
            query = next(fasta_iter)
            yield fasta, og
            # print(f'iter exists {fasta_iter}')

        # no more sequences in fasta file, and no more OGs to process - we are done, break out of
        # forever loop
        if not query:
            break

    # end of forever loop

    return None


class InterproscanManager:
    """=============================================================================================
    handle the sumbmission, polling, retrieval, and processing of searches
    ============================================================================================="""
    n = 0

    def __init__(self, opt, batch_limit=30, poll_time=30, poll_max=50, pkl='pkl'):
        """-----------------------------------------------------------------------------------------
        Attributes
            submitted = []      list of submitted Interpro objects
            finished = []       list of completed Interpro objects
            failed = []         list of failed Interpro objects
            save = []           list of jobs ready to pickled and save to disk
            opt                 command line options (namespace from argparse)
            batch_limit = 3     maximun number of jobs to submit in a batch
            poll_time = 30      sleep time after polling
            poll_max = 50       maximum times to poll before giving up (job is failed)
            pkl                 directory for pickled output files
        -----------------------------------------------------------------------------------------"""
        self.submitted = []
        self.finished = []
        self.failed = []
        self.save = []
        self.opt = opt
        self.batch_limit = batch_limit
        self.poll_time = poll_time
        self.poll_max = poll_max
        self.pkl = f'{pkl}/'
        if not exists(pkl):
            # make sure output directory exists
            mkdir(pkl)
        self.log_fh = InterproscanManager.getlog('og_interpro.log')

    def submit(self, fasta, og):
        """-----------------------------------------------------------------------------------------
        Submit jobs to interproscan service, up to self.batch_limit can be queued. interpro_setup()
        creates a new Interpro object.

        :param fasta: Fasta             sequence for submission
        :param og: string               orthogroup for this sequence
        :return: Bool                   False if skipped, True if submitted
        -----------------------------------------------------------------------------------------"""
        ip_submitted = self.submitted
        ip_failed = self.failed
        if len(ip_submitted) < self.batch_limit:
            ips = interpro_setup(fasta, og, self.pkl)

            if self.opt.skip and exists(ips.outputfile):
                # this sequence exists in the output, skip
                print(f'\tskipping {og}:{fasta.id}')
                return False

            print(f'\tJob {self.n} - {og}:{fasta.id}')
            self.n += 1  # total number of jobs submitted (class variable)

            tries = 1
            success = ips.submit(show_query=False)
            while not success and tries < 3:
                # try three times to submit with 5 seconds between tries
                time.sleep(5)
                success = ips.submit(show_query=False)
                tries += 1

            if success:
                # success
                self.log('SUBMIT', ips.title)
                ip_submitted.append(ips)
            else:
                # failure
                ip_failed.append(ips)
                sys.stderr.write(f'{ips.title} failed\n')
                self.log('FAIL-SUB', ips.title)

        return True

    def poll(self):
        """-----------------------------------------------------------------------------------------
        EBI requests that no further jobs be submitted until all have been finished. Keep polling
        until all the jobs in ip_submitted have been completed and moved to ip_finished. Note that
        you only have to poll until you find the first job that is still running.
        -----------------------------------------------------------------------------------------"""
        ip_submitted = self.submitted
        ip_finished = self.finished
        ip_failed = self.failed

        if not ip_submitted:
            # no jobs in submitted queue
            sys.stderr.write('No submitted jobs\n')
            return

        tries = 0
        while ip_submitted:
            ips = ip_submitted.pop()
            self.log('POLL', ips.title)
            tries += 1

            if ips.status() == 'finished':
                # remove job from ip_submitted and add to ip_finished
                ip_finished.append(ips)
                continue

            # if we reach here at least one job is not finished, push ips back on the stack
            ip_submitted.append(ips)

            if tries > self.poll_max:
                # polled too many times
                ip_failed.append(ips)
                self.log('FAIL-MAX', ips.title)

            time.sleep(self.poll_time)

        return

    def getresult(self):
        """-----------------------------------------------------------------------------------------

        :return:
        -----------------------------------------------------------------------------------------"""
        finished = self.finished
        if not finished:
            # sys.stderr.write('no finished jobs to retrieve\n')
            return

        # parse finished jobs and extract desired information
        while finished:
            thisjob = finished.pop()
            self.log('RETRIEVE', thisjob.title)
            thisjob.result()
            self.save.append(thisjob)
            self.save_as_pickle(thisjob)
            # print(thisjob.content)

        # a cooldown period appears to be necessary
        # time.sleep(75)
        return len(self.save)

    def log(self, message, jobtitle):
        """-----------------------------------------------------------------------------------------
        Write a message to the log file (self.log_fh). Format is
            Time/date\tmessage\tjob for example
            2024-03-11 10:24:54    SUBMITTED   OG0000066:ji|Acastr1|120024|AST1_006320-RA

        :param message: string      message content
        :param jobtitle: string     Interpro.title - job title
        :return: None
        -----------------------------------------------------------------------------------------"""
        timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        self.log_fh.write(f'{timestamp}\t{message:>8s}\t{jobtitle}\n')
        self.log_fh.flush()

        return None

    @staticmethod
    def save_as_pickle(ips):
        """-----------------------------------------------------------------------------------------
        write a pickle file with the interproscan result. the filename is set up in interpro_setup()
        to check whether it already exists when skipping files

        :param ips: Interpro        object with search result
        :return: True
        -----------------------------------------------------------------------------------------"""
        picklename = ips.outputfile
        with open(picklename, 'wb') as picklefile:
            pickle.dump(ips, picklefile, pickle.HIGHEST_PROTOCOL)

        return True

    @staticmethod
    def getlog(logfname):
        """-----------------------------------------------------------------------------------------
        try to open a new logfile. if file already exists add a header, if does not exist open in
        append mode. Exit with status==2 if file can't be opened

        :param logfname: string      path to logfile
        :return: filehandle          file open for append
        -----------------------------------------------------------------------------------------"""
        try:
            logfile = open(logfname, 'a')
        except IOError:
            sys.stderr.write(f'InterproscanManager:getlog - Error opening log file ({logfname}\n\n')
            exit(2)

        return logfile


# ==================================================================================================
# Main
# ==================================================================================================
if __name__ == '__main__':
    daytime = datetime.datetime.now()
    runstart = daytime.strftime('%Y-%m-%d %H:%M:%S')
    opt = process_command_line()

    sys.stderr.write(f'\nog_interpro.py {runstart}\n')
    sys.stderr.write(f'\tOG list: {opt.orthogroup}\n')
    sys.stderr.write(f'\tOutput directory: {opt.out}\n')

    # Read list of selected orthogroups, initially each orthogroup points to an empty list,
    # when we read the sequences these will be filled with SequenceInfo objects
    ogs = []
    og_n = read_list(opt, ogs)
    s = 's'
    if og_n == 1:
        s = ''
    sys.stderr.write(f'\n{og_n} orthogroup{s} read from {opt.orthogroup}\n')

    # set up the interproscan searches, all the searches can be done through a single object

    ips_manager = InterproscanManager(opt, batch_limit=20, pkl=opt.out)
    ips_manager.log('\nBEGIN', f'{runstart} Query={opt.orthogroup}')

    fasta = Fasta()
    og_old = ''
    og_n = 0
    og_job_n = 0    # number of jobs in OG
    for sequence, og in next_og_sequence(ogs, fasta):
        # print(f'sub:{len(ips_manager.submitted):}')
        # print(f'finished:{len(ips_manager.finished):}')
        # print(f'save:{len(ips_manager.save):}')
        # print(f'fail:{len(ips_manager.failed):}')
        if og != og_old:
            og_n += 1
            if og_job_n != 0:
                print(f'{og_job_n} sequences submitted for orthogroup {og}')
            print(f'\nStarting orthogroup {og_n}: {og}')
            og_old = og
            og_job_n = 0

        # send sequences to interproscan in groups of batch_size, waiting for batches to complete
        # use one Interpro object for each sequence
        og_job_n += ips_manager.submit(sequence, og)

        # only leave job submission if the desired number of jobs have been submitted
        if len(ips_manager.submitted) >= ips_manager.batch_limit:
            ips_manager.poll()
            ips_manager.getresult()

    # get any remaining job results
    if len(ips_manager.submitted):
        ips_manager.poll()
        ips_manager.getresult()

    sys.stderr.write(f'{ips_manager.n} jobs submitted\n')
    if ips_manager.failed:
        sys.stderr.write(f'Failed jobs\n')
        for ips in ips_manager.failed:
            sys.stderr.write(f'{ips.sequence.id}\n')
    else:
        sys.stderr.write(f'No Failed jobs\n')

    exit(0)
