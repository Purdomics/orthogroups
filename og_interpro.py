"""=================================================================================================
orthogroups:og_interpro.py

send selected orthogroups to interpro to find domains and GO terms

04 March 2024     gribskov
================================================================================================="""
import argparse
import time
import datetime
import sys
from os.path import basename

from sequence.fasta import Fasta
from interpro.interpro import Interpro


def process_command_line():
    """---------------------------------------------------------------------------------------------
    Read command line options
    :return:
    ---------------------------------------------------------------------------------------------"""
    cl = argparse.ArgumentParser(
        description='Run interpro on selected orthogroups',
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


def interpro_setup(fasta):
    """---------------------------------------------------------------------------------------------
    creates a new Interpro object and loads the parameters for one search.

    :param fasta: Fasta         query sequence object
    :return: Interpro           Interpro search object
    ---------------------------------------------------------------------------------------------"""
    # TODO move to manager or drop?
    ips = Interpro()
    ips.email = 'mgribsko@purdue.edu'
    ips.application_select(['PfamA', 'SMART', 'PrositePatterns', 'CDD', 'NCBIfam', 'PIRSF', 'SuperFamily'])
    ips.output_select('json')
    ips.parameter_select({'goterms': True, 'pathways': False})
    ips.title = fasta.id
    fasta.seq = fasta.seq.rstrip('*')  # interproscan doesn't like * at the end
    ips.sequence = fasta.format()

    return ips


def next_og_sequence(ogs, fasta):
    """---------------------------------------------------------------------------------------------
    Generator for the main loop, iterates over all OGs and all sequences in each OG.
    Yields a fasta object with the next query sequence, and its OG.
    If all orthogroups have
    been processed, returns None.

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
                yield fasta, og

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

    def __init__(self, batch_limit=30, poll_time=30, poll_max=50):
        """-----------------------------------------------------------------------------------------
        Attributes
            submitted = []      list of submitted Interpro objects
            finished = []       list of completed Interpro objects
            failed = []         list of failed Interpro objects
            batch_limit = 3     maximun number of jobs to submit in a batch
            poll_time = 30      sleep time after polling
            poll_max = 50       maximum times to poll before giving up (job is failed)
        -----------------------------------------------------------------------------------------"""
        self.submitted = []
        self.finished = []
        self.failed = []
        self.batch_limit = batch_limit
        self.poll_time = poll_time
        self.poll_max = poll_max

    def submit(self, fasta):
        """---------------------------------------------------------------------------------------------
        submit jobs to interproscan service, up to submit_max jobs can be queued

        :param fasta: Fasta             sequence for submission
        :return:
        ---------------------------------------------------------------------------------------------"""
        ip_submitted = self.submitted
        ip_failed = self.failed
        if len(ip_submitted) < self.batch_limit:
            ips = interpro_setup(fasta)
            print(f'submitting {og}:{fasta.id}')

            if not ips.submit(show_query=False):
                ip_failed.append(ips)
            else:
                ip_submitted.append(ips)

        return

    def poll(self):
        """-----------------------------------------------------------------------------------------
        EBI requests that no further jobs be submitted until all have been finished. Keep polling
        until all the jobs in ip_submitted have been completed and moved to ip_finished. Note that
        you only have to poll until you find the first job that is still running.
        -----------------------------------------------------------------------------------------"""
        ip_submitted = self.submitted
        ip_finished = self.finished

        if not ip_submitted:
            # no jobs in submitted queue
            return

        tries = 0
        i = 0
        while ip_submitted:
            ips = ip_submitted[i]
            tries += 1

            if ips.status() == 'finished':
                # remove job from ip_submitted and add to ip_finished
                ip_finished.append(ips)
                ip_submitted.remove(ips)
                i += 1
                continue

            if tries > self.maxtries:
                ip_finished.append(ips)
                ip_submitted.remove(ips)

            i = 0
            time.sleep(self.poll_time)

        return


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

    # setup the interproscan searches, all the searches can be done through a single object

    ips_manager = InterproscanManager(batch_limit=3)

    done = False
    fasta = Fasta()
    for sequence, og in next_og_sequence(ogs, fasta):
        s = sequence.format()
        print(og, s)

        # send sequences to interproscan in groups of batch_size, waiting for batches to complete
        # use one Interpro object for each sequence

        if len(ips_manager.submitted) < ips_manager.batch_limit:
            ips_manager.submit(sequence)

        ips_manager.interpro_poll()

        # if ips_finished:
        #     # parse finished jobs and extract desired information
        #     print('collecting result')
        #     ips.result()
        #     print(ips.content)

    exit(0)
