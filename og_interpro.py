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
    fh = None
    try:
        fh = open(filename, mode)
    except OSError:
        sys.stderr.write(f'\nUnable to open orthogroup list ({opt.orthogroup}\n\n')
        exit(1)

    return fh


def interpro_setup(title, fasta):
    """---------------------------------------------------------------------------------------------
    creates a new Interpro object and loads the parameters for one search.

    :param title: string        job title for interproscan
    :param fasta: Fasta         query sequence object
    :return: Interpro           Interpro search object
    ---------------------------------------------------------------------------------------------"""
    ips = Interpro()
    ips.email = 'mgribsko@purdue.edu'
    ips.application_select(['PfamA', 'SMART', 'PrositePatterns', 'CDD', 'NCBIfam', 'PIRSF', 'SuperFamily'])
    ips.output_select('json')
    ips.parameter_select({'goterms': True, 'pathways': False})
    ips.title = fasta.id
    ips.sequence = fasta.format()

    return ips


def interpro_submit(ips, batch_max, submitted):
    """---------------------------------------------------------------------------------------------
    If there is
    :param ips:
    :param batch_max:
    :param submitted:
    :return:
    ---------------------------------------------------------------------------------------------"""
    return


def next_og_sequence(ogs, fasta):
    """---------------------------------------------------------------------------------------------
    return a fasta object with the next query sequen e. If all orthogroups have been processed,
    return None

    :param ogs: list        strings with OG sequence file names
    :param fasta: Fasta     Fasta object with query sequences from current OG
    :return: Fasta          next query sequence, or None if all have been processed
    ---------------------------------------------------------------------------------------------"""
    query = None
    if fasta.fh:
        # already have an open file, try to read a sequence
        # query = fasta.next()
        if query:
            return query

    # fall through to hear if filehandle isn't open, or sequence can't be read (end of current file)

    # try the next orthogroup in ogs
    if not len(ogs):
        return None

    ogfilename = ogs.pop()
    if not ogfilename:
        # if all OGs are done, return None
        return None

    # if an OG was found open and read a sequence
    fasta = Fasta(filename=ogfilename, mode='r')
    og = basename(ogfilename)
    og = og[:og.rindex('.')]
    # next(fasta)

    # return fasta.next()
    return fasta


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

    ips_submitted = []
    ips_finished = []

    done = False
    fasta = Fasta()
    while not done:
        query = next_og_sequence(ogs, fasta)
        print(query)

    # # open each orthogroup file and read the sequences
    # for ogfilename in ogs:
    #     fasta = Fasta(filename=ogfilename, mode='r')
    #     og = basename(ogfilename)
    #     og = og[:og.rindex('.')]
    #
    #     while fasta.next():
    #         fasta.seq = fasta.seq.rstrip('*')
    #         print(fasta.format(linelen=60))
    #         print()
    #
    #         # send sequences to interproscan in groups of thirty, waiting for batches to complete
    #         ips = interpro_setup(fasta.id, fasta.format())
    #         print(f'submitting {og}:{fasta.id}')
    #         if not ips.submit(show_query=False):
    #             exit(1)
    #
    #         batch_size = 30
    #         poll_time = 30
    #         poll_count = 0
    #         poll_max = 50
    #         while ips.status() != 'finished':
    #             time.sleep(poll_time)
    #             poll_count += 1
    #             print(f'\t ... polling({poll_count}) = {ips.jobstatus}')
    #             if poll_count > poll_max:
    #                 break
    #
    #         print('collecting result')
    #         ips.result()
    #         print(ips.content)

    # parse results

    exit(0)
