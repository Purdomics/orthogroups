"""=================================================================================================
Run interproscan for multiple sequences using the EMBL interproscan server

Usage
    intropro_batch.py <Fasta file>

26 December 2018    Michael Gribskov
================================================================================================="""
import sys
import time
import argparse
from sequence.fasta import Fasta
from api.interpro.interpro import Interpro

import textwrap as _textwrap


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter):
    '''=============================================================================================
    Custom formatter for command line argument help
    ============================================================================================='''

    def _split_lines(self, text, width=60):
        '''-----------------------------------------------------------------------------------------
        Gracefully split lines in command line help. lines are split at 60 characters by default

        :param text: str, text to split
        :param width: int, width to split at
        :return:
        -----------------------------------------------------------------------------------------'''
        text = self._whitespace_matcher.sub(' ', text).strip()
        self._max_help_position = 30
        return _textwrap.wrap(text, width)


# end of class CustomFormatter


def arguments_get():
    '''---------------------------------------------------------------------------------------------
    Set up command line arguments, read from command line, and store in argparse.ArgumentParser
    object, cl

    :return: argparse.ArgumentParser object
    ---------------------------------------------------------------------------------------------'''
    minlen_default = 80
    loglevel_default = 1
    batch_limit_default = 20
    batch_wait_default = 60
    cl = argparse.ArgumentParser(description='Interproscan of ORF sequences',
                                 formatter_class=CustomFormatter)
    cl.add_argument('--logfile', type=argparse.FileType('w'), default=sys.stderr,
                    help='Output file for log information')
    cl.add_argument('-m', '--minlen', type=int, default=minlen_default,
                    help='Minimum length ORF to run')
    cl.add_argument('--batch_limit', type=int, default=batch_limit_default,
                    help='Number of sequences to submit per batch')
    cl.add_argument('--batch_wait', type=int, default=batch_wait_default,
                    help='Seconds to wait between polling batch')
    cl.add_argument('--loglevel', type=int, default=loglevel_default,
                    help='detail for reporting REST queries')
    cl.add_argument('fasta_in', type=argparse.FileType('r'))

    return cl.parse_args()  # parse_args  reads the command line


def poll_all(joblist, poll_time=61, poll_max=50):
    """---------------------------------------------------------------------------------------------
    Poll the jobs in the jobs_pending list until all have finished. Finished can be
        1) reached maximum number of polling attempts
        2) returned a status other than success or waiting
        3) success

    :param joblist: list of interproscan objects that have been submitted
    :param poll_time: int, seconds to wait between polling
    :param poll_max: int, maximum number of times to poll
    :return: int, number of jobs in list
    ---------------------------------------------------------------------------------------------"""

    not_all_finished = True
    n = 0
    while not_all_finished:
        n += 1
        not_all_finished = False
        time.sleep(poll_time)

        for ips in joblist:
            if ips.status() == 'FINISHED':
                joblist[ips] = 'finished'

            else:
                not_all_finished = True

    return n


def reformat(job):
    """---------------------------------------------------------------------------------------------
    Return a text string with the result of an interproscan job processed with ips.parse_json()
    An example of a callback function for parsing output

    :param job: interpro object, should be a finished job
    :return: string
    ---------------------------------------------------------------------------------------------"""
    str = ''

    parsed = job.parse_json()
    motifs = parsed['motifs']
    go = parsed['go']
    path = parsed['pathway']

    for m in motifs:
        str += '{}\t{}\t{}\n'.format(m['ipr_accession'],
                                     m['src_accession'],
                                     m['description'])
    for g in go:
        str += '{}\t{}\t{}\t{}\n'.format(g, go[g]['name'], go[g]['category'], go[g]['source'])

    for p in path:
        str += '{}\t{}\t{}\n'.format(p, path[p]['name'], path[p]['source'])

    return str


def save_finished(joblist, reformat=None, fh=None, remove=True):
    """---------------------------------------------------------------------------------------------
    Return the output of all finished jobs as a string.
    Reformat is a callback function used to reformat the output.  for instance, interpro.parse_json
    If fh is True, output is written to the filehandle after reformatting.
    If remove is true, jobs are delete from the list after saving

    :param joblist: dict, ips object is key, staus is value
    :param reformat: function, callback function for formatting job result, argument is ips object
    :param fh: filehandle for writable file
    :param remove: boolean, remove finished jobs after saving
    :return: string, text of job content
    ---------------------------------------------------------------------------------------------"""
    delete_list = []
    for job in joblist:
        if joblist[job] != 'finished':
            # skip unfinished jobs
            continue

        joblist[job] = job.result()  # retrieve the completed job
        # text = job.content
        if reformat:
            text = reformat(job)

        if fh:
            if text:
                fh.write('!{} - {}s\n'.format(job.jobname, job.jobid))
                fh.write('{}\n'.format(text))
            else:
                fh.write('!{} - {} no hits\n'.format(job.jobname, job.jobid))

        if remove:
            delete_list.append(job)

    for job in delete_list:
        del joblist[job]

    return text


# ==================================================================================================
# Main
# ==================================================================================================

args = arguments_get()
args.logfile.write('\ninterpro_batch - interproscan of ORF sequences\n')
args.logfile.write('\tinput ORF file: {}\n'.format(args.fasta_in.name))
args.logfile.write('\tminimum ORF length: {}\n\n'.format(args.minlen))

fasta = Fasta(fh=args.fasta_in)

# The job list keeps track of the ips object that have been created and their current status
# the joblist is a dictionary where the ips object is the key and the value is a status string
joblist = {}

# create a template for the jobs.  The template is an interpro object with the metadata added
template = Interpro(loglevel=1)
template.log_fh = args.logfile
template.email = 'gribskov@purdue.edu'
template.application_select(['Pfam', 'Panther', 'SignalP'])
template.output_select = 'json'
template.poll_time = 60
template.poll_max = 100

sequence_limit = 20
batch_limit = 5
n_sequence = 0
nskip = 60
s = 0
while fasta.next():
    while s < nskip:
        s += 1
        fasta.next()

    if fasta.length() < args.minlen: continue  # skip short sequences
    if n_sequence >= sequence_limit: break  # run up to sequence_limit queries
    n_sequence += 1

    # copy the template and add the sequence information
    ips = template.clone()
    ips.sequence = fasta.translate().seq.rstrip('*').format(linelen=60)
    ips.jobname = fasta.id
    ips.title = 'ORF{}'.format(n_sequence)
    joblist[ips] = 'new'

    # submit job
    if ips.run():
        joblist[ips] = 'submitted'

    if n_sequence % batch_limit and n_sequence < sequence_limit:
        # keep submitting until batch_limit is reached
        continue

    # polling - you only reach here if either the batch_limit or n_sequence limit has been reached
    poll_all(joblist, template.poll_time, template.poll_max)
    save_finished(joblist, reformat, sys.stdout, True)

# end of loop over all sequences

args.logfile.close()
exit(0)
