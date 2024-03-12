"""=================================================================================================
og_interpro.py does interpro searches on all the sequences in specific orthogroups. This program
gathers the resuls and summarizes by orthogroup

Michael Gribskov     12 March 2024
================================================================================================="""
import argparse
import sys
import time
import datetime
import glob
import os.path


class Og:
    """=============================================================================================
    Holds the summary information about the orthogroup
    ============================================================================================="""

    def __init__(self, name='', emin=1e-20):
        """-----------------------------------------------------------------------------------------

        -----------------------------------------------------------------------------------------"""
        self.name = name if name else ''
        self.n = 0
        self.members = []
        self.e_min = emin
        self.matches = []


class Match:
    """=============================================================================================
    One match from interproscan
    ============================================================================================="""

    def __init__(self, json):
        self.content = json


def process_command_line():
    """---------------------------------------------------------------------------------------------
    Read command line options
    :return:
    ---------------------------------------------------------------------------------------------"""
    cl = argparse.ArgumentParser(
        prog='og_analyze_ips_result.py',
        description='Analyze interpro result for orthogroups',
        formatter_class=lambda prog: argparse.HelpFormatter(prog, width=120, max_help_position=40)
        )

    cl.add_argument('inputfilename',
                    help="Intput file name (default: %(default)s)",
                    type=str,
                    nargs='?',
                    default='pklfiles/*.pkl')  # positional argument

    cl.add_argument('outputfilename',
                    help="Output file name (default: %(default)s)",
                    type=str,
                    nargs='?',
                    default='og.analysis.txt')  # positional argument

    return cl.parse_args()


def expand_input(filespec):
    """---------------------------------------------------------------------------------------------
    Exapand the input file spec to an explicit list of files and return as a dictionary of lists
    with the OG name as the key and a list of pickled files as the value

    :param filespec: string     Possible wildcard file path specification
    :return: list og Og         list of Og objects containing a list of files in each OG
    ---------------------------------------------------------------------------------------------"""
    ogfiles = set()
    oglist = []
    for f in glob.glob(filespec):
        dir, fname = os.path.split(f)
        under = fname.index('_')
        og = fname[:under]
        if not og in ogfiles:
            # assumes that files will be ordered by OG, this is true if written by og_interpro.py
            this_og = Og(og)
            ogfiles.add(og)
            oglist.append(this_og)

        this_og.members.append(fname)

    return oglist


# --------------------------------------------------------------------------------------------------
# Main program
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    runstart = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    opt = process_command_line()

    sys.stderr.write(f'\nog_analyze_ips_result.py {runstart}\n')
    sys.stderr.write(f'\tOG input: {opt.inputfilename}\n')
    sys.stderr.write(f'\tOutput directory: {opt.outputfilename}\n')

    oglist = expand_input(opt.inputfilename)
    for ogname in oglist:
        og = Og(ogname)
        og.members = oglist[ogname]

    exit(0)
