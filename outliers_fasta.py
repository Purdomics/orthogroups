"""=================================================================================================
orthogroups:outliers_fasta.py

Retrieve the fasta sequences for the expanded and contracted orthogroups

21 February 2024     gribskov
================================================================================================="""
import argparse
import datetime
import sys
import json

# import numpy as np
# from scipy.stats.mstats import trimmed_mean, trimmed_std
#
# from orthogroups import Orthogroup


def process_command_line():
    """---------------------------------------------------------------------------------------------

    :return:
    ---------------------------------------------------------------------------------------------"""
    cl = argparse.ArgumentParser(
        description='Retrieve orthogroup outlier FastA sequences',
        formatter_class=lambda prog: argparse.HelpFormatter(prog, width=120, max_help_position=40)
        )
    cl.add_argument('-j', '--json',
                    help='Output orthogroups from outlier.py',
                    type=str,
                    default='None supplied')

    cl.add_argument('-n', '--ntop',
                    help='number of top/bottom groups to report',
                    type=str,
                    default='all')

    cl.add_argument('-p', '--prefix',
                    help='prefix string for output file',
                    type=str,
                    default='outlier')

    return cl.parse_args()


def make_sequence_list(seqdict, group, n_max):
    """---------------------------------------------------------------------------------------------
    creates a list of the sequences needed from each sequence file

    seqdict = { species: {
                    'id':[seqid1, seqid2, ...],
                    'og_num'[og1, og2, ...] } }

    group = [{'og_num':int, 'z':float,
              'Members': [
                            [seqname1, seqname2 ... ],
                            ...
                        ] },
                ... ]
    where Members is a list over the organisms in the analysis

    :param seqdict: dict            list of sequence ids and corresponding ogs for each proteome file
    :param group: list of dict      dict holds the member sequences of each selected orthogroup
    :param n_max: int               maximum number of groups to add to seqdict
    :return: int, int               number of sequences, number of groups
    ---------------------------------------------------------------------------------------------"""
    ngroups = 0
    nseqs = 0
    for og in group:
        ngroups += 1
        og_num = og['og_num']
        for member in og['Members']:
            if not member:
                # skip blank groups
                break

            species = member[0]
            nseq = member[1]
            for i in range(2, nseq + 2):
                # print(f'{species}\t{og_num}\t{member[i]}')
                seqdict[species]['id'].append(member[i])
                seqdict[species]['og_num'].append(og_num)
                nseqs += 1
        if ngroups > n_max:
            break

    return ngroups, nseqs


# ==================================================================================================
# Main
# ==================================================================================================
if __name__ == '__main__':
    daytime = datetime.datetime.now()
    runstart = daytime.strftime('%Y-%m-%d %H:%M:%S')
    opt = process_command_line()
    sys.stderr.write(f'\noutliers.py {runstart}\n')
    sys.stderr.write(f'Output file prefix: {opt.prefix}')

    # json file with top orthogroups
    try:
        fp_json = open(opt.json, 'r')
    except OSError:
        sys.stderr.write(f'outliers_fasta - unable to open JSON file ({opt.json})\n')
        exit(1)

    sys.stderr.write(f'Top orthogroup file: {opt.json}')
    top = json.load(fp_json)
    ntop = len(top['Expanded'])
    nbottom = len(top['Reduced'])
    sys.stderr.write(f'\tgroups:{ntop + nbottom}\n')

    if opt.ntop == 'all':
        top_n = ntop
        bottom_n = nbottom
    else:
        top_n = min(ntop, int(opt.ntop))
        bottom_n = min(nbottom, int(opt.ntop))

    sys.stderr.write(f'Highest groups: {top_n}\n')
    sys.stderr.write(f'Lowest groups: {bottom_n}\n\n')

    # for each sequence, make a list of all the sequences
    # protein sequence file names are in top['source_data']
    # sequences for each orthogroup are in top['Expanded'|'Reduced']['Members']
    sequences = {species: {'id': [], 'og_num': []} for species in top['source_data']}
    nexp, nseqexp = make_sequence_list(sequences, top['Expanded'], top_n)
    nred, nseqred = make_sequence_list(sequences, top['Reduced'], bottom_n)
    sys.stderr.write(f'{nseqexp+nseqred} sequences will be read for {nexp+nred} groups\n')

    # set up output files for sequences
    # names are opt.<prefix>_<e|r><og_num>.fa, e.g. cold_e568.fa
    # files need to opened simultaneously because each sequence has entries in each orthogroup
    files = {}
    ngroup = 0
    for og in top['Expanded']:
        fname = f'opt.prefix + '_r' + og['og_num'] +

    exit(0)
