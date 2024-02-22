"""=================================================================================================
orthogroups:outliers.py

for the orthogroups data, find  which groups are expanded or contracted
method:
    normalize counts for all groups (trimmed standard normal deviate)
    for target group, find groups low/high deviation

The count data is a list of lists or lists: group refers to Orthogroup.group
group[i] gives the list of matching sequences for group[i]
group[i][j] give the list of matching sequences in each proteome (may be empty list)
    the groups are in the same order as Orthogroup.proteome
    the length of group[i][j] is the count of matching sequences in each proteome
group[i][j][k] are the individual matching sequences

18 February 2024     gribskov
================================================================================================="""
import argparse
import datetime
import json
import sys

import numpy as np
from scipy.stats.mstats import trimmed_mean, trimmed_std

from json_custom_encoder import CompactJSONEncoder
from orthogroups import Orthogroup


def process_command_line():
    """---------------------------------------------------------------------------------------------
    TODO add defaults to output
    :return:
    ---------------------------------------------------------------------------------------------"""
    cl = argparse.ArgumentParser(
        description='Find expanded and contracted orthogroups',
        formatter_class=lambda prog: argparse.HelpFormatter(prog, width=120, max_help_position=40)
        )
    cl.add_argument('-g', '--orthogroup',
                    help='Orthogroups.tsv file',
                    type=str,
                    default='Orthogroups.tsv')

    cl.add_argument('-n', '--ntop',
                    help='number of top/bottom groups to report',
                    type=int,
                    default=20)

    cl.add_argument('-f', '--fraction',
                    help='total fraction of observations to trim in mean/std.dev. calculation',
                    type=int,
                    default=0.5)

    cl.add_argument('-t', '--target',
                    help='comma delimited string with list of target organisms',
                    type=str,
                    default='')

    cl.add_argument('-j', '--json',
                    help='JSON output file',
                    type=str,
                    default='')

    cl.add_argument('-v', '--tsv',
                    help='file for normalized counts in TSV format',
                    type=str,
                    default='')

    return cl.parse_args()


def trimmed_stats(og, proportion, ignore_zero=True):
    """---------------------------------------------------------------------------------------------
    return the trimmed mean and standard deviation. proportion defines how much data is omitted on
    each side (proportion/2). After sorting the first index used is the first >= proportion/2 and
    the last is the one <= 1.0 - proportion/2

    :param og: list             Othogroup object, rows are orthogroups
    :param proportion:          proportion to trim proportion/2 trimmed on each side
    :param ignore_zero: bool    NOT IMPLEMENTED. if True, zeroes are omitted otherwise they are
                                treated as values
    :return: float, float       Z normalized data
    ---------------------------------------------------------------------------------------------"""
    # use_zero = not ignore_zero
    cut = proportion / 2

    # nil = 0
    # if ignore_zero:
    #     nil = np.nan
    data = og.group
    norm = np.full([len(data), len(og.proteome)], np.nan)

    # make count vector from original data, and store in norm using zero or NaN depending
    # on the state of ignore_zero
    ngroup = len(data)
    for g in range(ngroup):
        nspecie = len(data[g])
        norm[g] = [len(data[g][s]) for s in range(nspecie)]

    # if ignore_zero:
    #     norm[norm == 0] = np.nan

    # mean and standard deviation, values are trimmed until index+1 >= n*cut
    # for n=12 and cut = .125, two values are cut off

    ave = trimmed_mean(norm, limits=(cut, cut), inclusive=(False, False), axis=1)
    std = trimmed_std(norm, limits=(cut, cut), inclusive=(False, False), axis=1)
    std[std == 0] = 1

    return ((norm.T - ave) / std).T


def tsv_write(opt, col_labels, z):
    """---------------------------------------------------------------------------------------------
    if opt.tsv, write a tab delimited file with the normalized data

    :param opt: namespace       command line options
    :param col_labels: list     column labels (str)
    :param z: np.array          normalized data

    :return: int                number of line sritten
    ---------------------------------------------------------------------------------------------"""
    if not opt.tsv:
        return 0

    tsvfp = open(opt.tsv, 'w')

    # column lables
    tsvfp.write(f'Orthogroup')
    for label in col_labels:
        tsvfp.write(f'\t{label}')
    tsvfp.write('\n')

    # data matrix
    # TODO may not be writing the last column
    nline = 0
    for row in z:
        tsvfp.write(f'OG{nline:07d}')
        for col in range(1, len(row)):
            tsvfp.write(f'\t{row[col]}')

        tsvfp.write('\n')
        nline += 1

    tsvfp.close()
    return nline


# ==================================================================================================
# Main
# ==================================================================================================
if __name__ == '__main__':
    daytime = datetime.datetime.now()
    runstart = daytime.strftime('%Y-%m-%d %H:%M:%S')
    opt = process_command_line()

    sys.stderr.write(f'\noutliers.py {runstart}\n')
    sys.stderr.write(f'Trim fraction: {opt.fraction}\n')
    sys.stderr.write(f'Top groups: {opt.ntop}\n')
    sys.stderr.write(f'Orthogroups: {opt.orthogroup}\n')
    jsonfile = None
    if opt.json:
        sys.stderr.write(f'JSON Output: {opt.json}\n')
        jsonfile = open(opt.json, 'w')
    if opt.tsv:
        sys.stderr.write(f'TSV output: {opt.tsv}\n')

    target = opt.target.split(',')
    sys.stderr.write(f'\nTargets: {opt.target}\n\n')

    og = Orthogroup(opt.orthogroup)
    n_proteome = og.proteome_read()
    sys.stderr.write(f'{n_proteome} sequence file names read\n')
    n_seq = og.groups_read()
    sys.stderr.write(f'{n_seq} orthogroup sequences read\n\n')

    # make a trimmed list of genome names (part up to first underline)
    organism = {}
    orgidx = []
    for p in range(len(og.proteome)):
        field = og.proteome[p].split('_')
        organism[field[0]] = p
        orgidx.append(field[0])

    # normalize the count data so that each row is a standard normal deviate, using trimmed mean
    # and STD dev.

    cols = []
    for t in target:
        try:
            cols.append(organism[t])
        except KeyError:
            sys.stderr.write(f'Unknown organism ({t})')
            exit(2)

    # summary of source sequences for looking up fasta
    json_out = {'source_data': {}}
    print(f'Source Data:')
    for p in range(len(og.proteome)):
        print(f'{orgidx[p]}\t{og.proteome[p]}')
        json_out['source_data'][orgidx[p]] = og.proteome[p]

    print(f'\nExpanded/contracted Groups:')
    z = trimmed_stats(og, opt.fraction)
    selected = z[:, cols].mean(axis=1)
    tsv_write(opt, orgidx, z)

    n = 0
    ranked = sorted(range(len(selected)), key=lambda g: selected[g])

    print(f'\n{opt.ntop} most reduced in {target}')
    reduced = []
    for g in ranked[:opt.ntop]:
        print(f'\nOrthogroup {g:6d}\t{selected[g]:.3f}')
        this_group = {'og_num': g, 'z': selected[g], 'Members': []}
        reduced.append(this_group)
        members = this_group['Members']
        for i in range(n_proteome):
            print(f'\t{orgidx[i]}: {len(og.group[g][i])}\t{og.group[g][i]}')
            members.append([orgidx[i], len(og.group[g][i])] + og.group[g][i])
        json_out['Reduced'] = reduced

    print(f'\n{opt.ntop} most expanded in {target}')
    expanded = []
    for g in ranked[-opt.ntop:]:
        this_group = {'og_num': g, 'z': selected[g], 'Members': []}
        expanded.append(this_group)
        members = this_group['Members']
        print(f'\nOrthogroup {g:6d}\t{selected[g]:.3f}')
        for i in range(n_proteome):
            print(f'\t{orgidx[i]}: {len(og.group[g][i])}\t{og.group[g][i]}')
            members.append([orgidx[i], len(og.group[g][i])] + og.group[g][i])
            if not og.group[g][i]:
                members[-1].append('')

    json_out['Expanded'] = expanded
    if jsonfile:
        jsonfile.write(f'{json.dumps(json_out, indent=2, cls=CompactJSONEncoder)}\n')
        jsonfile.close()

    exit(0)
