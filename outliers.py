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
import sys
from math import floor, ceil
# from scipy.stats import trim_mean, tstd
import numpy as np
import statistics as stat
from orthogroups import Orthogroup


def trimmed_stats(data, norm, proportion, ignore_zero=True):
    """---------------------------------------------------------------------------------------------
    return the trimmed mean and standard deviation. proportion defines how much data is omitted on
    each side (proportion/2). After sorting the first index used is the first >= proportion/2 and
    the last is the one <= 1.0 - proportion/2

    :param data: list           Othogroup.group data, rows are orthogroups
    :param norm:
    :param proportion:          proportion to trim proportion/2 trimmed on each side
    :param ignore_zero: bool    if True, zeroes are omitted otherwise they are treated as values
                                note that the list
    :return: float, float       trimmed mean and standard deviation
    ---------------------------------------------------------------------------------------------"""
    use_zero = not ignore_zero
    cutmin = proportion / 2
    cutmax = 1.0 - proportion / 2

    nil = 0
    if ignore_zero:
        nil = np.nan

    # make count vector from original data, and store in norm using zero or NaN depending
    # on the state of ignore_zero
    ngroup = len(data)
    for g in range(ngroup):
        nspecie = len(data[g])
        norm[g] = [len(data[g][s]) for s in range(nspecie)]

    if ignore_zero:
        norm[norm == 0] = np.nan

    # mean and standard deviation ignoring NaN
    ave = np.nanmean(norm, 1)
    std = np.nanstd(norm, 1)
    std[std==0] = 1

    for g in range(ngroup):
        aave = ave[g]
        astd = std[g]
        norm[g] = (norm[g] - aave) /astd

    # start = floor(n * cutmin)
    # stop = ceil(n * cutmax)
    # print(f'n: {n}\tstart: {start}\tstop: {stop}')



    return ave, std


# ==================================================================================================
# Main
# ==================================================================================================
if __name__ == '__main__':
    target = ['Cryan3', 'Crymi1']
    ogdata = 'data/Orthogroups.tsv'
    og = Orthogroup(ogdata)
    n_proteome = og.proteome_read()
    print(f'{n_proteome} sequences file names read')

    n_seq = og.groups_read()
    print(f'{n_seq} sequences read')

    # make a trimmed list of genome names (part up to first underline)
    organism = []
    for p in og.proteome:
        field = p.split('_')
        organism.append(field[0])

    # normalize the count data so that each row is a standard normal deviate, using trimmed mean
    # and STD dev.
    # transformed = np.zeros(shape=(len(og.group), len(og.proteome)))
    transformed = np.full([len(og.group), len(og.proteome)], np.nan)
    stats = trimmed_stats(og.group, transformed, 0.25, ignore_zero=True)

    exit(0)
