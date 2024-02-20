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
# from math import floor, ceil
from scipy.stats.mstats import trimmed_mean, trimmed_std
# from scipy import stats
import numpy as np
from orthogroups import Orthogroup


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
    organism = {}
    for p in range(len(og.proteome)):
        field = og.proteome[p].split('_')
        organism[field[0]] = p

    # normalize the count data so that each row is a standard normal deviate, using trimmed mean
    # and STD dev.

    z = trimmed_stats(og, 0.5, ignore_zero=True)

    exit(0)
