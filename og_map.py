"""=================================================================================================
og_map

map orthogroups from to runs of orthofinder to each other

Michael Gribskov     01 April 2024
================================================================================================="""

# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    old_og_file = 'feb22_orthogroups.tsv'
    new_og_file = "mar31_orthogroups.tsv"

    # read the old orthogroups and make an index where the key is a sequence name and the value is
    # the old orthogroup
    old = open(old_og_file, 'r')
    header = old.readline()

    old_idx = {}
    og_n = 0
    seq_n = 0
    for line in old:
        og_n += 1
        # print(line)
        field = line.rstrip().split('\t')

        og = field[0]
        for sequence_list in field[1:]:
            if not sequence_list:
                continue
            for sequence in sequence_list.split(','):
                seq_n += 1
                old_idx[sequence] = og

    print(f'{seq_n} sequences read from {og_n} orthogroups')
    old.close()

    # for each new orthogroup, read in the list of sequences and see what their former group was

    new = open(new_og_file, 'r')
    header = new.readline()

    og_n = 0
    seq_n = 0
    og_count = {}
    for line in new:
        og_n += 1
        field = line.rstrip().split('\t')

        og = field[0]
        og_count[og] = {'all': 0, 'unknown': 0}

        for sequence_list in field[1:]:
            if not sequence_list:
                continue
            for sequence in sequence_list.split(','):
                seq_n += 1
                og_count[og]['all'] += 1
                try:
                    og_old = old_idx[sequence]
                    if og_old in og_count[og]:
                        og_count[og][og_old] += 1
                    else:
                        og_count[og][og_old] = 1
                except KeyError:
                    og_count[og]['unknown'] += 1
    new.close()

    # result
    out = open('og.map.out', 'w')
    out.write(f'new_og:count\told_og:count(pct) ...\n')
    warning = {}
    warning_threshold = 0.90
    for og in og_count:
        matches = sorted(og_count[og], key=lambda g: og_count[og][g], reverse=True)
        best = matches[1]
        n = og_count[og]['all']

        best_pct = og_count[og][matches[1]] / n
        if best_pct < warning_threshold:
            warning[og] = og_count[og]

        out.write(f"{og}:{n}")
        for other_og in matches[1:]:
            pct = og_count[og][other_og] / n
            out.write(f"\t{other_og}:{og_count[og][other_og]}({pct * 100:.1f}%)")

        out.write('\n')

    # OGs with biggest change
    print(f'{len(warning)} conversion warnings, change > {warning_threshold * 100:.1f}%)')
    for og in warning:
        matches = sorted(og_count[og], key=lambda g: og_count[og][g], reverse=True)
        best = matches[1]
        n = og_count[og]['all']

        best_pct = og_count[og][matches[1]] / n
        if best_pct < warning_threshold:
            warning[og] = og_count[og]

        print(f"{og}:{n}", end='')
        for other_og in matches[1:]:
            pct = og_count[og][other_og] / n
            print(f"\t{other_og}:{og_count[og][other_og]}({pct * 100:.1f}%)", end='')

        print()

    exit(0)
