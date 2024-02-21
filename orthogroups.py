import sys


class Orthogroup:
    """=============================================================================================
    Holds list of orthogroups, sequence files, and sequences in each group

    status:
        0 - success
        1 - input file could not be opened

    Michael Gribskov     09 February 2024
    ============================================================================================="""

    def __init__(self, filename='Orthogroups.tsv'):
        """-----------------------------------------------------------------------------------------
        attributes:
            fh          filehandle with orthogroups data
            proteome    list, proteome files used in orthofineder run
            group       list, index is the orthogroup number, value is a list of species
                            each species is a list of sequences present in the group
                            group[og][species][sequence]

        -----------------------------------------------------------------------------------------"""
        self.fh = None
        self.proteome = []
        self.group = []

        if filename:
            self.fh = self.open_safe(filename)

    def open_safe(self, filename, mode='r', status=1):
        """-----------------------------------------------------------------------------------------
        open file with error check, exit(status) if open fails.
        currently status==1 input file failure, status==2 output file failure

        :param filename: string     path to file
        :param mode: string         mode to open file
        :param status: int          exit status to use on failure
        :return: filehandle         open filehandle
        -----------------------------------------------------------------------------------------"""
        try:
            fh = open(filename, mode)
        except IOError:
            sys.stderr.write(f'Orthogroup.open_safe() - cannot open file ({filename})\n')
            exit(status)

        return fh

    def proteome_read(self):
        """-----------------------------------------------------------------------------------------
        Read the first line of data and save as list of proteome source files

        Format
        Orthogroup	Aurpu_var_nam1_GeneCatalog_proteins_20130418.aa	Cap6580_1_GeneCatalog_proteins_20221012.aa	...

        :return: int    number of proteome sequence files read
        -----------------------------------------------------------------------------------------"""
        line = self.fh.readline().rstrip()
        field = line.split('\t')
        for i in range(1, len(field)):
            seqfile = field[i]
            self.proteome.append(seqfile)

        return len(self.proteome)

    def groups_read(self):
        """-----------------------------------------------------------------------------------------
        Read in the orthogroups. Store in self.group where the index is the orthogroup number, and
        the value is a list sequences in it

        Format
        OG0000000		jgi|Cap6580_1|1008530|estExt_fgenesh1_pg.C_2540015, jgi|Cap6580_1|155246|CE155245_972 ...
        each row has a tab delimited list of species, and within each species a comma delimited list
            of sequences, the list of species follows the proteome order, but the list ends at the
            last species with sequences present

        :return: int    number of sequences read
        -----------------------------------------------------------------------------------------"""
        fh = self.fh
        group = self.group

        n_seq = 0
        for line in fh:
            field = line.rstrip().split('\t')
            row = [[] for _ in range(len(self.proteome))]
            group.append(row)
            s = 0
            for species in field[1:]:
                seq = species.replace(', ', ',').split(',')
                if seq[0] != '':
                    # species is present in this orthogroup
                    row[s] = seq

                s += 1
                n_seq += len(seq)

        return n_seq


# --------------------------------------------------------------------------------------------------
# testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    # basic setup, create object and read proteome (species) and orthogroup data
    og = Orthogroup('data/Orthogroups.tsv')
    n_proteome = og.proteome_read()
    print(f'{n_proteome} sequences file names read')

    n_seq = og.groups_read()
    print(f'{n_seq} sequences read')
    for g in range(5770, 5775):
        print(f'group {g}')
        for s in og.group[g]:
            print(f'\t{s}')

    # count the number of orthogroups that each species is in
    counts = [ 0 for _ in range(n_proteome)]
    for g in og.group:
        for s in range(n_proteome):
            if len(g[s]) > 0:
                counts[s] += 1

    print('\northogroups per species')
    for s in range(n_proteome):
        print(f'{counts[s]}\t{og.proteome[s]}')

    # count the number of sequences per species
    counts = [0 for _ in range(n_proteome)]
    for g in og.group:
        for s in range(n_proteome):
            counts[s] += len(g[s])

    print('\nsequences per species')
    for s in range(n_proteome):
        print(f'{counts[s]}\t{og.proteome[s]}')

    exit(0)
