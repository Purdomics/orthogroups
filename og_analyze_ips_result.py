"""=================================================================================================
og_interpro.py does interpro searches on all the sequences in specific orthogroups. This program
gathers the resuls and summarizes by orthogroup

Michael Gribskov     12 March 2024
================================================================================================="""
import argparse
# import time
import datetime
import glob
import json
import os.path
import pickle
import sys
from math import exp, log


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

    def __init__(self, accession, description, query):
        self.accession = accession
        self.description = description
        self.query = query
        self.evalue = 1
        self.querypos = [None, None]
        self.subjpos = [None, None]
        self.bounds = ''
        self.ipnumber = ''
        self.go = []

    def pfam(self, hit):
        """-----------------------------------------------------------------------------------------
        Parse pfam information out of the hit section of the interpro json output (converted to
        python)

        :param hit: dict        interpro match data structure
        :return: None
        -----------------------------------------------------------------------------------------"""
        # TODO what about repeated motifs?
        self.evalue = hit['locations'][0]['evalue']
        self.querypos = [hit['locations'][0]['start'], hit['locations'][0]['end']]
        self.subjpos = [hit['locations'][0]['hmmStart'], hit['locations'][0]['hmmEnd']]
        self.bounds = hit['locations'][0]['hmmBounds']

        return None

    def cdd(self, hit):
        """-----------------------------------------------------------------------------------------
        Parse cddinformation out of the hit section of the interpro json output (converted to
        python)

        :param hit: dict        interpro match data structure
        :return: None
        -----------------------------------------------------------------------------------------"""
        self.evalue = hit['locations'][0]['evalue']
        self.querypos = [hit['locations'][0]['start'], hit['locations'][0]['end']]
        self.subjpos = [None, None]
        self.bounds = None

        return None

    def prosite_pattern(self, hit):
        """-----------------------------------------------------------------------------------------
        Parse prosite patterns information out of the hit section of the interpro json output
        (converted to python)

        level of the prosite hit is returned as self.bounds

        :param hit: dict        interpro match data structure
        :return: None
        -----------------------------------------------------------------------------------------"""
        self.evalue = 0
        self.querypos = [hit['locations'][0]['start'], hit['locations'][0]['end']]
        self.subjpos = [None, None]
        self.bounds = hit['locations'][0]['level']

        return None


# ##################################################################################################
# End of class Match
# ##################################################################################################

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
    this_og = None
    for f in glob.glob(filespec):
        dir, fname = os.path.split(f)
        under = fname.index('_')
        og = fname[:under]
        if not og in ogfiles:
            # assumes that files will be ordered by OG, this is true if written by og_interpro.py
            this_og = Og(og)
            ogfiles.add(og)
            oglist.append(this_og)

        this_og.members.append(f)

    return oglist


def get_matches(ip):
    """---------------------------------------------------------------------------------------------
    Extract the match information from the interproscan result. Each query can have multiple
    matches

    query name: ip_parsed['results'][0]['xref'][0]['id']
    match list: ip_parsed['results'][0]['match]
    accession: ip_parsed['results'][0]['match][i]['signature']['accession']
    library type: hit['signature']['signatureLibraryRelease']['library']

    :param ip: 
    :return: 
    ---------------------------------------------------------------------------------------------"""
    all = []
    ip_parsed = json.loads(ip.content)
    # query = ip_parsed['results'][0]['xref'][0]['id']
    for hit in ip_parsed['results'][0]['matches']:
        # for each match found in this result
        accession = hit['signature']['accession']
        description = f"{hit['signature']['name']}"
        query = ip_parsed['results'][0]['xref'][0]['id']
        if hit['signature']['description'] != 'None':
            description += f" - {hit['signature']['description']}"

        libtype = hit['signature']['signatureLibraryRelease']['library']
        if libtype in ('PFAM', 'NCBIFAM', 'PIRSF', 'SMART'):
            m = Match(accession, description, query)
            all.append(m)
            m.pfam(hit)

        elif libtype == 'CDD':
            m = Match(accession, description, query)
            all.append(m)
            m.cdd(hit)

        elif libtype == 'PROSITE_PATTERNS':
            m = Match(accession, description, query)
            all.append(m)
            m.prosite_pattern(hit)

        else:
            sys.stderr.write(f'get_matches() - skipping motif library ({libtype}\n')
            sys.stderr.write(f'Tryng to parse with generic parser\n\n')
            m = Match(accession, description, query)
            all.append(m)
            m.pfam(hit)

        try:
            m.ipr_number = hit['signature']['entry']['accession']
        except TypeError:
            # no ipr accession number
            m.ipr_number = 'None'

        # GO terms
        if hit['signature']['entry']:
            for go in hit['signature']['entry']['goXRefs']:
                m.go.append({'id': go['id'], 'category': go['category'], 'description': go['name']})

    return all

def summarize(matches):
    """---------------------------------------------------------------------------------------------
    Summarize the matches to the orthogroup by annotated feature. One could call this on any list of
    objects, but logically it only makes sens to do it with all the matches in the OG in matches

    :param matches: list of Match
    :return:
    ---------------------------------------------------------------------------------------------"""
    accession = {}
    sequences = set()
    for m in matches:
        sequences.add(m.query)
        score = m.evalue
        if m.accession in accession:
            acc = accession[m.accession]
            acc['query'].append(m.query)
            acc['high'] = max(acc['high'], score)
            acc['low'] = min(acc['low'], score)
            try:
                ls = log(score)
            except Exception as e:
                ls = 0

            acc['sum'] += ls
        else:
            accession[m.accession] = {'query':[m.query], 'description':m.description,
                                      'high':score, 'low':score, 'sum':ls}

    for a in accession:
        acc = accession[a]
        acc['sum'] /= len(acc['query'])
        acc['sum'] = exp(acc['sum'])

        sys.stdout.write(f"\n\tFeature: {a}: {accession[a]['description']}\n")
        sys.stdout.write(f"\tNumber of hits:\t{len(accession[a]['query'])}\n")
        sys.stdout.write(f"\tLow score:\t{accession[a]['low']}\n")
        sys.stdout.write(f"\tHigh score:\t{accession[a]['high']}\n")
        sys.stdout.write(f"\tMean score:\t{accession[a]['sum']:.2g}\n")

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
    for og in oglist:
        matches = []
        interpro = None
        og_member_n = 0
        for member in og.members:
            og_member_n += 1
            try:
                fh = open(member, 'rb')

                # the pickled data is an Interpro object
                interpro = pickle.load(fh)
                # print(interpro.content)

            except Exception as e:
                print(e)
                sys.stderr.write(f'Error reading pickle:{e} - ({member})')

            if interpro:
                matches += get_matches(interpro)

        summarize( matches)
        match_n = len(matches)

    exit(0)
