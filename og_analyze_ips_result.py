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
        self.evalue = 1
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
    try:
        ip_parsed = json.loads(ip.content)
    except Exception as e:
        # sys.stderr.write( f'\n{e} - json error\n')
        return []

    # query = ip_parsed['results'][0]['xref'][0]['id']
    for hit in ip_parsed['results'][0]['matches']:
        # for each match found in this result
        accession = hit['signature']['accession']
        description = f"{hit['signature']['name']}"
        query = ip_parsed['results'][0]['xref'][0]['id']
        if hit['signature']['description']:
            # discription is None
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
            m.ipr_description = hit['signature']['entry']['description']
        except TypeError:
            # no ipr accession number
            m.ipr_number = 'None'
            m.ipr_description = ''

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
    species = set()
    go = {}

    # short name for gene ontology
    gshort = {'BIOLOGICAL_PROCESS': '(BP)',
              'MOLECULAR_FUNCTION': '(MF)',
              'CELLULAR_COMPONENT': '(CC)'}

    for m in matches:
        sequences.add(m.query)
        score = m.evalue
        try:
            # avoid log zero
            ls = log(score)
        except Exception as e:
            ls = 0
        if m.accession in accession:
            acc = accession[m.accession]
            acc['query'].append(m.query)
            acc['high'] = max(acc['high'], score)
            acc['low'] = min(acc['low'], score)
            acc['sum'] += ls
        else:
            accession[m.accession] = {'query': [m.query], 'description': m.description,
                                      'ipr':   f"{m.ipr_number} - {m.ipr_description}",
                                      'high':  score, 'low': score, 'sum': ls, 'go':{}}
        acc = accession[m.accession]
        for g in m.go:
            gid = g['id']
            if gid in acc['go']:
                acc['go'][gid]['n'] += 1
            else:
                acc['go'][gid] = {'n':           1,
                           'id':          g['id'],
                           'category':    gshort[g['category']],
                           'description': g['description']}

    for a in accession:
        acc = accession[a]
        acc['sum'] /= len(acc['query'])
        acc['sum'] = exp(acc['sum'])

        sys.stdout.write(f"\n\tFeature: {a}: {accession[a]['description']} ")
        sys.stdout.write(f"(IPR: {accession[a]['ipr']})\n")
        sys.stdout.write(f"\tNumber: {len(accession[a]['query'])}\t")
        sys.stdout.write(f"High score: {accession[a]['high']}\t")
        sys.stdout.write(f"Mean score: {accession[a]['sum']:.2g}\t")
        sys.stdout.write(f"Low score: {accession[a]['low']}\n")
        go = accession[a]['go']
        for gid in sorted(go, key=lambda g: go[g]['category']):
            goterm = go[gid]
            sys.stdout.write(f"\tNumber: {goterm['n']}\t{goterm['id']}: {goterm['category']} ")
            sys.stdout.write(f"{goterm['description']}\n")

    return None


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
        sys.stdout.write(f'\n{og.name}\t{len(og.members)} sequences')
        matches = []
        interpro = None
        og_member_n = 0
        for member in og.members:
            og_member_n += 1
            try:
                fh = open(member, 'rb')

                # the pickled data is an Interpro object
                interpro = pickle.load(fh)

            except Exception as e:
                sys.stderr.write(f'\nError reading pickle:{e} - ({member})\n\n')

            if interpro:
                match = get_matches(interpro)
                if match:
                    matches += match
                else:
                    print(f'json unreadable {member}')

        summarize(matches)
        match_n = len(matches)

    exit(0)
