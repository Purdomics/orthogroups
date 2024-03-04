import sys
import time
import json
import requests
from api.jobmanager_api import JobManagerAPI


class Interpro(JobManagerAPI):
    """=============================================================================================
    Interpro class for running interproscan.  See the end of the file for usage examples.
    Can be used by itself or from jobmanager.py

    25 December 2018    Michael Gribskov
    ============================================================================================="""
    # availble options/parameters taken from
    # https://www.ebi.ac.uk/Tools/services/rest/iprscan5/parameterdetails/appl

    # available outputs
    # from https://www.ebi.ac.uk/Tools/services/rest/iprscan5/resulttypes
    # log - The output from the tool itself
    # out - The results of the job (XML format)
    # tsv - The results of the job in text format, tab separated values
    # xml - The results of the job in XML
    # gff - The results of the job in GFF3 format
    # json - The results of the job in JSON format
    # htmltarball - The results of the job in a tarball zip file
    # sequence - Input sequence as seen by the tool
    # submission - The submission details which were submitted as a job
    # available = {'applications':['TIGRFAM', 'SFLD', 'Phobius', 'SignalP', 'SignalP_EUK',
    #                              'SignalP_GRAM_POSITIVE', 'SignalP_GRAM_NEGATIVE', 'SUPERFAMILY',
    #                              'Panther', 'Gene3d', 'HAMAP', 'PrositeProfiles',
    #                              'PrositePatterns', 'Coils', 'SMART', 'CDD', 'PRINTS', 'PfamA',
    #                              'MobiDBLite', 'PIRSF', 'TMHMM', ],
    available = {'applications': ['NCBIfam', 'SFLD', 'Phobius', 'SignalP', 'SignalP_EUK',
                                  'SignalP_GRAM_POSITIVE', 'SignalP_GRAM_NEGATIVE', 'SUPERFAMILY',
                                  'Panther', 'Gene3d', 'HAMAP', 'PrositeProfiles',
                                  'PrositePatterns', 'Coils', 'SMART', 'CDD', 'PRINTS', 'PfamA',
                                  'MobiDBLite', 'PIRSF', 'TMHMM','FunFam', 'PIRSR'],
                 'commands':    ['run', 'status', 'result'],
                 'outputs':     ['out', 'log', 'tsv', 'xml', 'gff', 'json',
                                 'htmltarball', 'sequence', 'submission']
                 }

    def __init__(self):
        """-----------------------------------------------------------------------------------------
        interpro query/response constructor

        self.message is inherited from JobManagerAPI

        loglevel   0 no log, 1 job submission/completion, 2 all
        -----------------------------------------------------------------------------------------"""

        self.email = ''  # user email (optional)
        self.title = ''  # title for job (optional)
        self.sequence = ''
        self.applications = []
        self.output = 'json'
        self.parameters = {}

        self.url = u'https://www.ebi.ac.uk/Tools/services/rest/iprscan5/'
        self.jobid = ''
        self.jobname = ''
        self.jobstatus = ''

        self.response = None
        self.content = ''

    def application_select(self, selected, keep=False):
        """-----------------------------------------------------------------------------------------
        Add a list of applications to be run.  Each application in the list is compared to the
        available applications and if not present a warning is issued.  The default is to run all
        applications so an empty selected list signifies the default.

        :param selected: list of strings, selected applications to run
        :param keep: Boolean, retain current applications, just add new ones
        :return: int, number of selected applications
        -----------------------------------------------------------------------------------------"""
        if not keep:
            self.applications = []

        for app in selected:
            if app == 'Pfam':
                app = 'PfamA'
            if app in Interpro.available['applications']:
                self.applications.append(app)
            else:
                self.message = {'type':    'not_available',
                                'text':    f'application={app}',
                                'loglevel':2}

        return len(self.applications)

    def output_select(self, selected):
        """-----------------------------------------------------------------------------------------
        select the output format.  Only one is allowed

        :param selected: string, one of the formats in self.output_avail
        :return: True if format is available
        -----------------------------------------------------------------------------------------"""
        # self.output = ''
        if selected in Interpro.available['outputs']:
            self.output = selected
        else:
            self.message = {'type':    'not_available',
                            'text':    f'output={selected}',
                            'loglevel':2}

            return False

        return True

    def parameter_select(self, select):
        """-----------------------------------------------------------------------------------------
        Select additional tag value pairs to add to parameters.  There is no checking so be correct

        :param select: dict
        :return: int number of parameters in dictionary
        -----------------------------------------------------------------------------------------"""
        for key in select:
            self.parameters[key] = select[key]

        return len(self.parameters)

    def parse_json(self):
        """-----------------------------------------------------------------------------------------
        Parse the contented returned from the server in JSON format. Three outputs are produced
            A list of dictionaries for each hit in the sequence
            A list of dictionaries listing  GO terms and what entries they were drawn from
                keys: 'name': gene ontology ID
                      'category': ontology category = BIOLOGICAL_PROCESS |
                                                      MOLECULAR_FUNCTION |
                                                      CELLULAR_COMPONENT,
                      'source': list of strings, motifs that are associated with this term
            A list of dictionaries listing pathways and the entry theies were drawn from
                keys: 'databaseName': the pathway database KEGG | Reactome | Metacyc
                      'id': UID of the pathway in its database
                      'name': text description of pathway
                      'source': list of strings, motifs that are associated with this term

        This is fairly specific for my purpose

        :return:
        -----------------------------------------------------------------------------------------"""
        pjson = json.loads(self.content)

        matches = pjson['results'][0]['matches']
        go_all = {}
        path_all = {}
        motifs = []
        for m in matches:
            # each match in the interproscan search
            # source_accession is the UID of the matching motif in the source database
            # source is the database and version listed at interpro
            # entry is the subtree of information about a match
            entry = m['signature']['entry']
            signature = m['signature']['signatureLibraryRelease']
            source_accession = m['signature']['accession']
            source = '{} {}'.format(signature['library'],
                                    signature['version'])

            if not entry:
                # panther subfamily entries have no ['signature']['entry']
                # e.g., PTHR21139:SF24
                continue

            # parse an entry, en entry is a hit vs a specific entry in a database
            motifs.append({'ipr_accession':entry['accession'],
                           'src_accession':source_accession,
                           'description':  entry['description'] or ''})

            # name = entry['name']
            # type = entry['type']

            if 'goXRefs' in entry:
                gostr = ''
                for go in entry['goXRefs']:
                    gostr += '{} ({}:{})'.format(go['id'], go['category'], go['name'])
                    if go['id'] in go_all:
                        go_all[go['id']]['source'].append(source_accession)
                    else:
                        go_all[go['id']] = {'name':  go['name'], 'category':go['category'],
                                            'source':[source_accession]}

            if 'pathwayXRefs' in entry:
                pathstr = ''
                for path in entry['pathwayXRefs']:
                    pathstr += '{} ({})'.format(path['id'], path['name'])
                    if path['databaseName'] == 'Reactome':
                        field = path['id'].split('-')
                        id = 'Reactome:{}'.format(field[2])
                    elif path['databaseName'] == 'MetaCyc':
                        id = 'Metacyc:{}'.format(path['id'])
                    elif path['databaseName'] == 'KEGG':
                        id = 'KEGG:{}'.format(path['id'])
                    else:
                        print('unknown pathway {} | {} | {}'.format(path['databaseName'],
                                                                    path['id'], path['name']))

                    if id in path_all:
                        if source_accession not in path_all[id]['source']:
                            path_all[id]['source'].append(source_accession)
                    else:
                        path_all[id] = {'name':  path['name'],
                                        'source':[source_accession]}

        return {'jobname':self.jobname, 'motifs':motifs, 'go':go_all, 'pathway':path_all}

    def submit(self, show_query=False):
        """-----------------------------------------------------------------------------------------
        Construct a REST command and dispatch the job to the server
        Any previously existing jobID is overwritten

        :param show_query: boolean, print query if true
        :return: logical, True = success, False = failure
        -----------------------------------------------------------------------------------------"""
        is_success = False

        # general fields for all queries
        param = {u'email': self.email, u'title':self.title, u'sequence':self.sequence,
                 u'output':self.output}

        if self.applications:
            # add selected applications
            param['appl'] = ','.join(self.applications)

        if self.parameters:
            for para in self.parameters:
                param[para] = self.parameters[para]

        command = self.url + 'run'
        self.response = requests.post(command, files=param, headers={'User-Agent':'ips-client'})

        if show_query:
            # print out query if requested
            print(self.response.request.headers, '\n')
            print(self.response.request.body, '\n')

        if self.response_is_error('submitting job'):
            self.jobstatus = 'failed'
        else:
            # success
            self.jobid = self.response.text
            self.jobstatus = 'submitted'
            self.message = {'type':    'submitted',
                            'text':    f'job_name={self.title};job_id={self.jobid}',
                            'loglevel':1}

            is_success = True

        return is_success

    def status(self, log=True):
        """-----------------------------------------------------------------------------------------
        Poll job status at the server. The job is polled only once so it you want to poll
        multiple times call this method in a loop

        :return: string, status of job at server
        -----------------------------------------------------------------------------------------"""
        command = self.url + 'status/' + self.jobid
        self.response = requests.get(command)
        response_text = self.response.text.rstrip()

        if 'RUNNING' in self.response.text:
            self.jobstatus = 'running'
            self.message = {'type':    'polling',
                            'text':    f'job_id={self.jobid};response={response_text}',
                            'loglevel':2}

        elif 'FINISHED' in self.response.text:
            if self.jobstatus != 'finished':
                # only print finished message once
                self.jobstatus = 'finished'
                self.message = {'type':    'finished',
                                'text':    f'job_id={self.jobid}',
                                'loglevel':1}

        return self.jobstatus

    def result(self):
        """-----------------------------------------------------------------------------------------
        Retrieve the result

        :return: string, 'retrieved' if successful, '' if unsuccessful (False)
        -----------------------------------------------------------------------------------------"""
        # get the final result
        command = self.url + 'result/' + self.jobid + '/' + self.output
        self.response = requests.get(command)
        if not self.response_is_error('retrieving result'):
            # success
            self.content = self.response.text
            self.message = {'type':    'retrieved',
                            'text':    f'job_id={self.jobid};output_len={len(self.output)}',
                            'loglevel':1}
            return 'retrieved'

        return ''

    def response_is_error(self, task):
        """-----------------------------------------------------------------------------------------
        Return true if the response code is other than 200. Write error message to stderr if
        loglevel > 1. Task is a string describing the task that failed for inclusion in the error
        message.  The most recent response is stored in self.response

        :param task: string, text description of response being tested for error message
        :return: logical True = error, False = no error
        -----------------------------------------------------------------------------------------"""
        if self.response.status_code == 200:
            # success
            is_error = False

        else:
            # error
            is_error = True
            self.jobstatus = 'error'
            self.message = {'type':    task,
                            'text':    f'job_id={self.jobid};status={self.response.status_code}',
                            'loglevel':1}

        return is_error


def json_test():
    """---------------------------------------------------------------------------------------------
    return some json for testing.  This is not part of the Interpro class

    :return:
    ---------------------------------------------------------------------------------------------"""
    return '''{
 "interproscan-version": "5.50-84.0",
"results": [ {
  "sequence" : "MAPSRKFFVGGNWKMNGRKQSLGELIGTLNAAKVPADTEVVCAPPTAYIDFARQKLDPKIAVAAQNCYKVTNGAFTGEISPGMIKDCGATWVVLGHSERRHVFGESDELIGQKVAHALAEGLGVIACIGEKLDEREAGITEKVVFEQTKVIADNVKDWSKVVLAYEPVWAIGTGKTATPQQAQEVHEKLRGWLKSNVSDAVAQSTRIIYGGSVTGATCKELASQPDVDGFLVGGASLKPEFVDIINAKQ",
  "md5" : "a8d44fc2c980a7677a3b54788d0fa323",
  "matches" : [ {
    "signature" : {
      "accession" : "PTHR21139",
      "name" : "TRIOSEPHOSPHATE ISOMERASE",
      "description" : null,
      "signatureLibraryRelease" : {
        "library" : "PANTHER",
        "version" : "15.0"
      },
      "entry" : {
        "accession" : "IPR000652",
        "name" : "Triosephosphate_isomerase",
        "description" : "Triosephosphate isomerase",
        "type" : "FAMILY",
        "goXRefs" : [ {
          "name" : "triose-phosphate isomerase activity",
          "databaseName" : "GO",
          "category" : "MOLECULAR_FUNCTION",
          "id" : "GO:0004807"
        } ],
        "pathwayXRefs" : [ {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-DME-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-GGA-70171"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-DME-70263"
        }, {
          "name" : "Inositol phosphate metabolism",
          "databaseName" : "KEGG",
          "id" : "00562"
        }, {
          "name" : "glycolysis IV",
          "databaseName" : "MetaCyc",
          "id" : "PWY-1042"
        }, {
          "name" : "glycolysis II (from fructose 6-phosphate)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-5484"
        }, {
          "name" : "Microbial metabolism in diverse environments",
          "databaseName" : "KEGG",
          "id" : "01120"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-MMU-70263"
        } ]
      }
    },
    "locations" : [ {
      "start" : 1,
      "end" : 249,
      "hmmStart" : 1,
      "hmmEnd" : 249,
      "hmmLength" : 249,
      "hmmBounds" : "COMPLETE",
      "envelopeStart" : 1,
      "envelopeEnd" : 249,
      "location-fragments" : [ {
        "start" : 1,
        "end" : 249,
        "dc-status" : "CONTINUOUS"
      } ]
    } ],
    "evalue" : 0.0,
    "familyName" : "Not available",
    "score" : 585.3,
    "model-ac" : "PTHR21139"
  }, {
    "signature" : {
      "accession" : "PS51440",
      "name" : "TIM_2",
      "description" : "Triosephosphate isomerase (TIM) family profile.",
      "signatureLibraryRelease" : {
        "library" : "PROSITE_PROFILES",
        "version" : "2019_11"
      },
      "entry" : {
        "accession" : "IPR000652",
        "name" : "Triosephosphate_isomerase",
        "description" : "Triosephosphate isomerase",
        "type" : "FAMILY",
        "goXRefs" : [ {
          "name" : "triose-phosphate isomerase activity",
          "databaseName" : "GO",
          "category" : "MOLECULAR_FUNCTION",
          "id" : "GO:0004807"
        } ],
        "pathwayXRefs" : [ {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-DME-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-GGA-70171"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-MMU-70263"
        } ]
      }
    },
    "locations" : [ {
      "start" : 6,
      "end" : 247,
      "score" : 81.047,
      "alignment" : "KFFVGGNWKMNGRKQSLGELIGTLNAaKVPAD--TEVVCAPPTAYIDFARQKLdPKIAVAAQNCYKVTNGAFTGEISPGMIKDCGATWVVLGHSERRHVFGESDELIGQKVAHALAEGLGVIACIGEKLDEREAGITEKVVFEQTKVIADNVKD-WSKVVLAYEPVWAIGTGKTATPQQAQEVHEKLRGWLKSNVSDAVAQSTRIIYGGSVTGATCKELASQPDVDGFLVGGASLK-PEFVDIINA",
      "location-fragments" : [ {
        "start" : 6,
        "end" : 247,
        "dc-status" : "CONTINUOUS"
      } ]
    } ],
    "model-ac" : "PS51440"
  }, {
    "signature" : {
      "accession" : "MF_00147_B",
      "name" : "TIM_B",
      "description" : "Triosephosphate isomerase [tpiA].",
      "signatureLibraryRelease" : {
        "library" : "HAMAP",
        "version" : "2020_05"
      },
      "entry" : {
        "accession" : "IPR022896",
        "name" : "TrioseP_Isoase_bac/euk",
        "description" : "Triosephosphate isomerase, bacterial/eukaryotic",
        "type" : "FAMILY",
        "goXRefs" : [ {
          "name" : "triose-phosphate isomerase activity",
          "databaseName" : "GO",
          "category" : "MOLECULAR_FUNCTION",
          "id" : "GO:0004807"
        }, {
          "name" : "glycolytic process",
          "databaseName" : "GO",
          "category" : "BIOLOGICAL_PROCESS",
          "id" : "GO:0006096"
        } ],
        "pathwayXRefs" : [ {
          "name" : "Microbial metabolism in diverse environments",
          "databaseName" : "KEGG",
          "id" : "01120"
        }, {
          "name" : "glycolysis IV",
          "databaseName" : "MetaCyc",
          "id" : "PWY-1042"
        } ]
      }
    },
    "locations" : [ {
      "start" : 5,
      "end" : 247,
      "score" : 58.311584,
      "alignment" : "Not available",
      "location-fragments" : [ {
        "start" : 5,
        "end" : 247,
        "dc-status" : "CONTINUOUS"
      } ]
    } ],
    "model-ac" : "MF_00147_B"
  }, {
    "signature" : {
      "accession" : "PS00171",
      "name" : "TIM_1",
      "description" : "Triosephosphate isomerase active site.",
      "signatureLibraryRelease" : {
        "library" : "PROSITE_PATTERNS",
        "version" : "2019_11"
      },
      "entry" : {
        "accession" : "IPR020861",
        "name" : "Triosephosphate_isomerase_AS",
        "description" : "Triosephosphate isomerase, active site",
        "type" : "ACTIVE_SITE",
        "goXRefs" : [ {
          "name" : "triose-phosphate isomerase activity",
          "databaseName" : "GO",
          "category" : "MOLECULAR_FUNCTION",
          "id" : "GO:0004807"
        } ],
        "pathwayXRefs" : [ {
          "name" : "gluconeogenesis II (Methanobacterium thermoautotrophicum)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6142"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-DDI-70171"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-MMU-70263"
        } ]
      }
    },
    "locations" : [ {
      "start" : 164,
      "end" : 174,
      "level" : "STRONG",
      "cigarAlignment" : "11M",
      "alignment" : "AYEPVWAIGTG",
      "location-fragments" : [ {
        "start" : 164,
        "end" : 174,
        "dc-status" : "CONTINUOUS"
      } ]
    } ],
    "model-ac" : "PS00171"
  }, {
    "signature" : {
      "accession" : "G3DSA:3.20.20.70",
      "name" : "Aldolase class I",
      "description" : null,
      "signatureLibraryRelease" : {
        "library" : "GENE3D",
        "version" : "4.2.0"
      },
      "entry" : {
        "accession" : "IPR013785",
        "name" : "Aldolase_TIM",
        "description" : "Aldolase-type TIM barrel",
        "type" : "HOMOLOGOUS_SUPERFAMILY",
        "goXRefs" : [ {
          "name" : "catalytic activity",
          "databaseName" : "GO",
          "category" : "MOLECULAR_FUNCTION",
          "id" : "GO:0003824"
        } ],
        "pathwayXRefs" : [ {
          "name" : "2-chloroacrylate degradation II",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7428"
        }, {
          "name" : "Geraniol degradation",
          "databaseName" : "KEGG",
          "id" : "00281"
        }, {
          "name" : "meleagrin biosynthesis",
          "databaseName" : "MetaCyc",
          "id" : "PWY-7607"
        } ]
      }
    },
    "locations" : [ {
      "start" : 2,
      "end" : 249,
      "hmmStart" : 3,
      "hmmEnd" : 247,
      "hmmLength" : 248,
      "hmmBounds" : "COMPLETE",
      "evalue" : 2.3E-108,
      "score" : 363.1,
      "envelopeStart" : 2,
      "envelopeEnd" : 249,
      "postProcessed" : true,
      "location-fragments" : [ {
        "start" : 2,
        "end" : 249,
        "dc-status" : "CONTINUOUS"
      } ]
    } ],
    "evalue" : 2.1E-108,
    "score" : 363.3,
    "model-ac" : "1r2rB00"
  }, {
    "signature" : {
      "accession" : "SSF51351",
      "name" : "Triosephosphate isomerase (TIM)",
      "description" : null,
      "signatureLibraryRelease" : {
        "library" : "SUPERFAMILY",
        "version" : "1.75"
      },
      "entry" : {
        "accession" : "IPR035990",
        "name" : "TIM_sf",
        "description" : "Triosephosphate isomerase superfamily",
        "type" : "HOMOLOGOUS_SUPERFAMILY",
        "goXRefs" : [ {
          "name" : "triose-phosphate isomerase activity",
          "databaseName" : "GO",
          "category" : "MOLECULAR_FUNCTION",
          "id" : "GO:0004807"
        } ],
        "pathwayXRefs" : [ {
          "name" : "1-butanol autotrophic biosynthesis (engineered)",
          "databaseName" : "MetaCyc",
          "id" : "PWY-6886"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-DME-70263"
        }, {
          "name" : "Biosynthesis of secondary metabolites",
          "databaseName" : "KEGG",
          "id" : "01110"
        } ]
      }
    },
    "locations" : [ {
      "start" : 3,
      "end" : 247,
      "hmmLength" : 257,
      "location-fragments" : [ {
        "start" : 3,
        "end" : 247,
        "dc-status" : "CONTINUOUS"
      } ]
    } ],
    "evalue" : 9.69E-96,
    "model-ac" : "0048198"
  }, {
    "signature" : {
      "accession" : "PF00121",
      "name" : "TIM",
      "description" : "Triosephosphate isomerase",
      "signatureLibraryRelease" : {
        "library" : "PFAM",
        "version" : "33.1"
      },
      "entry" : {
        "accession" : "IPR000652",
        "name" : "Triosephosphate_isomerase",
        "description" : "Triosephosphate isomerase",
        "type" : "FAMILY",
        "goXRefs" : [ {
          "name" : "triose-phosphate isomerase activity",
          "databaseName" : "GO",
          "category" : "MOLECULAR_FUNCTION",
          "id" : "GO:0004807"
        } ],
        "pathwayXRefs" : [ {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-DME-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-GGA-70171"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-MMU-70263"
        } ]
      }
    },
    "locations" : [ {
      "start" : 8,
      "end" : 245,
      "hmmStart" : 2,
      "hmmEnd" : 243,
      "hmmLength" : 243,
      "hmmBounds" : "C_TERMINAL_COMPLETE",
      "evalue" : 1.8E-86,
      "score" : 289.4,
      "envelopeStart" : 7,
      "envelopeEnd" : 245,
      "postProcessed" : true,
      "location-fragments" : [ {
        "start" : 8,
        "end" : 245,
        "dc-status" : "CONTINUOUS"
      } ]
    } ],
    "evalue" : 1.6E-86,
    "score" : 289.6,
    "model-ac" : "PF00121"
  }, {
    "signature" : {
      "accession" : "TIGR00419",
      "name" : "TIGR00419",
      "description" : "tim: triose-phosphate isomerase",
      "signatureLibraryRelease" : {
        "library" : "TIGRFAM",
        "version" : "15.0"
      },
      "entry" : {
        "accession" : "IPR000652",
        "name" : "Triosephosphate_isomerase",
        "description" : "Triosephosphate isomerase",
        "type" : "FAMILY",
        "goXRefs" : [ {
          "name" : "triose-phosphate isomerase activity",
          "databaseName" : "GO",
          "category" : "MOLECULAR_FUNCTION",
          "id" : "GO:0004807"
        } ],
        "pathwayXRefs" : [ {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-DME-70171"
        }, {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-GGA-70171"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-MMU-70263"
        } ]
      }
    },
    "locations" : [ {
      "start" : 8,
      "end" : 240,
      "hmmStart" : 1,
      "hmmEnd" : 228,
      "hmmLength" : 228,
      "hmmBounds" : "COMPLETE",
      "evalue" : 2.2E-77,
      "score" : 257.9,
      "envelopeStart" : 8,
      "envelopeEnd" : 240,
      "postProcessed" : false,
      "location-fragments" : [ {
        "start" : 8,
        "end" : 240,
        "dc-status" : "CONTINUOUS"
      } ]
    } ],
    "evalue" : 2.0E-77,
    "score" : 258.0,
    "model-ac" : "TIGR00419"
  }, {
    "signature" : {
      "accession" : "cd00311",
      "name" : "TIM",
      "description" : "TIM",
      "signatureLibraryRelease" : {
        "library" : "CDD",
        "version" : "3.18"
      },
      "entry" : {
        "accession" : "IPR000652",
        "name" : "Triosephosphate_isomerase",
        "description" : "Triosephosphate isomerase",
        "type" : "FAMILY",
        "goXRefs" : [ {
          "name" : "triose-phosphate isomerase activity",
          "databaseName" : "GO",
          "category" : "MOLECULAR_FUNCTION",
          "id" : "GO:0004807"
        } ],
        "pathwayXRefs" : [ {
          "name" : "Glycolysis",
          "databaseName" : "Reactome",
          "id" : "R-DME-70171"
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-SCE-70263"
        }, {
          "name" : "Gluconeogenesis",
          "databaseName" : "Reactome",
          "id" : "R-MMU-70263"
        } ]
      }
    },
    "locations" : [ {
      "start" : 7,
      "end" : 246,
      "sites" : [ {
        "description" : "catalytic triad",
        "numLocations" : 3,
        "siteLocations" : [ {
          "start" : 166,
          "end" : 166,
          "residue" : "E"
        }, {
          "start" : 96,
          "end" : 96,
          "residue" : "H"
        }, {
          "start" : 14,
          "end" : 14,
          "residue" : "K"
        } ]
      }, {
        "description" : "dimer interface",
        "numLocations" : 13,
        "siteLocations" : [ {
          "start" : 83,
          "end" : 83,
          "residue" : "M"
        }, {
          "start" : 50,
          "end" : 50,
          "residue" : "D"
        }, {
          "start" : 47,
          "end" : 47,
          "residue" : "A"
        }, {
          "start" : 53,
          "end" : 53,
          "residue" : "R"
        }, {
          "start" : 86,
          "end" : 86,
          "residue" : "D"
        }, {
          "start" : 46,
          "end" : 46,
          "residue" : "T"
        }, {
          "start" : 99,
          "end" : 99,
          "residue" : "R"
        }, {
          "start" : 65,
          "end" : 65,
          "residue" : "Q"
        }, {
          "start" : 48,
          "end" : 48,
          "residue" : "Y"
        }, {
          "start" : 87,
          "end" : 87,
          "residue" : "C"
        }, {
          "start" : 98,
          "end" : 98,
          "residue" : "E"
        }, {
          "start" : 15,
          "end" : 15,
          "residue" : "M"
        }, {
          "start" : 12,
          "end" : 12,
          "residue" : "N"
        } ]
      }, {
        "description" : "substrate binding site",
        "numLocations" : 9,
        "siteLocations" : [ {
          "start" : 172,
          "end" : 172,
          "residue" : "G"
        }, {
          "start" : 231,
          "end" : 231,
          "residue" : "L"
        }, {
          "start" : 233,
          "end" : 233,
          "residue" : "G"
        }, {
          "start" : 166,
          "end" : 166,
          "residue" : "E"
        }, {
          "start" : 234,
          "end" : 234,
          "residue" : "G"
        }, {
          "start" : 212,
          "end" : 212,
          "residue" : "S"
        }, {
          "start" : 96,
          "end" : 96,
          "residue" : "H"
        }, {
          "start" : 14,
          "end" : 14,
          "residue" : "K"
        }, {
          "start" : 12,
          "end" : 12,
          "residue" : "N"
        } ]
      } ],
      "evalue" : 1.8089E-130,
      "score" : 366.476,
      "location-fragments" : [ {
        "start" : 7,
        "end" : 246,
        "dc-status" : "CONTINUOUS"
      } ]
    } ],
    "model-ac" : "cd00311"
  }, {
    "signature" : {
      "accession" : "PTHR21139:SF24",
      "name" : "TRIOSEPHOSPHATE ISOMERASE",
      "description" : null,
      "signatureLibraryRelease" : {
        "library" : "PANTHER",
        "version" : "15.0"
      },
      "entry" : null
    },
    "locations" : [ {
      "start" : 1,
      "end" : 249,
      "hmmStart" : 1,
      "hmmEnd" : 249,
      "hmmLength" : 249,
      "hmmBounds" : "COMPLETE",
      "envelopeStart" : 1,
      "envelopeEnd" : 249,
      "location-fragments" : [ {
        "start" : 1,
        "end" : 249,
        "dc-status" : "CONTINUOUS"
      } ]
    } ],
    "evalue" : 0.0,
    "familyName" : "Not available",
    "score" : 585.3,
    "model-ac" : "PTHR21139:SF24"
  } ],
  "xref" : [ {
    "name" : "CAA49379.1 triosephosphate isomerase [Homo sapiens]",
    "id" : "CAA49379.1"
  } ]
} ]
}'''


# ==================================================================================================
# Testing
# ==================================================================================================
if __name__ == '__main__':
    testseq = []
    testseq.append('''>AAX90616.1 Src [Mus musculus]
MGSNKSKPKDASQRRRSLEPSENVHGAGGAFPASQTPSKPASADGHRGPSAAFVPPAAEPKLFGGFNSSD
TVTSPQRAGPLAGGVTTFVALYDYESRTETDLSFKKGERLQIVNNTRKVDVREGDWWLAHSLSTGQTGYI
PSNYVAPSDSIQAEEWYFGKITRRESERLLLNAENPRGTFLVRESETTKGAYCLSVSDFDNAKGLNVKHY
KIRKLDSGGFYITSRTQFNSLQQLVAYYSKHADGLCHRLTTVCPTSKPQTQGLAKDAWEIPRESLRLEVK
LGQGCFGEVWMGTWNGTTRVAIKTLKPGTMSPEAFLQEAQVMKKLRHEKLVQLYAVVSEEPIYIVTEYMN
KGSLLDFLKGETGKYLRLPQLVDMSAQIASGMAYVERMNYVHRDLRAANILVGENLVCKVADFGLARLIE
DNEYTARQGAKFPIKWTAPEAALYGRFTIKSDVWSFGILLTELTTKGRVPYPGMVNREVLDQVERGYRMP
CPPECPESLHDLMCQCWRKEPEERPTFEYLQAFLEDYFTSTEPQYQPGENL
''')
    testseq.append('''>QAB35645.1 BRI1 protein [Solanum tuberosum]
MKAHKTVFYQHPLSLNKLFFVLLLIFFLPPASPASVNGLFKDSQQLLSFKAALPPTPTLLQNWLPSTDPC
SFTGVSCKNSRVSSIDLSNTFLSVDFSLVTSYLLPLSNLESLVLKNANLSGSLTSAAKSQCGVSLDSIDL
AENTISGPISDISSFGVCSNLKSLNLSKNFLDPPGKEILKGATFSLQVLDLSYNNISGFNLFPWVSSMGF
GELEFFSLKGNKLAGSIPELDFKNLSHLDLSANNFSTVFPSFKDCSNLQHLDLSSNKFYGDIGSSLSSCG
KLSFLNLTNNQFVGLVPKLQSESLQYLYLRGNDFQGVYPNQLADLCKTVVELDLSYNNFSGMVPESLGEC
SSLELVDISNNNFSGKLPVDTLLKLSNMKTMVLSFNKFVGVLPDSFSNLLKLETLDVSSNNLTGVIPSGI
CKDPMNNLKVLYLQNNLFEGPIPDSLSNCSQLVSLDLSFNYLTRRIPSSLGSLSKLKDLILWLYQLSGEI
PQELMYLQALENLILDFNDLTGPIPASLSNCTKLNWISLSNNQLSGEIPASLGRLSNLAILKLGNNSISG
NIPAELGNCQSLIWLDLNTNFLSGSIPPPLFKQSGNIAVALLTGKRYVYIKNDGSKECHGAGNLLEFGGI
RQEQLGRISTRHPCNFTRVYRGITQPTFNHNGSMIFLDLSYNKLEGSIPKELGTMYYLSILNLGHNDLSG
MIPQDHGGLKNVAILDLSYNRFNGPIPNSLTSLTLLGEIDLSNNNLSGMIPESAPFDTFPDYRFANNSLC
GYPLPLPCSSGPKSDANQHQKSHRRQASLAGSVAMGLLFSLFCIFGLIIVAIETKKRRKKKEAALEAYMD
GHSHSATANSAWKFTSAREALSINLAAFEKPLRKLTFADLLEATNGFHNDSLVGSGGFGDVYKAQLKDGS
VVAIKKLIHVSGQGDREFTAEMETIGKIKHRNLVPLLGYCKVGEERLLVYEYMKYGSLEDVLHDRKKIGI
KLNWPARRKIAIGAARRLAFLHHNCIPHIIHRDMKSSNVLLDENLEARVSDFGMARLMSAMDTHLSVSTL
AGTPGYVPPEYYQSFRCSTKGDVYSYGVVLLELLTGKQPTDSADFGDNNLVGWVKLHAKGKITDVFDREL
LKEDPSIEIELLQHLKVACACLDDRHWKRPTMIQVMAMFKEIQAGSGMDSTSTIGADDVNFSAVEGGIEM
GINESIKEGNELSKHL''')

    #nucleotide queries are no longer supported at interpro
    nucleotide_test = '''>AB512659.1 Bos javanicus birmanicus x Bos indicus HBB gene for hemoglobin beta, complete cds, allele: HBB-B, specimen_voucher: personal:Cambodia Hybrid-CamII-2
ACAAACAGACACCATGCTGACTGCTGAGGAGAAGGCTGCCGTCACCGCCTTTTGGAGCAAGGTGCATGTG
GATGAAGTTGGTGGTGAGGCCCTGGGCAGGTAGGTATCCCACTTACAAGGCAGGTTTAAGGAGAGTGAAA
TGCACCTGGGCGTGTGAGGACAGAGCCGTCCCTGAGATTCTGAAAGCTGCTGGCTTCCTCTGACCTTGTG
CTGTTTTCTCCCCCTAGGCTGCTGGTTGTCTACCCCTGGACTCAGAGGTTCTTTGAGTCCTTTGGGGACT
TGTCCACTGCTGATGCTGTTATGAACAACCCTAAGGTGAAGGCCCATGGCAAGAAGGTGCTAGATTCCTT
TAGTAATGGCATGAAGCATCTCGATGACCTCAAGGGCACCTTTGCTGCGCTGAGTGAGCTGCACTGTGAT
AAGCTGCATGTGGATCCTGAGAACTTCAAGGTGAGTTTGTGGAATCCTCAGTGTTCTTCTTCTTTTTATG
GTCAAGCTCATGTCATGAGGAGAAAGCTGAATGGCAGGACACAGTTTAGAATGGAGAAGAGGTATTCTGG
TTAGATTACTAAGGACTCCTCAGAACCGTTTAGACTCTTTTAACCTCTTTGTTCACAACCAGCATTTCCT
CTGATTCATTCTTGTTCTCTGTTGTCTGCAATGTCCTCTTTTTAATTATATTTTTTATTTTGAGGGTTTA
ATTTGAAAAAAAATTATATATCAACTTTAAAAATTGTATCTAATATTTCCCCCTTATCTGTTTCTTTCAA
GGAATAAAATGTTCTATTGCTTTTTGAAATGATTCAAAATAATAAAAATAATAACAAGTTCTGGATTAAG
TTAGAAAGAGAGAAACATTTCTAAATATATATTCAGGAAGATATAGGTAGATTCACATCAGTAGTAACAA
CTTCACTTCAGTCATCTTTGTGCTTATATCCTACGGTCACAGCTTGGGATAAGACTGAAATACCCTGAAT
CTAACCTTGGATTTCCCTCATAGCTCAGTTGGTTAAGCATCTGCCTGCAATGCAAGAGATCCCAGTTCGA
TTCCTGGGTCGGGAAGGATGGCTGGAGAAGGGATAGGCACCCACTCCAGTATTCTTGGGTTTCCCTTGTG
GCTCAGCTGGTAAAGAATCTGCCTGCAATGTGGGAGACCCAGCTTCTATCCCTGAGTTTGGAGGATCCCC
TGGAGAAGGGAAAGGCTACCCACTCCAGTATTCTGGCCTGGAGAAATCTATGGACTGTATAGTCCATGGG
GTTGCAAAGAATCAGACACGATTGAGAGACTCTCACTTCACTCACCTGCACTAACCCTGCCCTTGCTTAA
TGTCTTTTCCACACAGCTCCTGGGCAACGTGCTAGTGGTTGTGCTGGCTCGCAATTTTGGCAACGAATTC
ACCCCGGTGCTGCAGGCTGACTTTCAGAAGGTGGTGGCTGGTGTGGCCAATGCCCTGGCCCACAGATATC
ATTAAGCTCCCTTTCCTGCTTTCCAGGAAAGGTTTTTTCATCCTCAGAGCCCAAAAATTGAATATGGAAA
AATTATGAAGTGTTTTGAGCATCTGGCCTCTGCCTAATAAAGACATTTATTTTCATTGCACTGGTGTATT
TAAATTATTTCACTGTCTCTTACTCAGATGGGCACATGGGAGGGCAAAACACTGAAGACATAAAGAAATG
AAGG
'''
    ips = Interpro()
    ips.email = 'gribskov@purdue.edu'
    ips.title = 'BRI1'
    ips.sequence = testseq[1]
    ips.application_select(['TIGRFAM', 'CDD', 'PfamA'])
    # ips.application_select(['Phobius', 'ProSitePatterns'])
    ips.output_select('json')
    ips.parameter_select({'goterms':True, 'pathways':True})

    print('submitting')
    if not ips.submit(show_query=True):
        exit(1)

    poll_time = 20
    poll_count = 0
    poll_max = 50
    while ips.status() != 'finished':
        time.sleep(poll_time)
        poll_count += 1
        print(f'\t ... polling({poll_count}) = {ips.jobstatus}')
        if poll_count > poll_max:
            break

    print('collecting result')
    ips.result()

    # parse and print the result, comment out the above and ncomment the next line to test parsing
    # the result without running a query
    # ips.content = json_test()
    parsed_result = ips.parse_json()

    for eachmotif in parsed_result['motifs']:
        print('{ipr_accession}\t{src_accession}\t{description}'.format(**eachmotif))
    for goterm in parsed_result['go']:
        go = parsed_result['go'][goterm]
        print('{}\t{}\t{}'.format(goterm, go['name'], go['source']))
    for path in parsed_result['pathway']:
        pathway = parsed_result['pathway'][path]
        print('{}\t{}\t{}'.format(path, pathway['name'], pathway['source']))

    print('done')

    exit(0)
