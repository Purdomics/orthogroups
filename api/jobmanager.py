import pickle as pkl
import sys
import time


class Jobmanager:
    """=============================================================================================
    jobmanager

    The main datastructure in the jobmanager is a dictionary of jobs.  Each element in the list
    is an object that conforms to the JobManager_API.  The keys are the objects and the values
    are the last determined status of the job

    Jobs communicate with the jobmanger using JobManager_API.message

    manage submission of jobs using defined APIs. APIs should implement
        submit  submit job to server
        status  poll job status
        result  retrieve result from server
        clone   used to duplicate from template, basic version inherited from JobManager_API
        poke    identify class, inherited from JobManager_API

    Michael Gribskov     02 April 2021
    ============================================================================================="""

    def __init__(self, api=None, loglevel=1, log_fh=sys.stderr, poll_delay=61, poll_pause=0,
                 poll_max=50):
        """-----------------------------------------------------------------------------------------
        api = class with JobManagerAPI type
        -----------------------------------------------------------------------------------------"""

        self.joblist = {}
        self.api = api
        self.template = None

        self.loglevel = loglevel
        self.log_fh = log_fh
        self.poll_delay = poll_delay  # seconds between polling cycles
        self.poll_pause = poll_pause  # seconds between polling queries
        self.poll_max = poll_max  # maximum number of times to poll
        self.poll_count = 0  # number of times this job has been polled

    def log_message(self, message):
        """-----------------------------------------------------------------------------------------
        Write a message to the log. Messages are only written if self.loglevel >= loglevel
        Messages are stored in the client message attribute as a dict with the following keys
        message = { type, text, loglevel}

        types:
            not_available (for submission options)
            submitted
            polling
            finished
            retrieved
            server_error

        :type: string, type of message
        :param text: string, text of message
        :param loglevel: int, minimum loglevel to write message to log
        :return: True if message was written
        -----------------------------------------------------------------------------------------"""
        if self.loglevel >= message['loglevel']:
            event_time = time.strftime('%d/%b/%G:%H:%M:%S', time.localtime(time.time()))
            self.log_fh.write('{:24s}\t{:<16s}\t{}\n'.format(event_time, message['type'],
                                                             message['text']))
            return True

        return False

    def new_job(self, *args, **kwds):
        """-----------------------------------------------------------------------------------------
        Return a new object from the current API.  pass positional and keyword args through to
        the object constructor

        :return: object from self.api
        -----------------------------------------------------------------------------------------"""
        return self.api(*args, **kwds)

    def new_job_from_template(self):
        """-----------------------------------------------------------------------------------------
        Return a copy of the template object.  The object must implement clone().
        If template has not been loaded, return the default object

        :return: object from self.template
        -----------------------------------------------------------------------------------------"""
        if self.template:
            return self.template.clone()
        else:
            return self.api()

    def poll_all(self, wait_all=True):
        """-----------------------------------------------------------------------------------------
        Poll all jobs in the jobs list and record their status. If wait_all is True, poll until all
        have finished.

        :param wait_all: boolean, if True poll until all jobs finish
        :return: int, number of finished jobs
        -----------------------------------------------------------------------------------------"""

        not_done = True
        n = len(self.joblist)
        while not_done:
            time.sleep(self.poll_delay)
            not_done = False
            n_finished = 0

            for job in self.joblist:
                self.joblist[job] = job.status()
                self.log_message(job.message)
                if job.status() == 'finished':
                    n_finished += 1
                time.sleep(self.poll_pause)

            if wait_all:
                not_done = n > n_finished

        return n_finished

    def save_all(self, reformat=None, fh=None, pickle=False, remove=True):
        """-----------------------------------------------------------------------------------------
        Process the output of all finished jobs.
        Reformat is a callback function used to reformat the output.  for instance,
        interpro.parse_json
        result is returned as a dictionary with the jobname as the ID
        If fh is True, output is written to the filehandle after reformatting.
        if pickle is true, the result is pickled before writing
        If remove is true, jobs are deleted from the list after saving

        :param joblist: dict, JobManager_API object is key, status is value
        :param reformat: function, callback function for formatting job result
        :param fh: filehandle for writable file
        :param pickle: boolean, True to pickle
        :param remove: boolean, remove finished jobs from joblist after saving
        :return: string, text of job content
        -----------------------------------------------------------------------------------------"""
        delete_list = []
        text = ''
        allresult = {}
        for job in self.joblist:
            if self.joblist[job] != 'finished':
                # skip unfinished jobs
                continue

            # self.joblist[job] = job.result()  # retrieve the completed job
            # text = job.content
            if reformat:
                text = reformat(job)
                if text:
                    allresult[job.jobname] = text
                else:
                    allresult[job.jobname] = ''

                if pickle:
                    if fh.mode == 'wb':
                        pkl.dump(text, fh)
                    else:
                        sys.stderr.write(
                            'Jobmanager:save_all - pickle file must be in opened in "wb" mode')
                elif fh:
                    if text:
                        fh.write('!{} - {}s\n'.format(job.jobname, job.jobid))
                        fh.write('{}\n'.format(text))
                    else:
                        fh.write('!{} - {} result is empty\n'.format(job.jobname, job.jobid))

            else:
                # raw result without reformatting
                allresult[job.jobname] = job.content

            if remove:
                # mark the job for deletion from joblist
                delete_list.append(job)

        for job in delete_list:
            # delete jobs in delete_list from joblist
            del self.joblist[job]

        return allresult

    def start(self, job):
        """-----------------------------------------------------------------------------------------
        Add a job to the joblist and submit via job.submit according to the JobManager_API.  The
        job object is the key and the value is the status

        :param job: object from api
        :return: int, number of jobs in the joblist
        -----------------------------------------------------------------------------------------"""
        self.joblist[job] = job.submit()
        self.log_message(job.message)

        return len(self.joblist)


# --------------------------------------------------------------------------------------------------
# testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    import sys
    from api.blast.blastncbi import BlastNCBI
    from api.interpro.interpro import Interpro

    testinterpro = False
    testblast = True

    src = '''>AAX90616.1 Src [Mus musculus]
    MGSNKSKPKDASQRRRSLEPSENVHGAGGAFPASQTPSKPASADGHRGPSAAFVPPAAEPKLFGGFNSSD
    TVTSPQRAGPLAGGVTTFVALYDYESRTETDLSFKKGERLQIVNNTRKVDVREGDWWLAHSLSTGQTGYI
    PSNYVAPSDSIQAEEWYFGKITRRESERLLLNAENPRGTFLVRESETTKGAYCLSVSDFDNAKGLNVKHY
    KIRKLDSGGFYITSRTQFNSLQQLVAYYSKHADGLCHRLTTVCPTSKPQTQGLAKDAWEIPRESLRLEVK
    LGQGCFGEVWMGTWNGTTRVAIKTLKPGTMSPEAFLQEAQVMKKLRHEKLVQLYAVVSEEPIYIVTEYMN
    KGSLLDFLKGETGKYLRLPQLVDMSAQIASGMAYVERMNYVHRDLRAANILVGENLVCKVADFGLARLIE
    DNEYTARQGAKFPIKWTAPEAALYGRFTIKSDVWSFGILLTELTTKGRVPYPGMVNREVLDQVERGYRMP
    CPPECPESLHDLMCQCWRKEPEERPTFEYLQAFLEDYFTSTEPQYQPGENL'''

    if testinterpro:
        joblist = Jobmanager(Interpro)
        joblist.poll_delay = 15
        joblist.poll_pause = 5
        joblist.loglevel = 2

        ips = joblist.new_job()
        print(f'interface is {ips.poke()}')
        ips.email = 'gribskov@purdue.edu'
        ips.title = 'mouse src'
        ips.sequence = src
        ips.application_select(['TIGRFAM', 'CDD', 'PfamA'])
        # ips.application_select(['Phobius', 'ProSitePatterns'])
        ips.output_select('json')
        ips.parameter_select({'goterms':True, 'pathways':True})

        # uncomment to run a new job
        # joblist.start(ips)
        # joblist.poll_all()
        # or use a known jobid to test
        ips.jobid = 'iprscan5-R20210512-150227-0752-35433138-p2m'
        ips.jobname = 'test'
        joblist.joblist[ips] = 'finished'

        ips.result()
        all = joblist.save_all(reformat=lambda x:x.response.text, remove=False, fh=sys.stdout)
        picklejar = open('jobmanager.test.pkl', 'wb')
        joblist.save_all(reformat=Interpro.parse_json, pickle=True, fh=picklejar)
        picklejar.close()
        picklejar = open('jobmanager.test.pkl', 'rb')
        unpickle = pkl.load(picklejar)
        picklejar.close()

    if testblast:
        jobs = Jobmanager(BlastNCBI)
        jobs.poll_delay = 60
        jobs.poll_pause = 5
        jobs.loglevel = 2

        blast = jobs.new_job()
        blast.email = 'gribskov@purdue.edu'
        blast.jobname= 'mouse src'
        blast.program = 'blastp'
        blast.database = 'pdb'
        blast.query = src

        # comment out to test with a known rid
        # jobs.start(blast)
        # jobs.poll_all()

        blast.jobid = '9U4Z1YDR013'
        jobs.joblist[blast] = 'finished'

        blast.result()
        all = jobs.save_all(reformat=lambda x:x.response.text, remove=False, fh=sys.stdout)
        all = jobs.save_all(BlastNCBI.parse_xml, remove=False, fh=sys.stdout)

    exit(0)
