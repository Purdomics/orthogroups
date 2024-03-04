"""=================================================================================================
run interproscan using rest API
see https://www.ebi.ac.uk/seqdb/confluence/display/WEBSERVICES/iprscan5_rest\
================================================================================================="""
import requests

run = 'https://www.ebi.ac.uk/Tools/services/rest/iprscan5/run/'
email = 'gribskov@purdue.edu'
title = 'globin'
sequence = '''>sp|P69905.2|HBA_HUMAN RecName: Full=Hemoglobin subunit alpha; AltName: Full=Alpha-globin; AltName: Full=Hemoglobin alpha chain
MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNA
VAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSK
YR'''

# send the initial query
command = {'email': email, 'title': title, 'sequence': sequence}
response = requests.post(run, command)
id = response.text
print('job {} submitted'.format(id))

# poll for job completion
import time
maxtries = 10
notready = 1

status = 'http://www.ebi.ac.uk/Tools/services/rest/iprscan5/status/'

while notready:
    response = requests.get(status + id)
    print('    polling... response->{}'.format(response.text))

    if 'FINISHED' in response.text:
        notready = 0
        break
    else:
        notready += 1
    if notready >= maxtries:
        break

    # don't poll too often
    time.sleep(20)

if notready > 0:
    # polling reached limit
    print('unable to find result () in {} tries'.format(rid, notready))
    exit(1)
else:
    print('interproscan {} finished'.format(id))

# get the final result
result = 'http://www.ebi.ac.uk/Tools/services/rest/iprscan5/result/'
result += id + '/{}'.format('xml')
response = requests.get(result)
print('\n', response.text)
