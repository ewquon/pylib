#!/usr/bin/env python
import sys
import A2e

if len(sys.argv) <= 1:
    print('Usage:\n ',sys.argv[0],
          'dataset_name start_datetime [end_datetime] [filetype]',
          '\nwhere datetimes are specified as YYYYMMDD[HHMMSS].',
          'Example:\n ',sys.argv[0],
          'mmc/radar.z01.00 20131108 20131109 nc')
    sys.exit()

dataset = sys.argv[1]
starttime = sys.argv[2].ljust(14,'0')
if len(sys.argv) > 3:
    endtime = sys.argv[3].ljust(14,'0')
else:
    endtime = starttime
if len(sys.argv) > 4:
    filetype = sys.argv[4]
else:
    filetype = None

filter_arg = {
    'Dataset': dataset,
    'date_time': {
        'between': [starttime, endtime]
    },
}
if filetype is not None:
    filter_arg['file_type'] = filetype

a2e = A2e.A2e()
a2e.setup_cert_auth()
files = a2e.search(filter_arg)
print(files)

a2e.download_files(files, path='.')

