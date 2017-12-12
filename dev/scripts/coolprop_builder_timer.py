from __future__ import print_function

import urllib, json
filehandle = urllib.urlopen('http://www.coolprop.dreamhosters.com:8010/json/builders')
jj = json.loads(filehandle.read())

times = []
for key in jj.keys():
    filehandle = urllib.urlopen('http://www.coolprop.dreamhosters.com:8010/json/builders/' + key + '/builds/-1')
    builder = json.loads(filehandle.read())
    elapsed_time = builder['times'][1] - builder['times'][0]
    times.append((elapsed_time, key))

print(sorted(times)[::-1])
