import json

jj = json.load(open('mixture_binary_pairs.json', 'r'))
CASpairs = [tuple(sorted([p['CAS1'], p['CAS2']])) for p in jj]

dupes = {}
for unique in set(CASpairs):
    if CASpairs.count(unique) > 1:
        dupes[unique] = CASpairs.count(unique) - 1  # number to remove

dupe_pairs = []
for i in range(len(jj) - 1, -1, -1):
    pair = tuple(sorted([jj[i]['CAS1'], jj[i]['CAS2']]))
    if pair in dupes:
        dupe_pairs.append(jj.pop(i))
        dupes[pair] -= 1
        if dupes[pair] == 0:
            dupes.pop(pair)

with open('old_BIP.json', 'w') as fp:
    fp.write(json.dumps(dupe_pairs, indent=2, sort_keys=True))

with open('mixture_binary_pairs.json', 'w') as fp:
    fp.write(json.dumps(jj, indent=2, sort_keys=True))

jj = json.load(open('mixture_binary_pairs.json', 'r'))
CASpairs = [(p['CAS1'], p['CAS2']) for p in jj]

dupes = {}
for unique in set(CASpairs):
    if CASpairs.count(unique) > 1:
        dupes[unique] = CASpairs.count(unique) - 1  # number to remove

print(dupes)
