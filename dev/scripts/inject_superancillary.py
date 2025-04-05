from pathlib import Path 
import json 
FLUIDS = Path(__file__).parent.parent / 'fluids'
print(FLUIDS)
assert(FLUIDS.exists())
FASTCHEB_OUT = Path(r"/Users/ianbell/Documents/Code/fastchebpure/output")
assert(FASTCHEB_OUT.exists())
for path in FASTCHEB_OUT.glob('*.json'):
    print(path)
    name = path.stem.replace('_exps','')
    contents = json.load(path.open())
    destpath = FLUIDS / f'{name}.json'
    dest = json.load(destpath.open())
    dest['EOS'][0]['SUPERANCILLARY'] = contents
    with destpath.open('w') as fp:
        fp.write(json.dumps(dest, indent=2))
