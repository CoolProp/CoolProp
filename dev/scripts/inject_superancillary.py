from pathlib import Path 
import json 
FLUIDS = Path(__file__).parent / 'fluids'
assert(FLUIDS.exists())
FASTCHEB_OUT = Path(r"D:\Code\fastchebpure")
for path in FASTCHEB_OUT.glob('*.json'):
    name = path.stem.replace('_exps','')
    contents = json.load(path.open())
    destpath = FLUIDS / f'{name}.json'
    dest = json.load(destpath.open())
    dest['EOS'][0]['SUPERANCILLARY'] = contents
    with destpath.open('w') as fp:
        fp.write(json.dumps(dest, indent=2))