import CoolProp
from chemspipy import ChemSpider
from chemspipy_key import key  # private file with the key (DO NOT COMMIT!!)
import glob, json
cs = ChemSpider(key)

# Map from name to Chemspider ID
backup_map = {
    'Propyne': 6095,
    'R236EA': 71342,
    'R245ca': 62827,
    'trans-2-Butene': 56442,
    'Oxygen': 952,
    'Fluorine': 22932,
    'Hydrogen': 762,
    'Deuterium': 22931,
    'HFE143m': 66577,
    'SulfurHexafluoride': 16425,
    'R114': 13853215
}

# Make sure the key works
c = cs.get_compound(2157)
assert(c.inchikey == 'BSYNRYMUTXBXSQ-UHFFFAOYAW')

for fname in glob.glob('../fluids/*.json'):
    with open(fname, 'r') as fp:
        jj = json.load(fp)

    fluid = jj['INFO']['NAME']

    def doset(result):
        # A CAS search can return an unrelated compound (this is how several
        # fluids once got wrong InChIKeys and ChemSpider IDs; see COO-171).
        # The stored InChIKey has been checked against PubChem and RDKit, so
        # it is the reference: never overwrite a record whose key disagrees.
        # Compare the *standard* InChIKey: older chemspipy returns a
        # non-standard key in .inchikey (e.g. ...-UHFFFAOYAW), which would
        # never match and so silently reject every update.
        # A record with no stored key has nothing to check against and is
        # written unchecked.
        stored_key = jj['INFO'].get('INCHI_KEY')
        hit_key = getattr(result, 'stdinchikey', None) or result.inchikey
        if stored_key and hit_key != stored_key:
            print('%s: ChemSpider hit %s has InChIKey %s, stored is %s; not updating'
                  % (fluid, result.csid, hit_key, stored_key))
            return False
        jj['INFO']['INCHI_STRING'] = result.inchi
        jj['INFO']['INCHI_KEY'] = hit_key
        jj['INFO']['CHEMSPIDER_ID'] = result.csid
        jj['INFO']['2DPNG_URL'] = result.image_url
        jj['INFO']['SMILES'] = result.smiles
        return True

    CAS = CoolProp.CoolProp.get_fluid_param_string(fluid, "CAS")
    if '.' not in CAS:
        results = cs.search(CAS)
        results.wait()
        if fluid in backup_map:
            results = cs.search(backup_map[fluid])
            results.wait()
            assert(len(results) == 1)
            doset(results[0])
        if len(results) == 1:
            doset(results[0])
        elif fluid in backup_map:
            results = cs.search(backup_map[fluid])
            results.wait()
            assert(len(results) == 1)
            doset(results[0])
        else:
            print('%s %s !!failure!! %s' % (fluid, CAS, len(results)))
            for result in results:
                spectra = cs.get_compound_spectra(result.csid)
                if spectra and '##CAS REGISTRY NO=' + CAS in spectra[0].data:
                    if doset(result):
                        print('GOT IT!!')
                        break
                print(result.common_name, result.inchikey, result.stdinchi, cs.get_extended_compound_info(result.csid))
            print('')

    with open(fname, 'w') as fp:
        json.dump(jj, fp, indent=2, sort_keys=True)

    del jj, fp
