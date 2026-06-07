"""Enumerate the CoolProp Python import tree (module presence + public names)
for the import-tree parity acceptance gate (bd CoolProp-r9sq.3).

  python import_tree.py            # dump current env's CoolProp tree as JSON
Compare nanobind vs legacy; import_tree_legacy.json pins stock 7.2.0 as the
'before' reference.  Deviations must be closed or explicitly approved by Ian.
"""
import importlib, pkgutil, json
ROOT = "CoolProp"
def names(mod):
    return sorted(n for n in dir(mod) if not n.startswith('_'))
tree = {}
root = importlib.import_module(ROOT)
tree[ROOT] = names(root)
# submodules one level down
for m in pkgutil.iter_modules(root.__path__):
    full = ROOT + "." + m.name
    try:
        sm = importlib.import_module(full)
        tree[full] = names(sm)
    except Exception as e:
        tree[full] = ["<import error: %s>" % str(e)[:50]]
print(json.dumps(tree, indent=0))
