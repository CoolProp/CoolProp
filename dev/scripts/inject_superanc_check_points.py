"""
Stamp SUPERANCILLARY freshness metadata into each dev/fluids/*.json from a
fastchebpure release: extended-precision `check_points` and an `eos_hash`.

Source of truth is fastchebpure's output/{fluid}_exps.json (for the eos_hash,
verbatim from `source_eos_hash`) and outputcheck/{fluid}_check.json (for a
dense T grid with (T, p_mp, rho'_mp, rho''_mp) columns computed at extended
precision, plus the (SA)/(mp) ratio columns that report how closely
fastchebpure's double-precision Chebyshev evaluation reproduces the
multi-precision reference at that T).

What gets injected per fluid:

    check_points: a small sample (default 3 points) from the fastchebpure
        check file. Picked in Theta = (T_crit_num - T) / T_crit_num space at
        values {0.5, 0.3, 0.1} — the expansion's sweet spot, clear of both
        the sub-microPascal heavy-oil regime near T_min and the last-piece
        endpoint near T_crit where production code switches to the
        critical-anchor spline. Each entry carries the fastchebpure columns
        "T / K", "p(mp) / Pa", "p(SA)/p(mp)", "rho'(mp) / mol/m^3",
        "rho'(SA)/rho'(mp)", "rho''(mp) / mol/m^3", "rho''(SA)/rho''(mp)".

    eos_hash: copied verbatim from fastchebpure's `source_eos_hash` field.
        fastchebpure computes this at fit time — FNV-1a structural hash of
        EOS[0] minus SUPERANCILLARY, using the byte-stream contract
        documented in eos_fnv1a_hex below — so it is the definitive stamp
        of the EOS fastchebpure actually fit the SA against. The C++ hash
        test and check_superanc_freshness.py recompute the same hash from
        the current master EOS and flag any later edit as SA drift.

Safety guard: before stamping, we compare the fluid JSON's current
jexpansions_p[0].coef[0] bit-exactly against the same entry in this
fastchebpure release's output/{fluid}_exps.json. Fluids whose SA has been
regenerated out-of-band (or simply not yet updated to the current release)
are skipped rather than stamped with mismatched metadata.

Belt-and-suspenders: when the SA guard passes we also recompute
eos_fnv1a_hex(current_master_eos) and assert it equals the emitted
source_eos_hash. Any disagreement there means either the algorithms on the
two sides have drifted, or master's EOS was edited after fastchebpure ran.
Either way it's worth a loud warning.

The fluid JSON is edited in place via text-splice (no json.dumps round-trip
for the surrounding file content), so numbers with preserved trailing zeros
or other stylistic choices elsewhere in the file are not perturbed.
Injection is idempotent — re-running overwrites any previous check_points
and eos_hash.

Usage:
    python inject_superanc_check_points.py --outputcheck /path/to/fastchebpure/outputcheck
"""
import argparse
import json
import re
import struct
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
FLUIDS = REPO / 'dev' / 'fluids'

PASSTHROUGH_KEYS = [
    'T / K',
    "p(mp) / Pa",
    "rho'(mp) / mol/m^3", "rho'(SA)/rho'(mp)",
    "rho''(mp) / mol/m^3", "rho''(SA)/rho''(mp)",
]

# Theta = (T_crit_num - T) / T_crit_num values at which to sample, spanning the
# SA expansion's usable range away from both the near-triple low-pressure regime
# and the near-critical spline takeover zone.
THETAS_DEFAULT = (0.5, 0.3, 0.1)


def pick_rows(rows, Tcrit_num, thetas):
    """For each Theta, return the index of the check-file row whose T is closest
    to T_crit_num * (1 - Theta). Rows are picked in T-space (not index-space)
    because the fastchebpure grid is log-spaced in Theta."""
    N = len(rows)
    if N < 2:
        return list(range(N))
    Ts = [r['T / K'] for r in rows]

    def nearest(target):
        return min(range(N), key=lambda i: abs(Ts[i] - target))

    return [nearest(Tcrit_num * (1.0 - theta)) for theta in thetas]


def extract_points(check_path, thetas):
    doc = json.loads(check_path.read_text())
    Tcrit_num = doc['meta']['Tcrittrue / K']
    data = doc['data']
    points = []
    for i in pick_rows(data, Tcrit_num, thetas):
        row = data[i]
        pt = {k: row[k] for k in PASSTHROUGH_KEYS}
        # The fastchebpure check file carries explicit (SA)/(mp) ratio columns for
        # rho' and rho'' but not for p, so we derive p(SA)/p(mp) from the two p columns.
        pt['p(SA)/p(mp)'] = row['p(SA) / Pa'] / row['p(mp) / Pa']
        points.append(pt)
    return points


def eos_fnv1a_hex(eos_obj):
    """FNV-1a 64-bit structural hash of an EOS JSON subtree, used to stamp
    "this is the EOS the superancillary was fit against". The SUPERANCILLARY
    subtree is excluded (it's the thing we're gating on).

    This implementation exists solely for cross-checking fastchebpure's
    `source_eos_hash` against the current master EOS; fastchebpure's emitted
    value is the source of truth and is stamped verbatim.

    --------------------------------------------------------------------------
    Why a structural hash instead of hashing a canonical JSON dump?

    We want the C++ test (CoolProp-Tests.cpp, "Superancillary eos_hash matches
    current EOS at bit level") to produce the same hash for the same EOS, so
    edits to a fluid file are caught cross-language. Hashing a JSON *string*
    fails that goal: Python's `repr(float)` and nlohmann::json::dump() use
    shortest-round-trip algorithms that occasionally disagree in the last
    digit (e.g., nlohmann emits "19673.920781104862" where Python emits
    "19673.92078110486" for the same double). Hashing the parsed VALUES
    (IEEE bits, not text) sidesteps that entirely — two implementations that
    parse the same JSON into the same doubles produce identical hashes.

    --------------------------------------------------------------------------
    Byte-stream contract (must stay in lockstep with the C++ TreeHasher and
    with fastchebpure's eos_hash.hpp):

        null      -> 'n'
        false     -> 'f'
        true      -> 't'
        integer   -> 'i' then int64 two's-complement bits as LE u64
        float     -> 'd' then IEEE-754 bits as LE u64
        string    -> 's' then LE u64 UTF-8 byte count then UTF-8 bytes
        array     -> 'a' then LE u64 length then each element walked
        object    -> 'o' then LE u64 size then, for each key in sorted order:
                     LE u64 UTF-8 byte count, UTF-8 bytes, walked value

    --------------------------------------------------------------------------
    Notes for future maintainers:

      * Type tags let us distinguish 0, 0.0, false, and "" — without them,
        multiple types could fold into the same byte sequence.
      * `isinstance(obj, bool)` MUST be checked before `isinstance(obj, int)`
        because Python's bool is a subclass of int; otherwise True would hash
        as integer 1.
      * Integers use `v & 0xFFFFFFFFFFFFFFFF`, mirroring C++'s cast of int64_t
        to uint64_t. For values in the int64 range (which covers every
        integer in CoolProp's fluid JSONs) the bit patterns match exactly.
      * Object keys are walked in `sorted(...)` order — lexicographic over
        UTF-8 byte sequences — matching C++'s std::map iteration order.
      * Length prefixes for arrays and strings prevent ambiguous parses:
        e.g. ["ab","c"] vs ["a","bc"] have distinct byte streams only because
        each entry is preceded by its length.
      * FNV-1a is NOT cryptographic; we only need determinism and reasonable
        distribution for change detection. Seed/prime are the standard 64-bit
        FNV-1a constants.
    """
    MASK = 0xffffffffffffffff
    PRIME = 0x100000001b3

    def mix_byte(h, b):
        h ^= b
        return (h * PRIME) & MASK

    def mix_bytes(h, data):
        for b in data:
            h = mix_byte(h, b)
        return h

    def mix_u64(h, v):
        v &= MASK
        for i in range(8):
            h = mix_byte(h, (v >> (i * 8)) & 0xff)
        return h

    def walk(h, obj):
        if obj is None:
            return mix_byte(h, ord('n'))
        if isinstance(obj, bool):  # must precede int: bool is a subclass of int
            return mix_byte(h, ord('t') if obj else ord('f'))
        if isinstance(obj, int):
            h = mix_byte(h, ord('i'))
            return mix_u64(h, obj)
        if isinstance(obj, float):
            # Feed raw IEEE-754 bits, not the textual representation. The LE u64
            # we emit is the same 8 bytes C++'s std::memcpy(&bits, &double, 8)
            # produces on a little-endian host.
            h = mix_byte(h, ord('d'))
            return mix_u64(h, struct.unpack('<Q', struct.pack('<d', obj))[0])
        if isinstance(obj, str):
            enc = obj.encode('utf-8')
            h = mix_byte(h, ord('s'))
            h = mix_u64(h, len(enc))
            return mix_bytes(h, enc)
        if isinstance(obj, list):
            h = mix_byte(h, ord('a'))
            h = mix_u64(h, len(obj))
            for v in obj:
                h = walk(h, v)
            return h
        if isinstance(obj, dict):
            h = mix_byte(h, ord('o'))
            h = mix_u64(h, len(obj))
            for k in sorted(obj):
                enc = k.encode('utf-8')
                h = mix_u64(h, len(enc))
                h = mix_bytes(h, enc)
                h = walk(h, obj[k])
            return h
        raise TypeError(f'unhashable JSON type: {type(obj).__name__}')

    stripped = {k: v for k, v in eos_obj.items() if k != 'SUPERANCILLARY'}
    return f'{walk(0xcbf29ce484222325, stripped):016x}'


def _find_superanc_close(text):
    """Return (start, end) of the SUPERANCILLARY object's enclosing braces, or
    None if absent. Uses a brace depth counter to find the matching close of the
    object that immediately follows the `"SUPERANCILLARY":` key."""
    m = re.search(r'"SUPERANCILLARY"\s*:\s*\{', text)
    if m is None:
        return None
    start = m.end() - 1  # index of the opening '{'
    depth = 0
    in_string = False
    escape = False
    for i in range(start, len(text)):
        c = text[i]
        if escape:
            escape = False
            continue
        if c == '\\' and in_string:
            escape = True
            continue
        if c == '"':
            in_string = not in_string
            continue
        if in_string:
            continue
        if c == '{':
            depth += 1
        elif c == '}':
            depth -= 1
            if depth == 0:
                return start, i
    return None


def inject_fluid(fluid_path, points, eos_hash):
    """Splice a check_points array (and an eos_hash stamp) into the SUPERANCILLARY
    block as a text edit, so the rest of the file (formatting, number
    representations, comments, etc.) is preserved byte-for-byte."""
    text = fluid_path.read_text()
    span = _find_superanc_close(text)
    if span is None:
        return False
    _, close_idx = span
    line_start = text.rfind('\n', 0, close_idx) + 1
    # Indent of the line containing SUPERANCILLARY's closing '}'
    close_indent = text[line_start:close_idx]
    entry_indent = close_indent + '  '
    item_indent = entry_indent + '  '

    # Drop any previous check_points / eos_hash so re-runs are idempotent.
    prefix = text[:close_idx]
    prefix = re.sub(r',\s*"check_points"\s*:\s*\[.*?\]\s*(?=[,}\s]*\Z)', '', prefix, count=1, flags=re.DOTALL)
    prefix = re.sub(r',\s*"eos_hash"\s*:\s*"[0-9a-f]+"\s*(?=[,}\s]*\Z)', '', prefix, count=1)
    # The existing last entry in SUPERANCILLARY has no trailing comma (it is the last
    # child of the object); snap off the whitespace that preceded the close brace so
    # we can place a comma immediately after that entry's own close brace.
    prefix = prefix.rstrip(' \t\n\r')

    lines = [',']
    if eos_hash is not None:
        lines.append(f'{entry_indent}"eos_hash": {json.dumps(eos_hash)},')
    lines.append(f'{entry_indent}"check_points": [')
    for i, pt in enumerate(points):
        lines.append(f'{item_indent}{{')
        kv = [f'{item_indent}  {json.dumps(k)}: {json.dumps(v)}' for k, v in pt.items()]
        lines.append(',\n'.join(kv))
        lines.append(f'{item_indent}}}{"," if i < len(points) - 1 else ""}')
    lines.append(f'{entry_indent}]')
    insertion = '\n'.join(lines) + '\n' + close_indent

    new_text = prefix + insertion + text[close_idx:]
    fluid_path.write_text(new_text)
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--outputcheck', required=True, type=Path,
                    help='Path to fastchebpure outputcheck/ directory')
    ap.add_argument('--output', type=Path,
                    help='Path to fastchebpure output/ directory (contains *_exps.json '
                         'with source_eos_hash). Defaults to the sibling "output" of '
                         '--outputcheck.')
    ap.add_argument('--thetas', type=float, nargs='+', default=list(THETAS_DEFAULT),
                    help=f'Theta = (Tc - T)/Tc values at which to sample (default: {list(THETAS_DEFAULT)})')
    args = ap.parse_args()

    if not args.outputcheck.is_dir():
        raise SystemExit(f'outputcheck directory not found: {args.outputcheck}')
    output_dir = args.output or args.outputcheck.parent / 'output'
    if not output_dir.is_dir():
        raise SystemExit(f'output directory not found: {output_dir}')

    injected = skipped_no_super = skipped_no_check = skipped_sa_mismatch = cross_check_warnings = 0
    for fluid_path in sorted(FLUIDS.glob('*.json')):
        check_path = args.outputcheck / f'{fluid_path.stem}_check.json'
        exps_path = output_dir / f'{fluid_path.stem}_exps.json'
        if not check_path.exists() or not exps_path.exists():
            skipped_no_check += 1
            continue
        fluid_doc = json.loads(fluid_path.read_text())
        if 'EOS' not in fluid_doc or not fluid_doc['EOS']:
            skipped_no_super += 1
            continue
        cur_sa = fluid_doc['EOS'][0].get('SUPERANCILLARY')
        if cur_sa is None:
            skipped_no_super += 1
            continue
        fc_sa = json.loads(exps_path.read_text())
        # Bit-exact guard: refuse to stamp metadata on a fluid whose SA has
        # been regenerated out-of-band and no longer matches this fastchebpure
        # release's output.
        if (not cur_sa.get('jexpansions_p')
                or not fc_sa.get('jexpansions_p')
                or cur_sa['jexpansions_p'][0]['coef'][0] != fc_sa['jexpansions_p'][0]['coef'][0]):
            skipped_sa_mismatch += 1
            print(f'  skipped {fluid_path.stem}: SA in fluid JSON differs from fastchebpure output '
                  f'(likely regenerated out-of-band; wait for the next fastchebpure release)')
            continue
        # fastchebpure's source_eos_hash is the authoritative stamp of the EOS
        # it fit against. Copy it verbatim; defense-in-depth: cross-check
        # against the hash of master's current EOS. Divergence means either
        # the two sides' hash algorithms have drifted (contract bug) or the
        # fluid EOS was edited after fastchebpure ran.
        source_eos_hash = fc_sa.get('source_eos_hash')
        if source_eos_hash is None:
            skipped_no_check += 1
            print(f'  skipped {fluid_path.stem}: fastchebpure output lacks source_eos_hash '
                  f'(need release with https://github.com/CoolProp/CoolProp/issues/2777 applied)')
            continue
        local_hash = eos_fnv1a_hex(fluid_doc['EOS'][0])
        if local_hash != source_eos_hash:
            cross_check_warnings += 1
            print(f'  WARN  {fluid_path.stem}: source_eos_hash={source_eos_hash} but local '
                  f'eos_fnv1a_hex={local_hash}; stamping fastchebpure value (freshness check '
                  f'will flag the drift)')
        points = extract_points(check_path, args.thetas)
        if inject_fluid(fluid_path, points, source_eos_hash):
            injected += 1
            print(f'  injected {len(points)} points (hash={source_eos_hash}): {fluid_path.stem}')
        else:
            skipped_no_super += 1

    print()
    print(f'injected: {injected}')
    print(f'skipped (no SUPERANCILLARY block): {skipped_no_super}')
    print(f'skipped (no fastchebpure data for this fluid): {skipped_no_check}')
    print(f'skipped (SA diverged from fastchebpure output): {skipped_sa_mismatch}')
    if cross_check_warnings:
        print(f'WARNINGS (source_eos_hash != local hash): {cross_check_warnings}')


if __name__ == '__main__':
    main()
