# Security policy

## Reporting a vulnerability

Please report security vulnerabilities privately by email to **ian.h.bell@gmail.com**, not in a public
issue, pull request or discussion.

Include what you can of:

- the CoolProp version or commit, and the interface you used (Python, the C API shared library, C++,
  another wrapper);
- a reproducer: the `PropsSI` / `AbstractState` call, fluid string, JSON, or file that triggers it;
- the sanitizer trace or crash output;
- a proposed fix, if you have one.

We aim to acknowledge reports within a week. Once a fix is ready it is merged to `master` and the
report is credited in the changelog, unless you would rather not be named.

## Supported versions

Releases are tagged from `master`, and fixes ship in the next release. Older releases are not patched
(the old `release` and `6.4.1` branches are not maintained); upgrade to the latest release.

## What counts as a security issue

In scope: memory-safety bugs (out-of-bounds reads or writes, use-after-free, double free), crashes,
hangs or unbounded memory use that can be triggered through CoolProp's inputs. Those inputs are fluid
strings, parameter names and state values passed to `PropsSI`, `HAPropsSI` or `AbstractState`;
arguments to the C API in `CoolPropLib.h`; fluid or mixture JSON (or REFPROP HMX.BNC departure-function
text) loaded at runtime; and the table cache files CoolProp reads back from disk.

Not security issues: inaccurate property values, NaN results, solver non-convergence that ends in an
error, or a thrown CoolProp exception (such as `CoolProp::ValueError`) for invalid input. That is how
CoolProp reports bad input. Please open a normal [issue](https://github.com/CoolProp/CoolProp/issues)
for these. A hang counts as a security issue only when the call never returns.

Configuration set by the embedding application (`set_config_*`, `set_config_as_json_string`,
environment variables, the REFPROP library path) is trusted. Behaviour that requires an attacker to
control it is not a vulnerability in CoolProp.
