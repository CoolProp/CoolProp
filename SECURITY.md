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

You should get an acknowledgement within a week. Once a fix is ready it is merged to `master` and the
report is credited in the changelog, unless you would rather not be named.

## Supported versions

Fixes go into `master` and ship in the next release. There are no maintenance branches, so older
releases are not patched; upgrade to the latest release.

## What counts as a security issue

In scope: memory-safety bugs (out-of-bounds reads or writes, use-after-free, double free), crashes,
hangs or unbounded memory use that can be triggered through CoolProp's inputs. Those inputs are fluid
strings, parameter names and state values passed to `PropsSI`, `HAPropsSI` or `AbstractState`;
arguments to the C API in `CoolPropLib.h`; fluid or mixture JSON loaded at runtime; and the table cache
files CoolProp reads back from disk.

Not security issues: inaccurate property values, NaN results, solver non-convergence, or a thrown
`CoolProp::ValueError` for invalid input. That is how CoolProp reports bad input. Please open a normal
[issue](https://github.com/CoolProp/CoolProp/issues) for these.

Configuration set by the embedding application (`set_config_string`, environment variables, the REFPROP
library path) is trusted. Behaviour that requires an attacker to control it is not a vulnerability in
CoolProp.
