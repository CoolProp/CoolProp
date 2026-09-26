#!/usr/bin/env python3
"""Decide whether an OBS build result list is a pass or a failure.

This exists because "osc results --fail-on-error" cannot be used for the job
CI needs done.  That flag is evaluated inside osc's polling loop:

    for results in get_package_results(...):
        ...
            if res['code'] in ('failed', 'broken', 'unresolvable'):
                failed = True

(osc/core.py, get_results).  The flag is never cleared, so a state left over
from BEFORE the run uploaded anything latches a failure on the first poll and
no later success can undo it.  A package sitting at "broken: no source
uploaded" therefore reports red even when every target ends up succeeded,
which is exactly what happened on the first real run of the CI loop.

So the verdict is taken once, from the final result list, rather than from
whatever was observed on the way there.

Exit codes:

    0   at least one target could build and none of them failed
    1   a target failed, or there was nothing to build at all
    2   the results have not settled yet, so the caller should re-read them

Reading the arguments: a repository that is "excluded" or "disabled" is a
target OBS was told not to build, which is a configuration choice rather than
a regression, so it is neither a failure nor evidence that anything was built.
A list where EVERY target is like that is a failure: a package with no enabled
repository, or a mistyped package name, would otherwise be a green tick for a
build that happened nowhere.
"""

import sys
import xml.etree.ElementTree as ET

# What OBS calls a broken build.  Same three codes osc uses, for the same
# reason: "unresolvable" means the build dependencies could not be satisfied,
# which is a packaging bug even though nothing was ever compiled.
FAILURE_CODES = ("failed", "broken", "unresolvable")

# Targets OBS deliberately did not build.  Not failures, and not proof of a
# build either.
SKIPPED_CODES = ("excluded", "disabled")

# States that mean the build is still moving.  "dirty" is the repository
# attribute OBS sets while it works out what to rebuild.
WAITING_CODES = ("blocked", "scheduled", "dispatching", "building",
                 "signing", "finished")


def main(argv):
    if len(argv) != 2:
        print("usage: obs-results.py <results.xml>", file=sys.stderr)
        return 1

    try:
        root = ET.parse(argv[1]).getroot()
    except (OSError, ET.ParseError) as exc:
        print(f"::error::could not read the OBS results: {exc}")
        return 1

    eligible = 0
    failures = []
    unsettled = []

    for result in root.findall("result"):
        repository = result.get("repository", "?")
        arch = result.get("arch", "?")
        target = f"{repository}/{arch}"

        # The repository itself can be marked dirty while the status inside it
        # still reads like the previous build.  Treat that as not settled
        # rather than reading the stale code.
        if result.get("dirty") == "true":
            unsettled.append(f"{target} (rebuilding)")
            continue

        for status in result.findall("status"):
            code = status.get("code", "")
            if code in SKIPPED_CODES:
                continue
            if code in WAITING_CODES:
                unsettled.append(f"{target} ({code})")
                continue

            eligible += 1
            if code in FAILURE_CODES:
                details = (status.findtext("details") or "").strip()
                failures.append(f"{target}: {code}"
                                + (f" ({details})" if details else ""))

    if unsettled:
        print("not settled yet: " + ", ".join(unsettled))
        return 2

    if failures:
        for failure in failures:
            print(f"::error::OBS build failed on {failure}")
        print(f"{len(failures)} of {eligible} target(s) failed")
        return 1

    if eligible == 0:
        print("::error::OBS built this package nowhere: every target is "
              "excluded or disabled, or there are no results at all. Check "
              "that the package exists and has build targets enabled.")
        return 1

    print(f"{eligible} target(s) built, none failed")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
