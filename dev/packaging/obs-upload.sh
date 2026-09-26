#!/usr/bin/env bash
#
# Push the packaging sources to an OBS package and wait for the builds.
#
# This is the "CI pushes with osc" half of the OBS loop: a release tarball plus
# the recipe files are committed to the OBS package, OBS rebuilds every target,
# and this script exits non-zero if any of them failed.  That exit code is what
# puts the OBS outcome into the GitHub checks list.
#
# It replaces the hand upload described in dev/packaging/README.md.  There used
# to be a second mechanism, an OBS webhook that re-ran the package's source
# services on a tag, and it has been removed: it never ran anything, and a
# webhook fires before CI has built the sources it would report on.
#
# Usage:
#   dev/packaging/obs-upload.sh --tarball dist/coolprop-8.0.1dev.tar.gz \
#                               --project home:jowr --package coolprop
#
# Credentials come from the environment, never the command line, so they stay
# out of shell history and out of /proc/<pid>/cmdline:
#
#   OSC_APIURL    defaults to https://api.opensuse.org
#   OSC_USERNAME  required
#   OSC_PASSWORD  required
#
# OSC_CONFIG=/dev/null is set below and is NOT optional.  osc looks for its
# config file before it ever looks at the environment (osc/conf.py: get_config
# resolves the file, and only then applies the OSC_* overrides), so on a runner
# with no ~/.config/osc/oscrc it announces that it is about to create one and
# waits for a username on stdin.  Pointing OSC_CONFIG at /dev/null makes osc
# start from an empty config, at which point the environment credentials are
# used and nothing is written to disk.

set -euo pipefail

tarball=""
project=""
package=""
# How long to wait for OBS to notice the new sources.  This is scheduling
# latency only, not build time; the build wait below is unbounded here and
# bounded by the CI job timeout instead.
schedule_timeout_s="${OBS_SCHEDULE_TIMEOUT_S:-600}"
schedule_poll_s=5
# How long to keep re-reading the results after --watch returns, for the case
# where the watch and the next request disagree about whether anything is still
# building.  Seconds, not build time.
settle_timeout_s="${OBS_SETTLE_TIMEOUT_S:-300}"

die() {
    echo "::error::$*" >&2
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --tarball)  tarball="${2:-}"; shift 2 ;;
        --project)  project="${2:-}"; shift 2 ;;
        --package)  package="${2:-}"; shift 2 ;;
        *) die "unknown argument: $1" ;;
    esac
done

[[ -n "${tarball}" ]] || die "--tarball is required"
[[ -n "${project}" ]] || die "--project is required"
[[ -n "${package}" ]] || die "--package is required"
[[ -f "${tarball}" ]] || die "tarball not found: ${tarball}"
# Unvalidated, this reaches $(( SECONDS + ... )) and produces a confusing
# arithmetic error rather than a clear complaint about the setting.
[[ "${schedule_timeout_s}" =~ ^[0-9]+$ ]] \
    || die "OBS_SCHEDULE_TIMEOUT_S must be a whole number of seconds, got '${schedule_timeout_s}'"
[[ "${settle_timeout_s}" =~ ^[0-9]+$ ]] \
    || die "OBS_SETTLE_TIMEOUT_S must be a whole number of seconds, got '${settle_timeout_s}'"

# Credentials are checked here rather than left to osc, so a missing secret
# fails with one clear line instead of an interactive password prompt that
# would hang the job until it times out.
[[ -n "${OSC_USERNAME:-}" ]] || die "OSC_USERNAME is not set"
[[ -n "${OSC_PASSWORD:-}" ]] || die "OSC_PASSWORD is not set"
export OSC_APIURL="${OSC_APIURL:-https://api.opensuse.org}"
export OSC_CONFIG="${OSC_CONFIG:-/dev/null}"

# A plain http:// API URL does not leak the password: osc keeps allow_http off
# by default, rewrites the request to https and then dies inside urllib3 with
# "Tried to open a foreign host", which is true but tells the reader nothing.
# Checked against the pinned osc rather than assumed.  So this is about the
# error message, not about secrecy, and it is still worth one line.
case "${OSC_APIURL}" in
    https://*) ;;
    *) die "OSC_APIURL must be an https:// URL, got '${OSC_APIURL}'. osc will not talk plain HTTP to an API host, and the error it raises instead is a urllib3 traceback about a foreign host." ;;
esac

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
recipe_dir="${repo_root}/dev/packaging/obs"
# Checked here rather than at the point of use, which is after the builds have
# run: discovering then that the verdict cannot be computed would waste the
# whole wait and report a bare exit 127.
results_script="${repo_root}/dev/packaging/obs-results.py"
[[ -f "${results_script}" ]] \
    || die "missing ${results_script}, which is what decides whether the OBS builds passed"
command -v python3 > /dev/null \
    || die "python3 is needed to read the OBS results (osc is itself a Python program, so this should not normally be possible)"
tarball="$(cd "$(dirname "${tarball}")" && pwd)/$(basename "${tarball}")"

work_dir="$(mktemp -d)"
trap 'rm -rf "${work_dir}"' EXIT

osc_cmd() { osc -A "${OSC_APIURL}" "$@"; }

# The state of the builds before anything is changed.  It is the reference the
# wait below compares against, so it has to be taken before the commit.
echo "==> reading the current build state of ${project}/${package}"
results_before="$(osc_cmd results "${project}" "${package}" --xml)"

echo "==> checking out ${project}/${package}"
# --output-dir keeps the checkout inside the temp dir, so a failed run leaves
# nothing behind in the repository.
osc_cmd checkout "${project}" "${package}" --output-dir "${work_dir}/pkg"
cd "${work_dir}/pkg"

echo "==> replacing the sources"
# Everything the repository owns is cleared out first, rather than only the
# tarball.  Copying on top would add and overwrite but never delete, so a
# recipe removed or renamed upstream would live on in the OBS package: rename
# debian.libcoolprop8.install to ...9.install at a soversion bump and both
# would ship, and dh_install would consume both.
#
# Two kinds of file are left alone.  _service belongs to the package on the
# server, where it is still what an 'osc service manualrun' would use, and the
# _service:* files are generated by OBS from it and are server owned, so
# deleting them here would have the commit rejected.
shopt -s nullglob
for existing in *; do
    [[ -f "${existing}" ]] || continue
    [[ "${existing}" == "_service" || "${existing}" == _service:* ]] && continue
    rm -f -- "${existing}"
done
shopt -u nullglob

cp -- "${tarball}" .
if [[ -f "${tarball}.sha256" ]]; then
    cp -- "${tarball}.sha256" .
fi

# The recipe files are the whole reason CI does this rather than a person: a
# spec change lands on OBS in the same push as the code it goes with.
#
# Everything in the recipe directory is copied rather than a named list.  A
# named list silently stops shipping any file added later, and the Debian side
# in particular is a growing set (debian.copyright, debian.*.install, ...) where
# one missing file breaks the build for a reason that is hard to see from the
# OBS log.
recipe_count=0
for recipe in "${recipe_dir}"/*; do
    [[ -f "${recipe}" ]] || continue
    # _service is the exception, for the reason given above.
    [[ "$(basename "${recipe}")" == "_service" ]] && continue
    cp -- "${recipe}" .
    recipe_count=$(( recipe_count + 1 ))
done
# A glob that matched nothing would otherwise upload a tarball with no spec,
# which OBS accepts and then fails to build in a confusing way.
(( recipe_count > 0 )) || die "no recipe files found in ${recipe_dir}"
echo "    ${recipe_count} recipe file(s) from dev/packaging/obs/"

osc_cmd addremove

# The status is captured into a variable rather than piped into a test.  Inside
# an "if" condition set -e does not apply, and with pipefail a pipeline is
# non-zero whether grep matched nothing OR osc itself failed, so
# "if osc status | grep -q ." reads any osc error as "nothing changed" and
# carries on to report the previous build as though it were this one.  A plain
# assignment keeps set -e in force, so a broken osc stops the run.
status_out="$(osc_cmd status)"

if [[ -n "${status_out}" ]]; then
    echo "${status_out}"
    osc_cmd commit -m "CI: ${GITHUB_REPOSITORY:-CoolProp/CoolProp}@${GITHUB_SHA:-unknown}"

    echo "==> waiting for OBS to pick up the new sources (up to ${schedule_timeout_s}s)"
    # Between the commit returning and OBS reacting there is a window where the
    # backend still reports the PREVIOUS build as finished.  Going straight to
    # --watch could read that stale "succeeded" and pass.
    #
    # The wait ends when OBS shows it has noticed: either the results are
    # marked dirty, or they differ in any way from the snapshot taken before
    # the commit.  Comparing against the snapshot is what makes this work on a
    # package that already has a failing target.  An earlier version broke out
    # of this loop on seeing any failed/broken/unresolvable code, which sounds
    # reasonable and is not: those codes are almost always left over from the
    # previous build, so the loop would exit on its first poll and the run
    # would report the old result against the new commit, every time.
    #
    # Running out of time is a FAILURE, not "probably fine". A rebuild that OBS
    # never scheduled is a result nobody has checked.
    deadline=$(( SECONDS + schedule_timeout_s ))
    while true; do
        results_now="$(osc_cmd results "${project}" "${package}" --xml)"

        if grep -q 'dirty="true"' <<<"${results_now}"; then
            echo "==> OBS has marked the package dirty"
            break
        fi
        if [[ "${results_now}" != "${results_before}" ]]; then
            echo "==> the build state changed after the commit"
            break
        fi
        if (( SECONDS >= deadline )); then
            echo "${results_now}"
            die "OBS did not react to the commit within ${schedule_timeout_s}s. The sources went up, so this is not a build failure, but nothing has verified them either. Check ${OSC_APIURL/api./build.}/package/show/${project}/${package}."
        fi
        sleep "${schedule_poll_s}"
    done
else
    # Nothing changed, so OBS will not schedule anything and waiting for it to
    # react would time out. Report the state of the last build rather than
    # pretending this run proved something new.
    echo "==> sources are unchanged; reporting the existing build results"
fi

echo "==> waiting for the builds to finish"
# --watch returns once nothing is dirty and nothing is in a waiting state.
#
# NOT --fail-on-error, which looks exactly right and is a trap.  osc evaluates
# it INSIDE its polling loop and never clears it:
#
#     for results in get_package_results(...):
#         ...
#             if res['code'] in ('failed', 'broken', 'unresolvable'):
#                 failed = True
#
# (osc/core.py, get_results).  So a state left over from before this run
# uploaded anything latches a failure on the very first poll, and no later
# success can undo it.  This package sat at "broken: no source uploaded", and
# the first real run of this loop went red with every single target reporting
# "succeeded".  Checked against osc 1.27.3, the version pinned in the workflow.
#
# The verdict is therefore taken once, from the FINAL results, below.
osc_cmd results "${project}" "${package}" --watch --verbose

# Re-read the results and judge them.  --watch should already have waited, so
# the loop normally runs once; it exists because --watch returning and the next
# request being served are two different moments, and reading a result list
# that is still moving would judge stale codes.
results_file="${work_dir}/results.xml"
settle_deadline=$(( SECONDS + settle_timeout_s ))
while true; do
    osc_cmd results "${project}" "${package}" --xml > "${results_file}"

    # 0 = built and clean, 1 = a target failed or nothing built, 2 = not
    # settled.  obs-results.py explains each in its own docstring.
    verdict=0
    python3 "${results_script}" "${results_file}" || verdict=$?

    (( verdict == 2 )) || break

    if (( SECONDS >= settle_deadline )); then
        die "OBS was still reporting builds in progress ${settle_timeout_s}s after the watch returned. Check ${OSC_APIURL/api./build.}/package/show/${project}/${package}."
    fi
    sleep "${schedule_poll_s}"
done

exit "${verdict}"
