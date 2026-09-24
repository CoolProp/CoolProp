#!/usr/bin/env bash
#
# Push the packaging sources to an OBS package and wait for the builds.
#
# This is the "CI pushes with osc" half of the OBS loop: a release tarball plus
# the recipe files are committed to the OBS package, OBS rebuilds every target,
# and this script exits non-zero if any of them failed.  That exit code is what
# puts the OBS outcome into the GitHub checks list.
#
# It replaces the hand upload described in dev/packaging/README.md.  The OBS
# webhook is a different mechanism and is NOT involved here; see the comment in
# .obs/workflows.yml for why the two cannot be combined.
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
# Two kinds of file are left alone.  _service belongs to the release path, and
# the _service:* files are generated by OBS from it and are server owned, so
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
# --watch returns once nothing is dirty or in a waiting state, and
# --fail-on-error exits 1 when any result is failed, broken or unresolvable
# (osc/core.py, get_results).  "excluded" and "disabled" are not failures,
# which is what we want: a repository that does not build a given architecture
# is a configuration choice, not a regression.  set -e carries that exit code
# out of the script.
osc_cmd results "${project}" "${package}" --watch --fail-on-error --verbose

# --fail-on-error decides from the rows it iterated, so it also exits 0 when
# there were no rows at all: a package with no repositories configured, or a
# mistyped package name, would otherwise look like a clean build.
results_final="$(osc_cmd results "${project}" "${package}" --xml)"
# Only targets that could actually build are counted.  A row that says
# "excluded" or "disabled" is a target OBS deliberately did not build, and
# --fail-on-error skips those too, so a package whose targets are ALL excluded
# or disabled produces rows, reports no failure, and would otherwise be a green
# tick for a package that built nowhere.
#
# The walk is over <status> tags rather than the whole document because a
# <result> element carries a code attribute of its own: a repository can be
# marked excluded while the status inside it is fine, and matching the raw text
# would subtract that one twice.
#
# awk rather than "grep -c ... || echo 0": grep -c already prints 0 when it
# matches nothing and exits 1 as well, so the fallback appends a second 0 and
# the arithmetic below then fails on "0\n0".  awk counts and exits 0 either way.
result_rows="$(awk '
    {
        s = $0
        while (match(s, /<status [^>]*>/)) {
            tag = substr(s, RSTART, RLENGTH)
            if (tag !~ /code="(excluded|disabled)"/) n++
            s = substr(s, RSTART + RLENGTH)
        }
    }
    END { print n+0 }' <<<"${results_final}")"
(( result_rows > 0 )) \
    || die "OBS built ${project}/${package} nowhere: every target is excluded, disabled, or there are no results at all. Check that the package exists and has build targets enabled."
echo "==> ${result_rows} target(s) built, none failed"
