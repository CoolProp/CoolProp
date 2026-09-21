#!/usr/bin/env bash
#
# Build the offline source tarball that every downstream packaging recipe
# consumes.  See GH #3388.
#
# The tarball is a plain export of the tree at a given git revision, plus:
#   - externals/cpm/, holding every CPM dependency (see vendor-deps.sh), so the
#     build needs no network;
#   - dev/gitrevision.txt, because dev/generate_headers.py falls back to that
#     file when it cannot run git, which is the case in every build chroot.
#
# Usage:
#   dev/packaging/make-release-tarball.sh [--ref <git-ref>] [--output-dir <dir>]
#
# Defaults: --ref HEAD, --output-dir dist/
#
# Output:
#   <output-dir>/coolprop-<version>.tar.gz
#   <output-dir>/coolprop-<version>.tar.gz.sha256
#
# The tarball is byte-for-byte reproducible for a given revision: entries are
# sorted, ownership is zeroed and timestamps come from the commit date, so
# rebuilding it does not churn the checksum that OBS and the .spec pin.

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
ref="HEAD"
output_dir="${repo_root}/dist"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --ref)        ref="$2"; shift 2 ;;
        --output-dir) output_dir="$2"; shift 2 ;;
        -h|--help)    sed -n '2,25p' "${BASH_SOURCE[0]}"; exit 0 ;;
        *)            echo "make-release-tarball.sh: unknown argument '$1'" >&2; exit 2 ;;
    esac
done

# Read one COOLPROP_VERSION_* component out of CMakeLists.txt.
#
# The 2014 packaging script did this with `cut -d " " -f 3` and started
# returning an empty string the day cmake-format removed the space in
# "set (COOLPROP_VERSION_MAJOR 5)".  Nobody noticed for years, because it
# failed silently.  Hence the regex tolerates both spellings AND the result is
# asserted to be a number before it is used anywhere.
read_version_component() {
    local name="$1" value
    value="$(sed -n -E \
        "s/^[[:space:]]*set[[:space:]]*\\([[:space:]]*${name}[[:space:]]+([0-9]+)[[:space:]]*\\).*/\\1/p" \
        "${repo_root}/CMakeLists.txt" | head -1)"
    if [[ ! "${value}" =~ ^[0-9]+$ ]]; then
        echo "error: could not read ${name} from CMakeLists.txt (got '${value}')." >&2
        echo "       Fix the parser here rather than shipping a mis-named tarball; GH #3388." >&2
        exit 1
    fi
    printf '%s' "${value}"
}

version_major="$(read_version_component COOLPROP_VERSION_MAJOR)"
version_minor="$(read_version_component COOLPROP_VERSION_MINOR)"
version_patch="$(read_version_component COOLPROP_VERSION_PATCH)"

# COOLPROP_VERSION_REVISION is "dev" between releases and empty on a release
# tag.  It is part of the library file name (libCoolProp.so.8.0.1dev), so it is
# part of the tarball name too, otherwise a snapshot tarball would claim to be
# the release.  Letters and digits only: it goes into an RPM Version tag.
version_revision="$(sed -n -E \
    's/^[[:space:]]*set[[:space:]]*\([[:space:]]*COOLPROP_VERSION_REVISION[[:space:]]*([A-Za-z0-9]*)[[:space:]]*\).*/\1/p' \
    "${repo_root}/CMakeLists.txt" | head -1)"

version="${version_major}.${version_minor}.${version_patch}${version_revision}"
if [[ -n "${version_revision}" ]]; then
    echo "note: this is a ${version_revision} snapshot, not a release tarball" >&2
fi

commit="$(git -C "${repo_root}" rev-parse "${ref}")"
commit_epoch="$(git -C "${repo_root}" show -s --format=%ct "${commit}")"

name="coolprop-${version}"
work_dir="$(mktemp -d)"
trap 'rm -rf "${work_dir}"' EXIT
stage="${work_dir}/${name}"

echo "==> Exporting ${ref} (${commit:0:12}) as ${name}"
mkdir -p "${stage}"
git -C "${repo_root}" archive --format=tar "${commit}" | tar -x -C "${stage}"

# generate_headers.py shells out to git for the revision and falls back to this
# file when there is no repository, which is the situation in every build
# chroot.
echo "${commit}" > "${stage}/dev/gitrevision.txt"

echo "==> Vendoring dependencies into the export"
"${stage}/dev/packaging/vendor-deps.sh"

echo "==> Creating the tarball"
mkdir -p "${output_dir}"
tarball="${output_dir}/${name}.tar.gz"

# --sort, --mtime, --owner/--group and --numeric-owner together make the
# archive reproducible; gzip -n keeps the timestamp out of the gzip header.
tar --create \
    --directory "${work_dir}" \
    --sort=name \
    --mtime="@${commit_epoch}" \
    --owner=0 --group=0 --numeric-owner \
    --file - \
    "${name}" \
  | gzip -9 -n > "${tarball}"

( cd "${output_dir}" && sha256sum "${name}.tar.gz" > "${name}.tar.gz.sha256" )

echo "==> ${tarball}"
echo "    $(cat "${output_dir}/${name}.tar.gz.sha256")"
echo "    $(du -h "${tarball}" | cut -f1)"
