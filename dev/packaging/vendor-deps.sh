#!/usr/bin/env bash
#
# Populate externals/cpm/ with the source of every CPM dependency, so that a
# later cmake configure can run with the network switched off.  See GH #3388.
#
# Distribution build systems -- sbuild, mock, OBS, makepkg in a clean chroot --
# all build inside a sandbox with no network at all, while CoolProp fetches ten
# dependencies at configure time.  Running this script once turns the working
# tree into something those builders can use: cmake/dependencies.cmake finds
# externals/cpm/<name>/ on its own and skips every download.
#
# Usage:
#   dev/packaging/vendor-deps.sh [--with-tests]
#
#   --with-tests   also vendor Catch2, so that the offline tree can build and
#                  run the test suite.  Off by default, because a packaging
#                  build does not need it and it is a large checkout.
#
# The dependency list is deliberately NOT written down in this script.  It
# configures the project once with the network available and then reads back
# which packages CPM resolved and where it put them, so a dependency added to
# cmake/dependencies.cmake is picked up here without anyone having to remember
# to edit this file.

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
vendor_dir="${repo_root}/externals/cpm"

with_tests=0
for arg in "$@"; do
    case "${arg}" in
        --with-tests) with_tests=1 ;;
        -h|--help)
            sed -n '2,25p' "${BASH_SOURCE[0]}"
            exit 0
            ;;
        *)
            echo "vendor-deps.sh: unknown argument '${arg}'" >&2
            exit 2
            ;;
    esac
done

work_dir="$(mktemp -d)"
trap 'rm -rf "${work_dir}"' EXIT

echo "==> Resolving dependencies (this needs network access)"

# Point the vendoring lookup at a directory that does not exist, so that this
# configure run always resolves from upstream.  Without it, an earlier run's
# output would be re-used and a changed pin in cmake/dependencies.cmake would
# never make it into externals/cpm/.
cmake_args=(
    -S "${repo_root}"
    -B "${work_dir}/resolve"
    -DCMAKE_BUILD_TYPE=Release
    -DCOOLPROP_SHARED_LIBRARY=ON
    -DCOOLPROP_VENDORED_DEPS_DIR="${work_dir}/deliberately-absent"
    -DCOOLPROP_REQUIRE_VENDORED_DEPS=OFF
)
if [[ ${with_tests} -eq 1 ]]; then
    cmake_args+=(-DCOOLPROP_CATCH_MODULE=ON -DBUILD_TESTING=ON)
fi

cmake "${cmake_args[@]}" > "${work_dir}/resolve.log" 2>&1 || {
    echo "error: cmake configure failed while resolving dependencies" >&2
    tail -40 "${work_dir}/resolve.log" >&2
    exit 1
}

cache="${work_dir}/resolve/CMakeCache.txt"
if [[ ! -f "${cache}" ]]; then
    echo "error: cmake produced no CMakeCache.txt at ${cache}" >&2
    exit 1
fi

# CPM records every package it resolved as
#   CPM_PACKAGE_<name>_SOURCE_DIR:INTERNAL=<path>
# Emit "<name> <path>" pairs, one per line.
mapfile -t entries < <(
    sed -n -E 's/^CPM_PACKAGE_(.*)_SOURCE_DIR:INTERNAL=(.*)$/\1 \2/p' "${cache}"
)

if [[ ${#entries[@]} -eq 0 ]]; then
    echo "error: CPM reported no packages in ${cache}." >&2
    echo "       Either the configure step resolved nothing, or CPM changed how it" >&2
    echo "       records packages.  Refusing to declare the tree vendored." >&2
    exit 1
fi

echo "==> Copying ${#entries[@]} dependencies into ${vendor_dir}"
mkdir -p "${vendor_dir}"

for entry in "${entries[@]}"; do
    name="${entry%% *}"
    src="${entry#* }"

    if [[ ! -d "${src}" ]]; then
        echo "error: package '${name}' points at '${src}', which is not a directory" >&2
        exit 1
    fi

    dest="${vendor_dir}/${name}"
    rm -rf "${dest}"
    mkdir -p "${dest}"
    cp -a "${src}/." "${dest}/"
    # Drop the git metadata of the packages CPM cloned.  It is dead weight in a
    # release tarball, and a nested .git confuses git archive and dpkg-source.
    rm -rf "${dest}/.git"
    printf '    %-20s %s\n' "${name}" "$(du -sh "${dest}" | cut -f1)"
done

echo "==> Verifying that the tree now configures with no downloads"

# This is the check that matters.  COOLPROP_REQUIRE_VENDORED_DEPS makes
# cmake/dependencies.cmake abort if any package resolves from anywhere other
# than externals/cpm/, so a dependency this script missed fails here rather
# than in somebody's build chroot.
verify_args=(
    -S "${repo_root}"
    -B "${work_dir}/verify"
    -DCMAKE_BUILD_TYPE=Release
    -DCOOLPROP_SHARED_LIBRARY=ON
    -DCOOLPROP_REQUIRE_VENDORED_DEPS=ON
)
if [[ ${with_tests} -eq 1 ]]; then
    verify_args+=(-DCOOLPROP_CATCH_MODULE=ON -DBUILD_TESTING=ON)
fi

if ! cmake "${verify_args[@]}" > "${work_dir}/verify.log" 2>&1; then
    echo "error: the vendored tree still wants to download something" >&2
    tail -40 "${work_dir}/verify.log" >&2
    exit 1
fi

echo "==> Done.  ${vendor_dir} now holds every dependency."
echo "    Build offline with:"
echo "      cmake -B build -S . -DCOOLPROP_SHARED_LIBRARY=ON -DCOOLPROP_REQUIRE_VENDORED_DEPS=ON"
