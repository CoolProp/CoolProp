# Shared git hooks

This directory holds tracked git hooks for the CoolProp repository.

To enable them, run once per local clone:

```sh
git config core.hooksPath .beads-hooks
```

After that, every commit/push runs the hooks below.

## Installed hooks

- **prepare-commit-msg, pre-commit, pre-push, post-merge, post-checkout** — beads
  integration shims (managed by `bd hooks install --shared`). They auto-export
  `.beads/issues.jsonl` and trigger bd's audit/sync logic on commit/push/pull.
  See https://github.com/steveyegge/beads/.

- **commit-msg** — auto-appends `[skip ci]` to a commit message when the only
  staged paths are under `.beads/`. GitHub Actions skips workflows when the
  commit message contains a `[skip ci]` (or `[ci skip]`, etc.) marker, so a
  beads-only sync commit no longer triggers the ~11 workflows that fire on
  every push to master. Mixed commits (beads + code) are unaffected.

## Updating the beads shims

If a future `bd` version updates its hook shims, re-run:

```sh
bd hooks install --shared --force
```

That regenerates the bd-managed sections inside the existing files; everything
outside the `BEGIN BEADS INTEGRATION` markers (including the `commit-msg` hook
above) is preserved.
