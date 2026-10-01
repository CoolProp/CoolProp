# Agent Instructions

## AI Assistance and Attribution

**AI assistance is permitted in CoolProp. Attribution is required.**

Using an AI coding assistant — Claude, GitHub Copilot, Cursor, ChatGPT, or
anything comparable — to help write code, documentation, tests, or data in this
repository is allowed and needs no prior approval. What is not optional is
saying so.

### When attribution is required

Attribute whenever an AI tool **generated or materially shaped** something that
lands in the repository: source code, comments, documentation, test cases, fluid
data, build/CI configuration, or the commit and PR text itself.

Attribution is **not** required when the assistant never touched what you
committed — editor autocomplete of a symbol name, asking a model to explain
existing code, or using one to search the codebase or the literature.

The line to apply: *would a reviewer reading this diff want to know an AI wrote
it?* When in doubt, attribute. Over-attribution costs nothing; undisclosed AI
authorship is the thing this policy exists to prevent.

### How to attribute

**1. Commit trailer (required).** Every commit whose content an AI **generated
or materially shaped** carries a `Co-authored-by:` trailer naming the tool, in
the trailer block at the end of the commit message after a blank line:

```
fix(HEOS): guard against negative density in the PT solve

Co-authored-by: Claude Opus 5 <noreply@anthropic.com>
```

Name the model or tool. Including the version is encouraged but not mandatory —
`Claude Opus 5`, `Claude`, `GitHub Copilot`, and `Cursor` are all acceptable
identities; an unattributed AI-written commit is not. Use the tool's documented
noreply address where it has one (Anthropic's is `noreply@anthropic.com`). If
several tools contributed, use one trailer line per tool.

This is already the de facto convention in this repository — `git log` shows
`Co-authored-by: Claude ...` on AI-assisted merges going back many releases.
This section makes it a rule rather than a habit.

**2. PR description (required).** In the pull request body, state briefly what
the AI did and what you verified yourself. The **Verification Process** section
of the [PR template](.github/PULL_REQUEST_TEMPLATE.md) is the natural place for
it. For example:

> Claude Opus 5 drafted the Chebyshev fitting loop and its unit tests. I derived
> the recurrence by hand, checked the fitted coefficients against REFPROP at 40
> states across the domain, and rewrote the error handling.

A commit trailer is invisible in the GitHub review UI; the prose note is what
the reviewer actually sees.

### What attribution does not change

Attribution is a disclosure, not a disclaimer. The human contributor is the
author of record and remains fully responsible for the change:

- **Correctness.** You must understand the code you submit and be able to
  explain and defend it in review. "The model wrote it" is not an answer to a
  review question.
- **Testing.** AI-generated code is held to the same standard as any other: it
  needs tests, and it must pass the project's quality gates.
- **Provenance and licensing.** By submitting, you assert the contribution is
  yours to give under CoolProp's MIT license. Verify that generated code is not
  a verbatim reproduction of incompatibly licensed source, and that any
  correlation or model taken from the literature is properly cited in
  `CoolPropBibTeXLibrary.bib`.
- **Scientific claims.** Thermophysical property models, fitted coefficients,
  and numerical results must be validated against a trusted reference — never
  accepted because a model produced them confidently.

Maintainers may request changes on, or close, a pull request that appears to be
unreviewed AI output. Review capacity is the scarce resource here; generated
volume that a human has not read first spends it without adding value.

## Issue Tracking

This project uses **bd** (beads) for issue tracking. Run `bd prime` for full workflow context.

### Quick Reference

```bash
bd ready              # Find available work
bd show <id>          # View issue details
bd update <id> --claim  # Claim work atomically
bd close <id>         # Complete work
bd dolt push          # Push beads data to remote
```

## Non-Interactive Shell Commands

**ALWAYS use non-interactive flags** with file operations to avoid hanging on confirmation prompts.

Shell commands like `cp`, `mv`, and `rm` may be aliased to include `-i` (interactive) mode on some systems, causing the agent to hang indefinitely waiting for y/n input.

**Use these forms instead:**
```bash
# Force overwrite without prompting
cp -f source dest           # NOT: cp source dest
mv -f source dest           # NOT: mv source dest
rm -f file                  # NOT: rm file

# For recursive operations
rm -rf directory            # NOT: rm -r directory
cp -rf source dest          # NOT: cp -r source dest
```

**Other commands that may prompt:**
- `scp` - use `-o BatchMode=yes` for non-interactive
- `ssh` - use `-o BatchMode=yes` to fail instead of prompting
- `apt-get` - use `-y` flag
- `brew` - use `HOMEBREW_NO_AUTO_UPDATE=1` env var

<!-- BEGIN BEADS INTEGRATION v:1 profile:minimal hash:ca08a54f -->
## Beads Issue Tracker

This project uses **bd (beads)** for issue tracking. Run `bd prime` to see full workflow context and commands.

### Quick Reference

```bash
bd ready              # Find available work
bd show <id>          # View issue details
bd update <id> --claim  # Claim work
bd close <id>         # Complete work
```

### Rules

- Use `bd` for ALL task tracking — do NOT use TodoWrite, TaskCreate, or markdown TODO lists
- Run `bd prime` for detailed command reference and session close protocol
- Use `bd remember` for persistent knowledge — do NOT use MEMORY.md files

## Session Completion

**When ending a work session**, you MUST complete ALL steps below. Work is NOT complete until `git push` succeeds.

**MANDATORY WORKFLOW:**

1. **File issues for remaining work** - Create issues for anything that needs follow-up
2. **Run quality gates** (if code changed) - Tests, linters, builds
3. **Update issue status** - Close finished work, update in-progress items
4. **PUSH TO REMOTE** - This is MANDATORY:
   ```bash
   git pull --rebase
   bd dolt push
   git push
   git status  # MUST show "up to date with origin"
   ```
5. **Clean up** - Clear stashes, prune remote branches
6. **Verify** - All changes committed AND pushed
7. **Hand off** - Provide context for next session

**CRITICAL RULES:**
- Work is NOT complete until `git push` succeeds
- NEVER stop before pushing - that leaves work stranded locally
- NEVER say "ready to push when you are" - YOU must push
- If push fails, resolve and retry until it succeeds
<!-- END BEADS INTEGRATION -->
