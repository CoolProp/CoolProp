# Beads archive (read-only)

This project tracked issues with [Beads](https://github.com/steveyegge/beads)
until 2026-09-26, when the tracker moved to **Linear** (team CoolProp,
issue key `COO`).  See "Issue Tracking — Linear" in `AGENTS.md`.

`issues.jsonl` is the final export and is kept only so that old
`CoolProp-xxx` ids in commit messages, plans and code comments still
resolve:

```bash
grep '"id":"CoolProp-p8ub"' .beads/issues.jsonl | python3 -m json.tool
```

Every issue still open at the cutover was either closed with a triage
reason (see its `close_reason`) or moved to Linear; a moved issue's close
reason names its `COO-<n>` id, and the Linear issue names its source ids.

Do not create or update issues with `bd`.
