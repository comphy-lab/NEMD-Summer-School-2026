# Vatsal Workflow

This fork is Vatsal's working copy for local testing and small course-material
edits. The official summer-school material remains the upstream NEMD repository;
do not point participant-facing instructions at this fork unless the organisers
explicitly decide to use it.

## Repositories

- Official course repository:
  `https://github.com/Non-Equilibrium-Molecular-Dynamics/NEMD-Summer-School-2026`
- Vatsal working fork:
  `git@github.com:comphy-lab/NEMD-Summer-School-2026.git`
- Local checkout:
  `/Users/vatsal/cowork-os/0-Projects/NEMD-2026-local/NEMD-Summer-School-2026`
- Cirrus checkout:
  `/work/tc075/tc075/vatsalsy/NEMD-Summer-School-2026`

## Edit Locally, Run On Cirrus

Make durable source edits on the Mac first:

```bash
cd /Users/vatsal/cowork-os/0-Projects/NEMD-2026-local/NEMD-Summer-School-2026
git status --short --branch
```

Commit and push to the working fork:

```bash
git add <files>
git commit -m "<focused message>"
git push origin main
```

Then pull and run on Cirrus:

```bash
ssh -Y cirrus-tc075
cdnemd
git pull --ff-only
nemd-python
```

Run scripts from the Cirrus `/work` checkout. Generated plots, temporary test
files, and large run outputs belong on Cirrus under `/work`, not in the local
Git repository unless they are intentionally small example artefacts.

## Remote Editing Rule

Do not make durable source edits directly on Cirrus. Temporary diagnostics are
fine, but any useful remote change should be reproduced locally before it is
committed.

Before pulling on Cirrus, check for local remote-only files:

```bash
cdnemd
git status --short --branch
```

Untracked generated files can usually stay in place. If they block a pull,
move them to a scratch/output folder under `/work/tc075/tc075/vatsalsy`.
