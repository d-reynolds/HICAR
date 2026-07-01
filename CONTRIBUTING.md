# Contributing to HICAR

Thanks for wanting to contribute to HICAR's development! This guide lays out some of the
steps to follow when contributing to the code base.

## Workflow

- Open an issue first (bug report or feature request) so the problem and the
  intended solution are documented before code is written.
- Make your changes on a branch and open a pull request against `main`
- Every PR is gated by the automated CI suite. See 
  [docs/ci_cd_pipeline.md](docs/ci_cd_pipeline.md) for what each
  gate checks and [docs/testing.md](docs/testing.md) for reproducing a lane
  locally before pushing.
- After review and pending a passing PR status, your code will be merged into main.

## Clean git history

A feature branch's commit history is part of the contribution — reviewers read it,
and so does anyone running `git blame` or `git bisect` later. Aim for a branch whose
commits each make one coherent change with a descriptive message, rather than a trail
of `"test fix"`, `"test fix fix"`, `"test fix 2"` checkpoints.

You don't have to commit cleanly _as you work_ — develop however you like, then **tidy
the branch before you open or update the PR**. A few techniques:

**Amend** the previous commit for a small follow-up, instead of adding a new "fix" commit:

```bash
git add <files>
git commit --amend
```

**Fix up** an earlier commit, then auto-squash every fix-up into its target in one step:

```bash
git commit --fixup <sha>            # records a "fixup! <subject>" commit
git rebase -i --autosquash main     # folds each fixup back into the commit it targets
```

A few guidelines:

- **One logical change per commit.** Ideally the build and tests pass at each commit.
  This is especially true for the history that you ultimately include as part of your PR.
- **Write specific, imperative subjects** — `"Fix halo-exchange off-by-one at tile edge"`,
  not `"fix"` — and add a body explaining the *why* when it isn't obvious from the diff.

## Test-driven Development

If you are not already among the converted, that's OK, it took us a while too.
But there is no better day to start than today. Depending on the scope and impact
of your contribution, we _really_ recommend either expanding an existing test under
`tests/` to cover it, or adding your own test. Unit tests are most preferable, but
integration and end-to-end testing can be useful for large feature additions. 

Under classical test driven development, the first step is to write a new test which
defines what your new feature/change should do, and gives a pass/fail outcome.
This test must fail at this point, since you haven't implemented anything. Your goal
now shifts from 'add feature X' to 'Make the test of feature X pass'. This helps
define your implementation goals early on, and gives a clear finish line to cross 🏁

## Documentation

If you are adding a feature, consider
adding documentation to cover this feature. If your change affects the user
experience, make sure to document this in relevant documentation section. 

To aid in keeping the documentation in sync with the code, HICAR ships with
an audimated docs audit agent that is described in [docs/developing.md](docs/developing.md).

## Branch & merge model

- HICAR follows a trunk-based model: `main` is the single long-lived branch. Create
  your `your_hot_branch` branch off the latest `main` (on your personal fork), do
  your work there, and open a PR back into `main`.
- A PR is most often integrated into `main` with a merge commit (not "Squash and merge" or
  "Rebase and merge", which rewrite the SHAs of your reviewed commits).

### Sync onto `main` before opening a PR

Rebase your feature branch onto the latest
`main` so the PR applies cleanly and your commits sit on top of current `main` (this
also surfaces any conflicts on your side, not in the PR):

```bash
# one-time: add the canonical repo as `upstream` (your fork remains `origin`)
git remote add upstream https://github.com/HICAR-Model/HICAR.git

git fetch upstream                  # get the latest main
git checkout your_hot_branch
git rebase upstream/main            # replay your commits on top of main
# if there are conflicts: edit the files, `git add <files>`, then `git rebase --continue`
git push
```

At this point, it is a good idea to run some of the tests locally before opening a PR.
You can select the [tests](docs/testing.md) which you think are relevant to your changes, 
or just run them all if you want to be really great. Then, open a new PR on github.

