# How to issue a HICAR release

This is the maintainer checklist for publishing a tagged HICAR release. Releases are
cut from **`main`** (the trunk — see [CONTRIBUTING.md](CONTRIBUTING.md)).

> You need **admin/maintainer** rights on `HICAR-Model/HICAR` (to publish the release
> and to re-bless the regression baseline).

---

## How versioning works

HICAR's version is the git tag.
At CMake-configure time the build stamps the nearest tag into the executable:

- `GIT_TAG` = `git describe --tags --match 'v*' --always --dirty` (nearest `v*` release tag + commits-since + `-dirty`)
- `GIT_COMMIT` = `git rev-parse --short=12 HEAD`

These are compiled in (`CMakeLists.txt`, the `-DGIT_TAG=` / `-DGIT_COMMIT=` defines),
printed in the run header banner, and written into every output file's
metadata. So once you create the `vX.Y.Z` tag, the next build reports that version
automatically.

**Tag format: `vX.Y.Z`** X for big, big changes, Y for major releases, and Z for 
minor releases. This formatting was previously not followed, but from now on should be:

```
    As soon as the generals and the politicos
    can predict the motions of your mind,
    lose it. Leave it as a sign
    to mark the false trail, the way
    you didn’t go. Be like the fox
    who makes more tracks than necessary,
    some in the wrong direction.

    - Wendell Berry
```

### The next version is resolved from PR labels

`release-drafter` computes the next release text from the labels on the PRs merged since
the last release (config: `.github/release-drafter.yml`):

| Bump | PR labels |
|---|---|
| **major** (`X`) | discretionary |
| **minor** (`Y`) | `breaking`, `enhancement`, `feature`, `changed` |
| **patch** (`Z`) | `bug`, `fix`, `documentation`, `ci`, `tooling` (default) |

So **label your PRs** as they merge (this also drives the release-note categories:
Added / Changed incl. model output / Fixed / Removed-Deprecated / Docs-CI-tooling).
The resolver only *suggests* the version — you confirm or override it at publish time.

---

## The release pipeline at a glance

```
PRs merge to main
   └─> release-drafter.yml  ──> keeps a DRAFT GitHub Release up to date
                                 (categorized notes + resolved next version)
maintainer:
   1. finalize CHANGELOG.md         (hand-curated narrative)
   2. tag the release commit         (annotated: git tag -a vX.Y.Z + push)
   3. review + PUBLISH the draft  ──> published on the vX.Y.Z tag
                                         │
        ┌────────────────────────────────┼───────────────────────────────┐
   release: published              the new tag                      next build
        │                               │                                │
   notify-hicar-mcp.yml          Zenodo archives the              git describe stamps
   ─> hicar-mcp rebuilds          tagged source + mints           vX.Y.Z into the run
      its artifacts               a version DOI                   header + output files
```

---

## Pre-release checklist

-  **`main` is at the commit you want to release and is green.** Every merge into
      `main` is gated by the four required checks (`hicar-full-test`, `valgrind`,
      `gpu-check`, `snow-parity`) and the post-merge push run re-signs `main`'s HEAD —
      see [docs/ci_cd_pipeline.md §2.6](docs/ci_cd_pipeline.md). (`gpu-check` on `main`
      is a maintainer dispatch: `gh workflow run gpu.yml --ref main` if you want it
      signed on the release commit.)
- **Note: there is still no validation test in place for the model, only a regression test.
  Once a validation test is added, it should be a mandatory part of the release pipeline and
  documented** 
- **Decide the version** (`vX.Y.Z`)
- **`CHANGELOG.md` is finalized** for this version (next section).

---

## Issuing the release

### 1. Finalize the changelog

`CHANGELOG.md` is hand-curated (Keep-a-Changelog style) and complements the
auto-generated release notes with the narrative.

#### First make sure our local main is up to date

```bash
git fetch --all
git checkout main
git pull origin main
```

- rename the `## [Unreleased]` section to `## [vX.Y.Z] - YYYY-MM-DD`, copy the 
  draft release's text to this version section. Edit the draft release text now, since
  this becomes cannon after this step.
- add a fresh empty `## [Unreleased]` block above it.

### 2. Tag and release

Now commit this updated `CHANGELOG.md` as the release commit, containing 
**the version tag as the commit message**:

```bash
git commit -a -m "HICAR release vX.Y.Z"
git push origin main
git tag -a vX.Y.Z -m "vX.Y.Z"
git push origin tag vX.Y.Z
```

Confirm the version stamp resolves to your tag:

```bash
git describe --tags --match 'v*' --always --dirty     # -> vX.Y.Z
```

### 3. Review the draft release

On GitHub → **Releases**, release-drafter has a **draft** for the next version:

- Set the draft's **tag** to `vX.Y.Z` — pointing at the tag you pushed in step 2
  (confirm the **target** is the intended commit on `main`).
- Copy the same release notes from the changelog that you edited/fixed in the last step here.
- Edit the generated notes for readability; reconcile them with `CHANGELOG.md`.

### 4. Publish

Press **Publish release** (on the `vX.Y.Z` tag — yours from step 2). 
This fires the release-published automation:

- **`notify-hicar-mcp.yml`** dispatches `hicar-release` to `HICAR-Model/hicar-mcp`,
  which regenerates its bundled knowledge artifacts (namelist schema, scheme registry,
  variable catalog, docs, semantic index) from the release ref — see
  [docs/mcp.md](docs/mcp.md). *(Requires the `HICAR_MCP_DISPATCH_TOKEN` repo secret.)*
- **Zenodo** archives the tagged source and mints a **version DOI** under the concept
  DOI, if the Zenodo–GitHub integration is enabled for the repo (the DOI badge in the
  README indicates it is).

---

## Post-release

- **Re-bless the regression baseline if the release changed model output.** The
      regression lane diffs against the most recent *blessed* commit (no stored file —
      see [docs/ci_cd_pipeline.md §2.3/§3](docs/ci_cd_pipeline.md)). If this release
      intentionally changed output, bless the release commit so later regressions
      compare against it:
      run **`bless-baseline.yml`** (Actions → *Bless regression commit* →
      `workflow_dispatch`) with the release SHA and a reason. Admin/maintainer only.
- **Verify the automation.** The version + DOI badges update; the `notify-hicar-mcp`
      run succeeded and `hicar-mcp` rebuilt.
- **Announce / cite.** Update any DOI references if you cite the new version.

---

## Prerequisites & secrets (one-time, for reference)

| Thing | Why |
|---|---|
| Admin/maintainer role | publish the release; dispatch `bless-baseline.yml` |
| `HICAR_MCP_DISPATCH_TOKEN` repo secret | lets `notify-hicar-mcp.yml` trigger `hicar-mcp` (fine-grained PAT / App token with Contents r/w + Actions on `HICAR-Model/hicar-mcp`) |
| Zenodo–GitHub integration enabled | archives each published release + mints the DOI |
