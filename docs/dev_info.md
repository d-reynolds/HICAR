# Developer Information

Here is an assortment of concepts useful to be familiar with if you are working
with the model's source code. If you are planning to contribute to the development,
make sure you have read the [contributing](../CONTRIBUTING.md) section first. For 
a general map of the code, see [Code overview](code_overview.md).

## Contribution workflow

- Make changes on a feature branch off `main`, and open pull requests against the **`main`** branch (trunk-based model).
- Every PR is gated by the automated CI suite (build, unit/invariant tests,
  smoke runs, decomposition/restart reproducibility, and bit-for-bit
  regression). See the [CI/CD pipeline](ci_cd_pipeline.md) for what each gate
  checks and [Testing](testing.md) for how to reproduce a lane locally before
  pushing.
- For substantial additions, please get in touch first so we can help fit the
  change into the model's structure.

## Adding a new output variable

HICAR's variable handling is **metadata-driven**: once a variable is registered
and described, allocation and output happen automatically — there is no
per-variable code to add in the domain or output objects. To add a new variable, three things are needed:

1. **Give it an index.** Add an entry for the variable in the `kVARS` structure
   in `src/constants/icar_constants.F90`. Every model variable is referred to
   throughout the code by its `kVARS%name` index.

2. **Describe it.** Add a branch for the new `kVARS` index to the `get_varmeta`
   function in `src/io/default_output_metadata.F90`. This is where you set the
   output `name`, the `dimensions` (e.g. `three_d_t_dimensions`), the `units`
   and other netCDF `attributes`, and the valid `minval`/`maxval` range. The
   dimensions you set here are what drives allocation, so this step is required —
   a `kVARS` entry with no metadata branch returns an empty name and is silently
   skipped.

3. **Request it.** Have the relevant physics module ask for the variable in its
   `*_var_request` routine (e.g. `adv_std_var_request` in
   `src/physics/advect.F90`, `lsm_var_request` in `lsm_driver.F90`). Use
   `options%alloc_vars([...])` to have it allocated and `options%restart_vars([...])`
   to have it carried through restarts and written to output:

   ```fortran
   call options%alloc_vars(   [kVARS%my_new_var, ...] )
   call options%restart_vars( [kVARS%my_new_var, ...] )
   ```

That's all. The domain object's `create_variables` loops over the requested
variables, reads each one's metadata, and allocates the correctly-shaped array;
the output path writes whatever is in the output set. (`vars_to_allocate`
controls which `kVARS` are allocated; `vars_for_output` controls which are
written.)

If the variable you want is *already* registered, simply request it (step 3) in
the module that needs it — or, for output, it is exposed through the namelist's
output-variable selection.

## Domain decomposition conventions

When writing physics or I/O code, note the index conventions:

- `ims:ime` — memory bounds (include halos)
- `its:ite` — tile/interior bounds (this rank's computed cells)
- `ids:ide` — global domain bounds

Halo (ghost-cell) data is exchanged via the domain's halo methods. Operations
that need up-to-date neighbor data must run *after* the relevant halo exchange.
Be especially careful with anything whose result depends on the local subdomain
extent (local reductions used as global parameters, search windows sized from
local bounds) — these break bit-for-bit reproducibility across different MPI rank
counts.

## Building the documentation

These docs are built with [MkDocs](https://www.mkdocs.org/). To preview them
locally:

```bash
pip3 install mkdocs
mkdocs serve          # or: python3 -m mkdocs serve
```

Then open <http://127.0.0.1:8000> in a browser. The site layout (page order and
titles) is defined in `mkdocs.yml`. Use `mkdocs build --strict` to catch broken
intra-doc links before committing.

## Auditing the documentation

The documentation makes many claims about the code — file and routine
names, namelist option names, defaults, and runnable commands — and these drift
as the model evolves. The hidden `docs/.audit/` folder holds an automated
**documentation audit** that keeps the docs honest: it reads every page, checks
each factual claim against the current source, and exercises the commands the
docs tell you to run.

The audit is a Claude Code workflow (`docs/.audit/doc-verify.js`). It runs a number of agents, so
it runs **from a Claude Code session in the repo, not as a stand-alone script**
(`node doc-verify.js` will not work). Either ask in plain language — *"run the
docs audit"* — or invoke the workflow tool directly:

```text
Workflow({ scriptPath: "docs/.audit/doc-verify.js", args: { execMode: "static" } })
```

Per page it extracts the verifiable claims and the runnable command blocks,
verifies each claim against `src/` / `CMakeLists.txt` / `helpers/` / `.github/`
(citing `file:line` evidence), has a second agent adversarially re-check every
flag to drop false positives, checks (and, in the heavier modes, runs) the
documented steps, and writes a report of what it flagged with a suggested edit
for each. **It only flags problems; it never edits the docs for you.**

### **This agent, especially when run in 'full' mode, will burn your tokens. Thank you for your sacrifice/You have been warned.**

### Modes

`execMode` controls how far it goes when checking the *runnable* steps:

| `execMode` | What it does | Cost |
|---|---|---|
| `static` *(default)* | Never builds or runs anything. Confirms referenced scripts, flags, make targets, and suite names exist and that commands are well-formed. | cheap, safe |
| `cheap` | Also runs the build-independent steps (e.g. `gen_HICAR_dir.sh`, or `HICAR-tester` / `make check` against an existing `build/`). | moderate |
| `full` | Also runs example commands that require building HICAR | slow, heavy |

`static` is the default on purpose — it is safe to run anytime and never triggers
a multi-hour build. Reach for `full` deliberately.

### Incremental runs

The audit is incremental. After each run it writes a cache,
`docs/.audit/doc-verify-state.json`, recording for every page the commit it was
checked at, the source files its claims depend on, and its findings. The next run
only re-audits a page whose text *or* one of its referenced source files has
changed since that commit (with a periodic full sweep as a backstop); unchanged
pages are served from the cache. **Commit `doc-verify-state.json`** together with
your documentation changes so the cache is shared across machines and CI.

### Output

The report is written to `docs/.audit/doc_verification_report.md`, containing a summary, the
flagged statements grouped by page (line numbers, the offending claim, the
contradicting code, and a suggested fix), and a table of documented-step results.
Work through it and apply the edits you agree with.