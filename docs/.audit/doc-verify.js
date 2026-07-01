export const meta = {
  name: 'doc-verify',
  description: 'Audit HICAR docs against the code: flag unsupported/false statements and run documented steps',
  whenToUse: 'When you want to fact-check the documentation against the actual source and confirm the documented setup/run/domain-gen steps still work. Produces a report of flagged passages with line refs and suggested edits. Incremental by default — only re-checks docs whose text or referenced source files changed since the last clean pass.',
  phases: [
    { title: 'Plan', detail: 'discover docs + decide which need re-checking (dependency-aware incremental)' },
    { title: 'Extract', detail: 'pull verifiable claims + runnable steps from each doc' },
    { title: 'Verify', detail: 'check each claim against the source tree' },
    { title: 'Recheck', detail: 'adversarially re-verify only the flagged claims (kill false positives)' },
    { title: 'Execute', detail: 'run or static-check each documented step' },
    { title: 'Report', detail: 'write the flagged-statements + step-results report' },
    { title: 'State', detail: 'update the committed per-doc verification manifest' },
  ],
}

// ---------------------------------------------------------------------------
// Config (override per-run via the Workflow `args` object)
//   execMode : 'static' | 'cheap' | 'full'   (default 'static' — safe, never builds)
//     static -> never execute; just confirm referenced scripts/flags/targets/
//               suite-names exist and commands are well-formed.
//     cheap  -> static-check everything AND run the cheap, build-independent steps.
//     full   -> also configure+compile HICAR, conda, docker, run a real test case.
//   incremental : true  -> only re-check docs whose text OR referenced source files
//                          changed since the last clean pass (per-doc manifest).
//                 false -> check every doc (a full sweep).
//   force          : true -> force a full sweep this run even with incremental on.
//   fullSweepEvery : every Nth run, force a full sweep regardless of diffs (catches
//                    dependencies the citation heuristic missed). Default 10.
//   docsGlob   : which docs to audit (description string used by the plan agent)
//   reportPath : where the markdown report is written
//   manifestPath : the committed per-doc verification state (commit this file)
//   repo       : repo root relative to the agent cwd (defaults to '.', the working dir — portable)
// ---------------------------------------------------------------------------
// Args may arrive already parsed as an object OR — depending on how the workflow
// is invoked — as a JSON-encoded STRING. If we read `args.execMode` off a string
// it is always undefined, so execMode:'full' would silently fall through to the
// safe 'static' default and never build anything. Normalize to an object first.
const A = (typeof args === 'string')
  ? (() => { try { return JSON.parse(args) } catch (e) { return {} } })()
  : (args || {})
const REPO        = A.repo        || '.'   // agents run in the repo root; '.' keeps prompts portable
// SAFE DEFAULT: 'static' never builds/executes. A lost `args` on resume must not
// silently escalate to building the whole model. Pass execMode:'full' explicitly.
const EXEC_MODE   = A.execMode     || 'static'
const INCREMENTAL = A.incremental === undefined ? true : !!A.incremental
const FORCE       = !!A.force
const FULL_SWEEP_EVERY = A.fullSweepEvery || 10
const DOCS_GLOB   = A.docsGlob     || 'every *.md under docs/ (EXCLUDING the docs/.audit/ folder, which holds this audit tooling), plus README.md at the repo root'
const REPORT_PATH = A.reportPath   || 'docs/.audit/doc_verification_report.md'
const MANIFEST_PATH = A.manifestPath || 'docs/.audit/doc-verify-state.json'
const RUN_BATCH   = A.runBatch   || 3   // build-dependent run steps grouped per batch agent (load spread across the shared build)
const INDEP_BATCH = A.indepBatch || 5   // build-independent / static steps grouped per batch agent (fewer agents = fewer tokens)

// --------------------------------- schemas ---------------------------------
const PLAN_SCHEMA = {
  type: 'object', additionalProperties: false,
  required: ['headSha', 'fullSweep', 'reason', 'priorRunCount', 'priorExecMode', 'toCheck', 'unchanged'],
  properties: {
    headSha:       { type: 'string', description: 'git rev-parse HEAD' },
    fullSweep:     { type: 'boolean', description: 'true if every doc is being checked this run' },
    reason:        { type: 'string', description: 'why this is/ isn\'t a full sweep (no baseline, force, periodic, escalation, or incremental diff)' },
    priorRunCount: { type: 'integer' },
    priorExecMode: { type: 'string', description: 'execMode recorded in the manifest, or empty' },
    toCheck: {
      type: 'array',
      items: {
        type: 'object', additionalProperties: false,
        required: ['path', 'title', 'reason'],
        properties: { path: { type: 'string' }, title: { type: 'string' }, reason: { type: 'string' } },
      },
    },
    unchanged: {
      type: 'array',
      description: 'docs whose CONTENT is unchanged (fingerprints match) so their cached verdicts are reused. Do NOT decide here whether their steps need re-running — the orchestrator compares cachedExecMode to the requested mode to split these into execute-only vs fully-cached.',
      items: {
        type: 'object', additionalProperties: false,
        required: ['path', 'reason', 'cachedExecMode'],
        properties: {
          path: { type: 'string' },
          reason: { type: 'string' },
          cachedExecMode: { type: 'string', description: "the execMode this doc's cached verdicts were last validated at (its manifest entry's execMode; '' if none)" },
        },
      },
    },
  },
}

const EXTRACT_SCHEMA = {
  type: 'object', additionalProperties: false,
  required: ['claims', 'steps'],
  properties: {
    claims: {
      type: 'array',
      description: 'Each independently-verifiable factual assertion about HICAR (a file/function/type exists, a default, a rule, a naming convention, a numeric value, a behavior).',
      items: {
        type: 'object', additionalProperties: false,
        required: ['id', 'text', 'lineStart', 'lineEnd', 'category'],
        properties: {
          id:        { type: 'string', description: 'docpath#N, stable within this doc' },
          text:      { type: 'string', description: 'the claim, quoted/paraphrased to be checkable on its own' },
          lineStart: { type: 'integer' },
          lineEnd:   { type: 'integer' },
          category:  { type: 'string', enum: ['code-entity', 'default-or-rule', 'numeric', 'behavior', 'convention', 'other'] },
        },
      },
    },
    steps: {
      type: 'array',
      description: 'Each runnable command block a user is told to execute.',
      items: {
        type: 'object', additionalProperties: false,
        required: ['id', 'description', 'commands', 'lineStart', 'lineEnd', 'requiresBuild', 'cost', 'kind'],
        properties: {
          id:           { type: 'string' },
          description:  { type: 'string' },
          commands:     { type: 'array', items: { type: 'string' } },
          lineStart:    { type: 'integer' },
          lineEnd:      { type: 'integer' },
          requiresBuild:{ type: 'boolean', description: 'true if it needs a compiled HICAR / populated build dir' },
          cost:         { type: 'string', enum: ['cheap', 'expensive'] },
          kind:         { type: 'string', enum: ['build', 'run', 'gen', 'env', 'docker', 'check', 'config', 'misc'] },
        },
      },
    },
  },
}

const VERIFY_SCHEMA = {
  type: 'object', additionalProperties: false,
  required: ['verdicts'],
  properties: {
    verdicts: {
      type: 'array',
      items: {
        type: 'object', additionalProperties: false,
        required: ['claimId', 'status', 'evidence', 'explanation'],
        properties: {
          claimId:  { type: 'string' },
          status:   { type: 'string', enum: ['SUPPORTED', 'UNSUPPORTED', 'CONTRADICTED'] },
          evidence: {
            type: 'array',
            items: {
              type: 'object', additionalProperties: false,
              required: ['file', 'lines', 'note'],
              properties: { file: { type: 'string' }, lines: { type: 'string' }, note: { type: 'string' } },
            },
          },
          explanation: { type: 'string' },
        },
      },
    },
  },
}

const RECHECK_SCHEMA = {
  type: 'object', additionalProperties: false,
  required: ['claimId', 'flagStands', 'foundEvidence', 'finalStatus', 'reasoning', 'suggestedFix'],
  properties: {
    claimId:       { type: 'string' },
    flagStands:    { type: 'boolean', description: 'true if the doc statement really is unsupported/wrong after an honest search to PROVE it correct' },
    foundEvidence: {
      type: 'array',
      items: {
        type: 'object', additionalProperties: false,
        required: ['file', 'lines', 'note'],
        properties: { file: { type: 'string' }, lines: { type: 'string' }, note: { type: 'string' } },
      },
    },
    finalStatus:  { type: 'string', enum: ['SUPPORTED', 'UNSUPPORTED', 'CONTRADICTED'] },
    reasoning:    { type: 'string' },
    suggestedFix: { type: 'string', description: 'concrete edit to make the doc correct, or empty if flag does not stand' },
  },
}

const EXEC_SCHEMA = {
  type: 'object', additionalProperties: false,
  required: ['stepId', 'executed', 'outcome', 'commandsRun', 'evidence', 'suggestedFix'],
  properties: {
    stepId:      { type: 'string' },
    executed:    { type: 'boolean' },
    outcome:     { type: 'string', enum: ['pass', 'fail', 'static-ok', 'static-fail', 'skipped'] },
    commandsRun: { type: 'array', items: { type: 'string' } },
    evidence:    { type: 'string', description: 'key output / what proved pass or fail (trim long logs)' },
    suggestedFix:{ type: 'string', description: 'doc edit if the documented step is wrong/broken, else empty' },
  },
}

const REPORT_SCHEMA = {
  type: 'object', additionalProperties: false,
  required: ['reportPath', 'summary', 'flagCount', 'stepsExecuted', 'stepsFailed'],
  properties: {
    reportPath:    { type: 'string' },
    summary:       { type: 'string' },
    flagCount:     { type: 'integer' },
    stepsExecuted: { type: 'integer' },
    stepsFailed:   { type: 'integer' },
  },
}

const MANIFEST_WRITE_SCHEMA = {
  type: 'object', additionalProperties: false,
  required: ['written', 'manifestPath', 'docsTracked'],
  properties: {
    written:      { type: 'boolean' },
    manifestPath: { type: 'string' },
    docsTracked:  { type: 'integer' },
  },
}

// --------------------------------- prompts ---------------------------------
const CTX = `You are auditing the HICAR atmospheric model repo${REPO === '.' ? '; your current working directory IS the repository root' : ` rooted at ${REPO}`}. ` +
  `Fortran source is under src/ (main/ objects/ physics/ io/ utilities/ constants/); ` +
  `build system is CMakeLists.txt; helper scripts are under helpers/; tests under tests/; CI under .github/. ` +
  `Use Read/Grep/Glob/Bash to inspect files. Report evidence as file paths with line numbers.`

const planPrompt = () => `${CTX}

PLAN which documentation files need (re-)checking this run, and HOW. Steps:

1. headSha = \`git rev-parse HEAD\`.
2. Candidate docs = ${DOCS_GLOB}. (Repo-relative paths; never include files under build*/, external deps, .claude/, or ANYTHING under docs/.audit/ — that folder holds this audit tooling, its cache, and its report.)
3. Read the manifest at ${MANIFEST_PATH} if it exists. Shape:
   { "lastCheckedCommit": "<sha>", "runCount": <int>, "execMode": "<mode>",
     "docs": { "<path>": { "commit","docHash","refHash","execMode","flagged","refFiles":[...],"findings":[...],"stepResults":[...] }, ... } }
   (For PLANNING you need each doc's docHash/refHash/refFiles/execMode; findings/stepResults are cached analysis the report reuses — leave them alone.)
   If it is missing or unparseable, treat priorRunCount=0, priorExecMode="", and this run as a FULL SWEEP.

   For each doc that HAS a manifest entry, also note its "cachedExecMode" = that entry's execMode (the mode its cached verdicts/steps were last validated at; use "" if the doc has no entry). You do NOT act on it — just report it. The orchestrator (not you) compares cachedExecMode to the requested mode "${EXEC_MODE}" to decide whether an unchanged doc's documented steps still need to be (re-)run.

DECIDE fullSweep = TRUE if ANY of:
   - no valid manifest / no lastCheckedCommit, OR
   - force flag is set: ${FORCE}, OR
   - periodic: (priorRunCount + 1) is a multiple of ${FULL_SWEEP_EVERY}.
   (NOTE: execMode escalation does NOT force a full sweep anymore — it is handled per-doc by the execute-only path below.)

IF fullSweep: toCheck = ALL candidate docs (reason e.g. "full sweep: <why>"); unchanged = [].

ELSE (incremental) — classify EACH candidate doc as EITHER toCheck (changed/new) OR unchanged:

   A) Is the doc CONTENT-CHANGED since it was last analysed? Decide with CONTENT FINGERPRINTS, not git state (the working tree may carry uncommitted edits that are nonetheless identical to what was last audited):
      - Doc has a manifest entry WITH docHash & refHash: run
            bash docs/.audit/_hash.sh <docPath> <refFile1> <refFile2> ...
        (pass the doc path then ALL of that entry's refFiles). It prints "<docHashNow> <refHashNow>".
        UNCHANGED iff docHashNow == entry.docHash AND refHashNow == entry.refHash. Else CHANGED.
      - Doc has NO manifest entry: NEW (treat as changed).
      - Doc's entry LACKS docHash/refHash (legacy): fall back to git — CHANGED if its path or any refFile is in
            \`git diff --name-only <lastCheckedCommit> HEAD\` ∪ \`git diff --name-only HEAD\` ∪ \`git diff --name-only --cached\`.
      (You may script the fingerprint loop over all candidate docs in one bash invocation.)

   B) ROUTE each doc into EXACTLY ONE of:
      - CHANGED or NEW  ->  toCheck     (must re-extract + re-verify; reason names the doc or the changed dep file).
      - UNCHANGED       ->  unchanged[] (cached verdicts still valid; set cachedExecMode = its entry.execMode; reason e.g. "content unchanged; last validated at <cachedExecMode>"). Do NOT decide execute-only vs fully-cached — that is the orchestrator's job.

Return headSha, fullSweep, reason, priorRunCount, priorExecMode, toCheck[], unchanged[].`

const extractPrompt = (doc) => `${CTX}

Read ${doc.path} in full (with line numbers). Extract TWO lists:

1. CLAIMS — every independently-verifiable factual assertion about HICAR: a file/function/type/module exists, a default value, a configuration rule, a numeric constant, a naming convention, or a described behavior. Skip pure prose, opinions, and external links. Each claim must be phrased so it can be checked on its own, and tagged with the doc line range it came from.

2. STEPS — every runnable command block a reader is instructed to execute (bash fences, inline \`./script ...\` invocations, make targets, cmake/mpiexec/docker/conda commands). For each, record the exact commands, whether it needs a compiled HICAR or populated build/ dir (requiresBuild), a cost estimate (cheap vs expensive=multi-minute build/download/run), and a kind.

Be thorough — a missed claim is a missed bug. Return the structured object.`

const verifyPrompt = (doc, claims) => `${CTX}

These claims were extracted from ${doc.path}. For EACH, search the actual source to determine whether the code supports it. Do real searches (Grep for the named symbols/files/flags, Read the relevant code, check CMakeLists.txt and .github/ for build/CI claims).

For each claim return a verdict:
- SUPPORTED   — code clearly backs it; cite file:lines.
- UNSUPPORTED — you could not find code that backs it (the claim may be stale or aspirational).
- CONTRADICTED— the code says something different; cite the contradicting file:lines and say what it actually does.

Do not guess. "I didn't look hard enough" must not become UNSUPPORTED — search by symbol, by filename, and by behavior before concluding. Cite concrete evidence (file:lines) for EVERY verdict, including SUPPORTED ones — these citations become the doc's dependency set for incremental re-checks, so list every source file a claim depends on.

CLAIMS (JSON):
${JSON.stringify(claims)}`

const recheckPrompt = (doc, verdict) => `${CTX}

A first reviewer flagged this statement from ${doc} as ${verdict.status}. Your job is the OPPOSITE: try hard to PROVE the documentation correct. Search exhaustively — by symbol name, by filename, by synonym, across src/, CMakeLists.txt, helpers/, .github/, tests/. Renamed symbols, submodule (_h/_obj) splits, macro/cpp guards, and CMake-generated names are common reasons a real feature looks "missing".

Claim: "${verdict.claimId}" — ${verdict.explanation}
First reviewer's evidence: ${JSON.stringify(verdict.evidence)}

Decide: does the flag STAND (the doc really is unsupported/wrong) or FALL (you found supporting code — false positive)? Set flagStands accordingly, give finalStatus, cite any evidence you found, and if the flag stands write a concrete suggestedFix for the doc.`

const execPrompt = (doc, step) => {
  const policy =
    EXEC_MODE === 'static'
      ? `MODE=static: DO NOT execute anything. Statically verify the command(s): confirm every referenced script/file exists, every flag/make-target/CMake option/test-suite name is real (grep the repo / CMakeLists / Makefiles), and the invocation is well-formed. outcome = static-ok or static-fail.`
    : EXEC_MODE === 'cheap'
      ? `MODE=cheap: If the step is cheap AND build-independent (requiresBuild=${step.requiresBuild}, cost=${step.cost}) — e.g. gen_HICAR_dir.sh into a fresh temp dir, or running HICAR-tester / make check against the EXISTING ${REPO}/build — actually RUN it and report pass/fail with output. If it needs a full build, conda env, docker image, or forcing data, DO NOT run it: static-check it instead (outcome=static-ok/static-fail) and note why it was skipped.`
      : `MODE=full: Actually EXECUTE this step end-to-end and confirm it works. (Steps that need a compiled HICAR are handled separately by a build-once runner — this step does NOT need to build the model.) Use temp dirs for anything that writes outside the repo. If the step needs forcing data / Test_Cases that are absent, say so explicitly (outcome=skipped) rather than faking a pass. Capture the decisive output.`

  return `${CTX}

${policy}

Step from ${doc} (lines ${step.lineStart}-${step.lineEnd}) — ${step.description}
Commands:
${step.commands.map((c) => '    ' + c).join('\n')}

Report what you actually did, the outcome, the decisive evidence (trim long logs to the relevant lines), and — if the documented step is wrong or broken — a concrete suggestedFix for the doc.`
}

// Build/run split: ONE agent builds the exe in the SHARED working tree; multiple
// BATCH agents then run the build-dependent steps against that single build, and
// build-independent / static steps are batched too — so we neither rebuild per step
// nor spawn one agent per step.
const STEP_RESULT_ITEM = {
  type: 'object', additionalProperties: false,
  required: ['stepId', 'executed', 'outcome', 'commandsRun', 'evidence', 'suggestedFix'],
  properties: {
    stepId:      { type: 'string' },
    executed:    { type: 'boolean' },
    outcome:     { type: 'string', enum: ['pass', 'fail', 'static-ok', 'static-fail', 'skipped'] },
    commandsRun: { type: 'array', items: { type: 'string' } },
    evidence:    { type: 'string' },
    suggestedFix:{ type: 'string' },
  },
}
const BUILD_SCHEMA = {
  type: 'object', additionalProperties: false,
  required: ['buildOutcome', 'buildEvidence', 'binPath', 'results'],
  properties: {
    buildOutcome:  { type: 'string', enum: ['built', 'reused-existing', 'build-failed', 'skipped'] },
    buildEvidence: { type: 'string', description: 'how the build was produced / key build output (trim logs)' },
    binPath:       { type: 'string', description: 'path to the runnable HICAR exe other agents can invoke (e.g. ./bin/HICAR); empty if no usable build' },
    results:       { type: 'array', items: STEP_RESULT_ITEM },
  },
}
const BATCH_SCHEMA = {
  type: 'object', additionalProperties: false,
  required: ['results'],
  properties: { results: { type: 'array', items: STEP_RESULT_ITEM } },
}

// One agent's ONLY job: produce the shared HICAR exe (and validate the documented
// build step). Runs in the MAIN tree (no worktree) so other agents can run the exe.
const buildPrompt = (buildSteps) => {
  const policy = EXEC_MODE === 'full'
    ? `Build HICAR EXACTLY ONCE, in the MAIN working tree at ${REPO} (you are NOT in a throwaway worktree — the build must persist on disk for OTHER agents to run). Follow the documented build step(s) below; if none are listed, do a standard default-flags build (\`mkdir -p build && cd build && cmake ../ && make -j && make install\`). On success: buildOutcome="built", binPath = the runnable exe (e.g. ${REPO}/bin/HICAR). On failure: buildOutcome="build-failed", binPath="", put the error in buildEvidence.`
    : `Do NOT compile. Confirm the EXISTING build at ${REPO}/build and installed ${REPO}/bin/HICAR are usable: buildOutcome="reused-existing", binPath=${REPO}/bin/HICAR. If unusable/missing: buildOutcome="skipped", binPath="".`
  return `${CTX}

Your ONLY job is to make the HICAR executable available for OTHER agents to run. Do NOT run any model/test/run steps yourself.

${policy}

Also record a result for EACH documented BUILD step below (did the documented build instructions actually work as written?): stepId, executed, outcome, commandsRun, decisive evidence (trim logs), suggestedFix if the instruction is wrong/broken.

BUILD STEPS (JSON; may be empty — then just build the default way):
${JSON.stringify(buildSteps.map((s) => ({ id: s.step.id, doc: s.doc, lines: `${s.step.lineStart}-${s.step.lineEnd}`, description: s.step.description, commands: s.step.commands })))}

Return buildOutcome, buildEvidence, binPath, and results[] (one per build step).`
}

// One agent runs a BATCH of steps in order (vs one agent per step). Build-dependent
// batches reuse the already-built exe; they never compile HICAR from scratch.
const batchRunPrompt = (group, ctx) => {
  const policy = EXEC_MODE === 'static'
    ? `MODE=static: DO NOT execute anything. For EACH step, statically verify it: every referenced script/file exists, every flag/make-target/CMake option/test-suite name is real (grep repo / CMakeLists / Makefiles), invocation well-formed. outcome = static-ok or static-fail.`
    : `MODE=${EXEC_MODE}: actually RUN each step and confirm it works (outcome pass/fail), capturing decisive output.${ctx.useBuild ? ` A release build of HICAR ALREADY EXISTS at ${ctx.binPath || REPO + '/bin/HICAR'} — do NOT compile HICAR yourself to satisfy a step; invoke the existing exe. Run each model invocation from its OWN temp run dir so parallel batches don't collide (a documented \`make test_*\` target that builds its own debug variant is fine — that is the step's job).` : ''} Use temp dirs for writes outside the repo. If a step needs forcing data / Test_Cases / a conda env / docker image that is absent, mark outcome="skipped" rather than faking a pass.`
  return `${CTX}

${policy}

Work through this BATCH of ${group.length} documented step(s) IN ORDER, all inside THIS one agent (they are grouped to avoid spawning an agent per step). For EACH record: stepId, executed, outcome, commandsRun, decisive evidence (trim long logs), suggestedFix if the documented instruction is wrong/broken.

STEPS (JSON — id, source doc, line range, description, commands):
${JSON.stringify(group.map((s) => ({ id: s.step.id, doc: s.doc, lines: `${s.step.lineStart}-${s.step.lineEnd}`, description: s.step.description, commands: s.step.commands })))}

Return results[] (one per step, same order).`
}

const manifestWritePrompt = (headSha, runCount, perDocChecked, perDocExecOnly, skippedPaths) => `${CTX}

Update the per-doc verification cache/manifest at ${MANIFEST_PATH} (create it if absent). This file caches the audit so incremental runs share state across machines/CI — do NOT git-commit it yourself, just write it.

Each doc entry carries CONTENT FINGERPRINTS ("docHash","refHash") used for incremental detection. Compute them with the shared helper:
    bash docs/.audit/_hash.sh <docPath> <refFile1> <refFile2> ...   ->  prints "<docHash> <refHash>"
(pass the doc path then ALL of that doc's refFiles; this is the SAME helper the planner uses, so fingerprints stay consistent).

1. Read the existing ${MANIFEST_PATH} if present.
2. Produce the merged JSON:
   { "lastCheckedCommit": "${headSha}", "runCount": ${runCount}, "execMode": "${EXEC_MODE}",
     "docs": { "<path>": { "commit","docHash","refHash","execMode","flagged","refFiles":[...],"findings":[...],"stepResults":[...] }, ... } }

   - CHECKED docs (freshly re-audited; data below): REPLACE each entry wholesale — commit="${headSha}", execMode="${EXEC_MODE}", the given flagged / refFiles (sorted, de-duped) / findings / stepResults, and docHash+refHash computed via the helper from the doc + its refFiles.
   - EXECUTE-ONLY docs (verdicts reused, steps freshly run; data below is path + fresh stepResults): start from the doc's EXISTING manifest entry and PRESERVE its findings, refFiles, and flagged UNCHANGED (its content was identical this run); UPDATE commit="${headSha}", execMode="${EXEC_MODE}", stepResults=the provided fresh results, and RECOMPUTE docHash+refHash via the helper from the doc + its (preserved) refFiles.
   - SKIPPED docs (${JSON.stringify(skippedPaths)}): carry their EXISTING entries forward UNCHANGED (preserve commit, docHash, refHash, findings, stepResults). If somehow missing, omit them.
   - Drop entries for docs that no longer exist.
3. Write the JSON pretty-printed (2-space indent) and confirm.

CHECKED docs (full per-doc data to cache) JSON:
${JSON.stringify(perDocChecked)}

EXECUTE-ONLY docs (path + fresh stepResults; preserve their cached findings/refFiles) JSON:
${JSON.stringify(perDocExecOnly)}`

// ============================ orchestration ============================

phase('Plan')
const plan = await agent(planPrompt(), { label: 'plan-incremental', phase: 'Plan', schema: PLAN_SCHEMA })
const headSha = (plan && plan.headSha) || 'UNKNOWN'
const docs = (plan && plan.toCheck) || []
const unchanged = (plan && plan.unchanged) || []
// Deterministic split — NOT left to the planner's judgment (it kept mis-ranking
// static<full and dumping escalated docs into "skipped"). An unchanged doc whose
// cached verdicts were validated at a LOWER execMode than requested still needs its
// documented steps run at the higher mode -> execute-only; otherwise fully cached.
const MODE_RANK = { static: 0, cheap: 1, full: 2 }
const reqRank = MODE_RANK[EXEC_MODE] ?? 0
const execOnlyDocs = unchanged.filter((d) => (MODE_RANK[d.cachedExecMode] ?? 0) < reqRank)
const skipped = unchanged.filter((d) => (MODE_RANK[d.cachedExecMode] ?? 0) >= reqRank)
const newRunCount = ((plan && plan.priorRunCount) || 0) + 1
log(`Plan: ${docs.length} doc(s) to re-audit, ${execOnlyDocs.length} execute-only (verdicts reused, steps run at ${EXEC_MODE}), ${skipped.length} fully cached. ` +
    `${plan && plan.fullSweep ? 'FULL SWEEP' : 'incremental'} — ${plan ? plan.reason : ''} | execMode=${EXEC_MODE}`)

if (docs.length === 0 && execOnlyDocs.length === 0) {
  log('Nothing to do — every doc is fully cached and no steps need (re-)execution.')
  return { headSha, docsAudited: 0, execOnly: 0, skipped: skipped.length, flagCount: 0, stepsExecuted: 0, stepsFailed: 0, execMode: EXEC_MODE, incremental: INCREMENTAL }
}

// Per-doc pipeline: extract -> verify -> adversarial recheck of flagged claims.
// No barrier between stages: a doc's claims start verifying the moment they're
// extracted, and its flags get rechecked the moment they're verified.
const docResults = docs.length === 0 ? [] : await pipeline(
  docs,
  // stage 1: extract claims + steps
  (doc) => agent(extractPrompt(doc), { label: `extract:${doc.path}`, phase: 'Extract', schema: EXTRACT_SCHEMA })
    .then((ex) => ({ doc, extracted: ex })),
  // stage 2: verify all of this doc's claims against the code
  (prev) => agent(verifyPrompt(prev.doc, prev.extracted.claims), { label: `verify:${prev.doc.path}`, phase: 'Verify', schema: VERIFY_SCHEMA })
    .then((v) => ({ ...prev, verify: v })),
  // stage 3: adversarially re-check only the flagged claims for this doc
  async (prev) => {
    const verdicts = (prev.verify && prev.verify.verdicts) || []
    const claimById = new Map((prev.extracted.claims || []).map((c) => [c.id, c]))
    const flagged = verdicts.filter((v) => v.status !== 'SUPPORTED')
    const rechecks = await parallel(
      flagged.map((v) => () =>
        agent(recheckPrompt(prev.doc.path, v), { label: `recheck:${v.claimId}`, phase: 'Recheck', schema: RECHECK_SCHEMA })
          .then((rc) => {
            // enrich with the claim's line range + text so cached findings are self-contained
            const c = claimById.get(rc.claimId) || {}
            return { ...rc, claimText: c.text || '', lineStart: c.lineStart, lineEnd: c.lineEnd }
          })
      )
    )
    const confirmedFlags = rechecks
      .filter(Boolean)
      .filter((rc) => rc.flagStands)
      .map((rc) => ({ doc: prev.doc.path, ...rc }))
    // dependency set for incremental: every source file any verdict cited
    const refFiles = [...new Set(verdicts.flatMap((v) => (v.evidence || []).map((e) => e.file)).filter(Boolean))].sort()
    return { ...prev, confirmedFlags, refFiles }
  }
)

const live = docResults.filter(Boolean)
const allConfirmedFlags = live.flatMap((r) => r.confirmedFlags || [])

// Execute-only docs: their CONTENT is unchanged so we reuse their cached verdicts
// (no extract/verify/recheck of claims), but we still need their step DEFINITIONS to
// run at this (escalated) execMode. A cheap extract recovers the steps; claims ignored.
const execOnlyResults = execOnlyDocs.length
  ? (await parallel(execOnlyDocs.map((d) => () =>
      agent(extractPrompt(d), { label: `extract-steps:${d.path}`, phase: 'Extract', schema: EXTRACT_SCHEMA })
        .then((ex) => ({ doc: d, steps: (ex && ex.steps) || [] }))
    ))).filter(Boolean)
  : []
log(`Verification done: ${allConfirmedFlags.length} statement(s) flagged across ${live.length} re-audited doc(s); ` +
    `${execOnlyResults.length} doc(s) reused cached verdicts (steps still run this mode).`)

// Collect + dedup the runnable steps across ALL checked docs (barrier is justified:
// dedup over the full set avoids running the same build/test block N times).
phase('Execute')
const stepMap = new Map()
const collectSteps = (docPath, stepsArr) => {
  for (const s of (stepsArr || [])) {
    const key = JSON.stringify(s.commands)
    if (!stepMap.has(key)) stepMap.set(key, { doc: docPath, step: s })
  }
}
for (const r of live) collectSteps(r.doc.path, r.extracted.steps)
for (const r of execOnlyResults) collectSteps(r.doc.path, r.steps)
const steps = [...stepMap.values()]

// Build/run decoupling + batching:
//  - ONE build agent produces the exe in the SHARED working tree (so others can run it).
//  - build-dependent run steps fan out across BATCH agents that reuse that one build.
//  - build-independent / static steps are BATCHED (not one agent per step).
const chunk = (arr, n) => { const o = []; const k = Math.max(1, n); for (let i = 0; i < arr.length; i += k) o.push(arr.slice(i, i + k)); return o }
const ISO_KINDS = new Set(['gen', 'docker', 'env'])
const buildSteps = steps.filter((s) => s.step.kind === 'build')
const runSteps   = steps.filter((s) => s.step.kind !== 'build' && s.step.requiresBuild)
const indepSteps = steps.filter((s) => s.step.kind !== 'build' && !s.step.requiresBuild)

const execResults = []
const byStepId = new Map(steps.map((s) => [s.step.id, s]))
const absorb = (results) => {
  for (const res of (results || [])) {
    const s = byStepId.get(res.stepId)
    if (s) execResults.push({ doc: s.doc, step: s.step, res })
  }
}
const skipSteps = (list, why) => list.forEach((s) =>
  execResults.push({ doc: s.doc, step: s.step, res: { stepId: s.step.id, executed: false, outcome: 'skipped', commandsRun: [], evidence: why, suggestedFix: '' } }))

if (EXEC_MODE === 'static') {
  const batches = chunk(steps, INDEP_BATCH)
  log(`Static-checking ${steps.length} step block(s) across ${batches.length} batch agent(s).`)
  const out = await parallel(batches.map((b) => () =>
    agent(batchRunPrompt(b, { useBuild: false }), { label: 'exec-batch:static', phase: 'Execute', schema: BATCH_SCHEMA })))
  out.filter(Boolean).forEach((o) => absorb(o.results))
} else {
  // build-independent steps run concurrently with the build; isolate the mutating kinds.
  const isoIndep   = indepSteps.filter((s) => ISO_KINDS.has(s.step.kind))
  const plainIndep = indepSteps.filter((s) => !ISO_KINDS.has(s.step.kind))
  const indepThunks = [
    ...chunk(plainIndep, INDEP_BATCH).map((b) => () =>
      agent(batchRunPrompt(b, { useBuild: false }), { label: 'exec-batch:indep', phase: 'Execute', schema: BATCH_SCHEMA })),
    ...chunk(isoIndep, INDEP_BATCH).map((b) => () =>
      agent(batchRunPrompt(b, { useBuild: false }), { label: 'exec-batch:iso', phase: 'Execute', schema: BATCH_SCHEMA, isolation: 'worktree' })),
  ]

  const buildThenRun = async () => {
    let build = null
    if (buildSteps.length || runSteps.length) {
      build = await agent(buildPrompt(buildSteps), { label: 'exec:build-once', phase: 'Execute', schema: BUILD_SCHEMA })
      if (build) absorb(build.results)  // documented build step(s)' own outcomes
    }
    if (!runSteps.length) return
    const buildOK = build && (build.buildOutcome === 'built' || build.buildOutcome === 'reused-existing')
    if (!buildOK) { skipSteps(runSteps, `build ${build ? build.buildOutcome : 'unavailable'} — run steps not executed`); return }
    const binPath = build.binPath || `${REPO}/bin/HICAR`
    const groups = chunk(runSteps, RUN_BATCH)
    log(`Build ${build.buildOutcome}; running ${runSteps.length} build-dependent step(s) across ${groups.length} batch agent(s) sharing ${binPath}.`)
    const out = await parallel(groups.map((g) => () =>
      agent(batchRunPrompt(g, { useBuild: true, binPath }), { label: 'exec-batch:run', phase: 'Execute', schema: BATCH_SCHEMA })))
    out.filter(Boolean).forEach((o) => absorb(o.results))
  }

  log(`Execute (${EXEC_MODE}): 1 build agent + ${chunk(runSteps, RUN_BATCH).length} run-batch(es) sharing the exe; ${indepThunks.length} independent batch(es) concurrent.`)
  await parallel([buildThenRun, () => parallel(indepThunks).then((out) => out.filter(Boolean).forEach((o) => absorb(o.results)))])
}

const stepsExecuted = execResults.filter((e) => e.res && e.res.executed).length
const stepsFailed = execResults.filter((e) => e.res && (e.res.outcome === 'fail' || e.res.outcome === 'static-fail')).length

// group step results by the doc they came from, for the per-doc results cache
const stepsByDoc = {}
for (const e of execResults) {
  ;(stepsByDoc[e.doc] ||= []).push({
    stepId: e.step.id, lines: `${e.step.lineStart}-${e.step.lineEnd}`,
    outcome: e.res.outcome, executed: e.res.executed, suggestedFix: e.res.suggestedFix || '',
  })
}

// Synthesize the report (an agent writes the file; scripts have no fs access).
phase('Report')
const reportInput = {
  reportPath: REPORT_PATH,
  execMode: EXEC_MODE,
  incremental: INCREMENTAL && !(plan && plan.fullSweep),
  headSha,
  docsChecked: docs.map((d) => d.path),
  docsExecOnly: execOnlyDocs,
  docsSkippedCached: skipped,
  flaggedStatements: allConfirmedFlags,
  stepResults: execResults.map((e) => ({ doc: e.doc, stepId: e.step.id, lines: `${e.step.lineStart}-${e.step.lineEnd}`, ...e.res })),
}
const report = await agent(
  `${CTX}

Write a documentation-audit report to ${REPORT_PATH} (create parent dirs if needed). Use clear Markdown:

1. Summary: counts of docs freshly re-audited, docs execute-only (verdicts reused, steps run), docs fully cached/skipped, statements flagged, steps executed / failed, the execMode, and whether this was an incremental run or a full sweep (headSha ${headSha}).
2. If this was incremental, account for the docs that were NOT freshly re-audited so the report stays COMPLETE. Read ${MANIFEST_PATH} and pull each such doc's cached \`findings\` (computed at its \`commit\`):
   (a) SKIPPED (fully cached) docs — list them with their reason and include their cached \`findings\` AND cached \`stepResults\`, each marked "(cached — not re-verified this run; from commit <commit>)".
   (b) EXECUTE-ONLY docs — their verdicts were REUSED because their content was unchanged. Include their cached \`findings\` marked "(verdicts reused — content unchanged since commit <commit>)", BUT their documented-step results were FRESHLY EXECUTED this run and are in the DATA below (stepResults) — use those fresh results for them, NOT cached step results.
3. "Flagged statements" — grouped by doc file. For each: the doc line range, the claim, why it is wrong/unsupported (with the supporting code evidence file:lines), and the concrete suggested edit. These are edit targets, NOT auto-applied.
4. "Documented step results" — a table: doc, step (lines), outcome, and any suggested fix for broken instructions.
5. Note that flagged items were confirmed by an adversarial second reviewer (false positives already removed).

Do not modify any doc file other than writing the report. After writing, return the path and the headline counts.

DATA (JSON):
${JSON.stringify(reportInput)}`,
  { label: 'write-report', phase: 'Report', schema: REPORT_SCHEMA }
)

// Update the committed per-doc manifest so the next run can skip unchanged docs.
phase('State')
const perDocChecked = live.map((r) => ({
  path: r.doc.path,
  flagged: (r.confirmedFlags || []).length > 0,
  refFiles: r.refFiles || [],
  findings: r.confirmedFlags || [],     // cached flagged statements (self-contained)
  stepResults: stepsByDoc[r.doc.path] || [],  // cached documented-step outcomes
}))
// Execute-only docs keep their cached findings/refFiles; only their step outcomes are fresh.
const perDocExecOnly = execOnlyResults.map((r) => ({
  path: r.doc.path,
  stepResults: stepsByDoc[r.doc.path] || [],
}))
const manifest = await agent(
  manifestWritePrompt(headSha, newRunCount, perDocChecked, perDocExecOnly, skipped.map((s) => s.path)),
  { label: 'update-manifest', phase: 'State', schema: MANIFEST_WRITE_SCHEMA }
)

return {
  reportPath: (report && report.reportPath) || REPORT_PATH,
  manifestPath: (manifest && manifest.manifestPath) || MANIFEST_PATH,
  headSha,
  fullSweep: !!(plan && plan.fullSweep),
  docsAudited: docs.length,
  execOnly: execOnlyResults.length,
  skipped: skipped.length,
  flagCount: allConfirmedFlags.length,
  stepsExecuted,
  stepsFailed,
  execMode: EXEC_MODE,
  incremental: INCREMENTAL,
}
