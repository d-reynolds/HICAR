# AI assistant (MCP server)

HICAR ships a companion **Model Context Protocol (MCP) server**,
[`hicar-mcp`](https://github.com/HICAR-Model/hicar-mcp), that gives AI coding
assistants — Claude Code, Cursor, and any other MCP client —
knowledge of HICAR. Once connected, you can ask your
assistant questions like:

- *"What does the `mp` namelist option do, and what are its valid values?"*
- *"Which microphysics schemes does HICAR support?"*
- *"Generate a minimal namelist using Morrison microphysics and the YSU PBL, then validate it."*
- *"What units is the `pressure_i` output variable in?"*
- *"Where is the iterative wind solver implemented?"*

The server's knowledge is extracted directly from the HICAR source (the same
`get_nml_var_metadata`, physics-scheme constants, and `get_varmeta` variable
catalog used by the [self-documenting executable](namelist_options.md)), so it
stays consistent with the model. It exposes every namelist option (group, type,
default, allowed values, range, description), the physics-scheme registry, the
full model-variable catalog, the documentation, and lexical + semantic code
search.

## Install

The recommended install is an isolated tool install (with the optional
offline semantic-search model):

```bash
pipx install "hicar-mcp[semantic]"
# or
uv tool install "hicar-mcp[semantic]"
# or, lightweight (lexical search only):
pip install hicar-mcp
```

The server ships with HICAR's knowledge **bundled in**, so it works with **zero
configuration** — you do not need a HICAR checkout or a compiled binary to look
up namelist options, schemes, variables, or docs.

## Connect your assistant

**Claude Code** (CLI):

```bash
claude mcp add hicar -- hicar-mcp serve
```

**Claude Desktop / Cursor** — add to the client's `mcpServers` config:

```json
{
  "mcpServers": {
    "hicar": { "command": "hicar-mcp", "args": ["serve"] }
  }
}
```

Then start (or restart) your client and ask it a HICAR question. In Claude Code
you can run `/mcp` to see the `hicar` server and its tools.

## Modes

The server has two modes, selected automatically:

| Mode | When | What you get |
|------|------|--------------|
| **Bundled** (default) | no configuration | Namelist / scheme / variable / docs lookup, static validation, and namelist generation — all from the bundled knowledge. |
| **Live** | `HICAR_REPO` points at a checkout | Everything above, always reflecting your working tree, **plus** source-aware tools: `code_search`, `find_symbol`, `read_source`. |

To enable live mode, pass the checkout path:

```bash
claude mcp add hicar -e HICAR_REPO=/abs/path/to/HICAR -- hicar-mcp serve
```

Optional upgrades light up automatically: build a HICAR binary (or set
`HICAR_BINARY`) and `validate_namelist`/`generate` can additionally delegate to
the model's own `--check-nml` / `--gen-nml`. Run `hicar-mcp doctor` at any time
to see what was discovered and which features are active.

## What it can do

Tools cover namelist options (`list`/`get`/`search_namelist_options`,
`validate_namelist`, `generate_namelist_tool`, `explain_namelist`), physics
schemes (`get_physics_schemes`, `resolve_physics_scheme`), model variables,
example namelists, documentation, and code (`code_search`, `find_symbol`,
`read_source`, `semantic_search`). It also provides `hicar://…` resources and a
few guided prompts (`configure_hicar_run`, `debug_namelist_error`,
`explain_physics_choice`).

## How the knowledge stays current

The bundled knowledge is regenerated from HICAR source automatically: a HICAR
release triggers the `hicar-mcp` repository to rebuild its artifacts and publish
a new version. See the [`hicar-mcp` repository](https://github.com/HICAR-Model/hicar-mcp)
for details, source, and issue tracking.
