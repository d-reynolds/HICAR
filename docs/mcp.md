# AI assistant (MCP server)

HICAR ships a companion **Model Context Protocol (MCP) server**,
[`hicar-mcp`](https://github.com/HICAR-Model/hicar-mcp), that gives AI coding
assistants — Claude Code, Cursor, any other MCP client —
knowledge of HICAR. Once connected, you can ask your
assistant questions like:

- *"What does the `mp` namelist option do, and what are its valid values?"*
- *"What input data do I need to run the model?"*
- *"Generate a minimal namelist using Morrison microphysics and the YSU PBL, then validate it."*
- *"What units is the `pressure_i` output variable in?"*
- *"Where is the iterative wind solver implemented?"*

The server's knowledge is extracted directly from the HICAR source, so it
stays consistent with the model. Obviously, LLM output is not the final word (i.e. they can be wrong),
and for any critical information you should consult the documentation and code yourself.
Asking an assistant where to start looking is a great use of the technology.

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

The server ships with HICAR's knowledge **bundled in**, so it works with zero
configuration — you do not need a HICAR checkout or a compiled binary to look
up namelist options, schemes, variables, or docs. That said, connecting it to
an [existing HICAR repo](#modes) will make it more effective.

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
| **Live** | `HICAR_REPO` points at a checkout | Everything above, always reflecting your working tree, **plus** source-aware tools. |

To enable live mode, pass the checkout path:

```bash
claude mcp add hicar -e HICAR_REPO=/abs/path/to/HICAR -- hicar-mcp serve
```

Optional upgrades light up automatically: build a HICAR binary (or set
`HICAR_BINARY`) and `validate_namelist`/`generate` can additionally delegate to
the model's own `--check-nml` / `--gen-nml`. Run `hicar-mcp doctor` at any time
to see what was discovered and which features are active.

## How the knowledge stays current

The bundled knowledge is regenerated from HICAR source automatically: a HICAR
release triggers the `hicar-mcp` repository to rebuild its artifacts and publish
a new version. See the [`hicar-mcp` repository](https://github.com/HICAR-Model/hicar-mcp)
for details, source, and issue tracking.
