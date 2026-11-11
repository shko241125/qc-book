# Brief: sync Claude Docs chapters to the corrected markdown files

Goal: each chapter's Claude Docs artifact must contain the same text as the corrected local `chNN_*.md`
(the md is the source of truth; review fixes were applied to the md files today).

Artifacts (project links):
- ch01 https://claude.ai/artifact/V2EyxZP6bjzs3zWdpzafnh  ↔ /home/claude/qc-book/ch01/ch01_qubit_statevector.md
- ch02 https://claude.ai/artifact/ETe286oNy5PBg6oAbW5pWs  ↔ /home/claude/qc-book/ch02/ch02_single_qubit_gates.md
- ch03 https://claude.ai/artifact/QTdrWfsU69DNv2WCGnGJf8  ↔ /home/claude/qc-book/ch03/ch03_multi_qubit.md
- ch04 https://claude.ai/artifact/Qmg9NNsH2LBfDmP6QZBBAJ  ↔ /home/claude/qc-book/ch04/ch04_universality.md
- ch05 https://claude.ai/artifact/Ax7UxFHNcr6grSBrk4dP3X  ↔ /home/claude/qc-book/ch05/ch05_measurement.md
- ch06 https://claude.ai/artifact/W5tTcZQknc48hQv44hUgL4  ↔ /home/claude/qc-book/ch06/ch06_circuit_library.md
- ch07 https://claude.ai/artifact/RdGqLrLD4rRFaWSzUWFLRb  ↔ /home/claude/qc-book/ch07/ch07_physical_layer.md
- ch08 https://claude.ai/artifact/FrV8bjWVFLCvW5P1GNW328  ↔ /home/claude/qc-book/ch08/ch08_oracles.md
- ch09 https://claude.ai/artifact/XRVuuvTdJfqhP84mge2Qwt  ↔ /home/claude/qc-book/ch09/ch09_grover.md
- ch10 https://claude.ai/artifact/9qLtAgPpzjiZp73ZDNPtGL  ↔ /home/claude/qc-book/ch10/ch10_qpe.md
- ch11 https://claude.ai/artifact/7ZSmA5WJGEkn8gsk7VhXnZ  ↔ /home/claude/qc-book/ch11/ch11_shor.md
- ch12 https://claude.ai/artifact/TPLepDA6G5GB9qh9KWMny9  ↔ /home/claude/qc-book/ch12/ch12_hamiltonian.md
- ch13 https://claude.ai/artifact/B3WKfoNib8vnNfQhamB54m  ↔ /home/claude/qc-book/ch13/ch13_variational.md
- ch14 https://claude.ai/artifact/6sqcF6tj1h6CdJuyCW6hUK  ↔ /home/claude/qc-book/ch14/ch14_hardware.md
- ch15 https://claude.ai/artifact/MQEBfjDhrj8oEoZi8SzRyL  ↔ /home/claude/qc-book/ch15/ch15_sources.md
- ch16 https://claude.ai/artifact/CyJgVgPwZSrCXXyqUBh13R  ↔ /home/claude/qc-book/ch16/ch16_primitives.md
- ch17 https://claude.ai/artifact/3PUkFPWS3DYPAn5dtXnfg7  ↔ /home/claude/qc-book/ch17/ch17_discovery.md
- ch18 https://claude.ai/artifact/3ktJRgQyxwPSc8h43GcPs5  ↔ /home/claude/qc-book/ch18/ch18_search_tools.md
- ch19 https://claude.ai/artifact/UcvZWQi3VLAZDxfQ8r9EcS  ↔ /home/claude/qc-book/ch19/ch19_candidates.md

Tip: decode the export's base64 with a python script that reads it from a file you save, never by hand-transcribing
(earlier runs mis-transcribed Korean characters). Put the Korean replacement text literally in the update payload.

Procedure per chapter (Claude_Docs MCP tools; load them with ToolSearch "select:mcp__Claude_Docs__read,mcp__Claude_Docs__export,mcp__Claude_Docs__update" first, and call
`mcp__Claude_Docs__guide` with items ["topic.index","topic.editing"] once before any write):
1. `read` with ref {object:"project", id:"<artifact id from the URL>"} → gives container project id (artifact.slug), the file id (files[0].id) and the node id (files[0].content.id).
2. `export` (container {kind:"project", id:<slug>}, file <file id>, format "markdown") → base64 → decode to a temp file in /tmp/claude-0/.../scratchpad (write a small python script; the base64 is in the tool result — save it via a heredoc).
3. Diff exported markdown vs the local md (python difflib on paragraphs). Ignore: the byline line (doc has a date chip/mention; md has "2026-10-29 · Gosh"), the `![...](chNN_results.png)` figure line (md only), and escaping differences (`\_`, `\[`, `\|`, `\~`).
4. For every differing paragraph/table cell/code block, update the doc so it matches the md:
   - prose: `update` with payload {ops:[{op:"replace", target:{kind:"find", text:"<old plain text, NO markdown markers like ** or backticks>"}, with:{from:{kind:"inline", content:"<new text>"}, as:"text"}}]}. `find` matches plain text (first match); keep the old text short but unique (a distinctive phrase inside the changed sentence), and replace just that span.
   - a whole code block or table cell: find its block id with `read` payload {kind:"search", text:"<distinctive words in it>"} then `update` with target {kind:"blocks", ids:[<id>]} and with {as:"markdown", from:{kind:"inline", content:"```python\n...\n```"}} (include ifRev from the last read/update result).
   - Table cells: avoid `|` characters inside the replacement text (they truncate the cell); rewrite without pipes.
5. Re-export and re-diff; stop when the only differences are the ignorable ones. Report per chapter: number of spans changed, anything you could not sync (and why).

Do NOT edit the local md files. Do not change doc titles or the byline. Work in English; the doc content is Korean — copy it verbatim from the md.
