# Brief: enforce the "one other SDK" rule (author decision Q2, 2026-11-11)

Book: Korean quantum-computing textbook at `/home/claude/qc-book/`. Chapter text is `chNN/chNN_<slug>.md`
(source of truth), html is rebuilt with
`cd chNN && python3 ../tools/md_to_html.py <slug-without-ext> "<title>" chNN_results.png [extra png]`
(title = the `<title>` of the existing html; ch17 also takes `ch17_checklist.png`). Each chapter also has a
Claude Docs copy that must end up identical to the md (see "Docs sync" below).

## The decision
The plan's rule is: **"다른 SDK 병기는 하나뿐이고, 왜 그 SDK인지 한 문장이 있다."** ch01–06 follow it
(section 6 opens with "이 장에서 X를 고른 이유는 …"). From ch07 on, section 6 became a 3–5 column comparison
table and some chapters carry snippets in three SDKs. The author decided: **cut back to the rule.**

Target state of every chapter you are assigned:
1. Exactly **one** SDK besides the chapter's main tool appears as code and as a table column. The main tool is
   Qiskit, except ch07 (main: QuTiP/qutip-qip) and ch17 (main: PennyLane).
2. Section 6 opens with one sentence saying why that SDK was chosen for this chapter (pattern of ch01–06:
   "이 장에서 X를 고른 이유는 …이기 때문이다."). Look at ch04 or ch06 section 6 for the tone.
3. The section-6 table keeps the main-tool column and the kept-SDK column only (for ch07, whose table has one row
   per SDK, keep the two rows). Rows that become empty/meaningless for the two remaining columns are dropped.
4. Snippets in any other SDK are removed (anywhere in the chapter, not only section 6), together with prose that
   exists only to explain them. Prose elsewhere in the chapter (sections 1, 4, 5, 7, 9 학습 노트, 10 출처와 환경)
   that refers to a removed SDK/snippet/function must be fixed or removed so nothing dangles: grep the chapter for
   the removed SDK names and API names. Version pins of removed SDKs leave "10. 출처와 환경". References and
   doc links that only supported a removed snippet are removed from section 10.
5. The section-6 heading names the kept SDK only (e.g. "## 6. 다른 SDK 대응 — PennyLane").
6. Other chapters that point at something you removed (grep all `ch*/ch*.md` for e.g. "N장 6절", the removed
   function names) are fixed too — report these; do not edit chapters outside your assignment except for such a
   one-sentence pointer fix, and tell me exactly what you changed there.

Keep mapping (decided; do not change it):

| ch | main tool | keep as the one other SDK | remove |
| --- | --- | --- | --- |
| 07 | QuTiP / qutip-qip (section 4) | Qiskit Aer (section 6 snippet) | Cirq, PennyLane, CUDA-Q rows |
| 08 | Qiskit | Cirq (section 6 snippet) | PennyLane column |
| 09 | Qiskit | PennyLane (section 6 snippet) | Cirq, CUDA-Q columns |
| 10 | Qiskit | PennyLane (section 6 snippet) | Cirq column |
| 11 | Qiskit | Cirq (section 6 snippet) | PennyLane column |
| 12 | Qiskit | QuTiP (snippet currently in section 4: independent check that integrates the Schrödinger equation without Trotter) | PennyLane snippet in section 6, PennyLane/Cirq/CUDA-Q columns. Move the QuTiP snippet into section 6 if section 4 reads fine without it; otherwise leave it in section 4 and point to it from section 6 |
| 13 | Qiskit | PennyLane (snippet currently in section 4; the differentiable-circuit tool of the chapter) | Cirq snippet in section 6, Cirq/CUDA-Q/QuTiP columns. Same move-or-point choice as ch12 |
| 14 | Qiskit + IBM Runtime | Cirq (section 6 snippet) | PennyLane, CUDA-Q, QuTiP columns |
| 15 | Qiskit Aer | Cirq (section 6 snippet) | PennyLane, CUDA-Q, "전용 도구" columns (Stim/quimb may stay as a one-line mention in prose if section 7 or 10 cites them) |
| 16 | Qiskit | PennyLane (snippet currently in section 4) | Qualtran snippet in section 6, Cirq+Qualtran/CUDA-Q/QuTiP columns. Same move-or-point choice |
| 17 | PennyLane (`qp.estimator`, section 4) | Qualtran (Cirq-based; section 6 snippet) | Qiskit, CUDA-Q, Azure QRE columns. The Qualtran snippet's comment refers to "16장 SelectPauliLCU" which no longer exists in ch16 — reword so it stands alone |

Write Korean literally (never `\uXXXX` escapes). In md running text `~ _ [ * |` are written `\~ \_ \[ \* \|`
(follow the file's existing style); a literal `|` inside a table cell breaks the table. PennyLane is imported as
`import pennylane as qp`. Snippet first lines carry "— 작성됨, 미실행" (PyPI is blocked; no SDK can be run).
Do not change any numbers, the numpy scripts, or sections you have no reason to touch.

## After editing each chapter
- Table sanity: every table row has the same cell count as its header. No garbled Hangul.
- Rebuild the html; check a distinctive phrase of your edit is in it.
- Keep a log `/tmp/claude-0/-home-claude/fd63146d-511d-5ae9-8dab-3100fb85f0db/scratchpad/q2/<chNN>_log.md`:
  what was removed, the new reason sentence, pointer fixes.

## Docs sync
Follow `/home/claude/qc-book/tools/docs_sync_brief.md` (artifact URLs, tool usage; step "Do NOT edit the local md
files" there applies to the sync step only). Extra notes for this job:
- Removing a block (code block, paragraph, table): `read` with payload {kind:"search", text:"<words in it>"} to get
  its block id, then `update` with payload {ops:[{op:"delete", target:{kind:"blocks", ids:["<id>"]}}]}.
- Replacing a whole table (columns dropped): get the table's block id and its `h` hash and the doc `rev` from a
  read, then {op:"replace", target:{kind:"blocks", ids:["<table id>"]}, with:{as:"markdown", from:{kind:"inline",
  content:"| a | b |\n| --- | --- |\n| … | … |"}}, ifHash:"<h>", ifRev:<rev>}.
- Heading text change: find/replace on the heading's plain text works.
- Verify with the ready-made differ after each chapter:
  `python3 /tmp/claude-0/-home-claude/fd63146d-511d-5ae9-8dab-3100fb85f0db/scratchpad/docdiff.py <saved export result file> <local md>`
  (an `export` result too large for the context is saved to a file automatically and its path is printed; pass
  that path). Done when it prints `hunks: 0`.
- If a Docs call fails with a rate limit (429), wait a little by doing other work, then retry; do not hammer it.

## Final message
Per chapter: kept SDK, reason sentence, what was removed (snippets/columns/prose), pointer fixes in other chapters,
html rebuilt (yes/no), Docs diff result (hunks: N). Anything you could not do and why.
