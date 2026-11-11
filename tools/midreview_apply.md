# Mid-review: apply accepted findings to ch01–ch17 (md, py, html)

You are applying reviewer findings to a Korean quantum-computing textbook. Four reviewer reports are in
`/tmp/claude-0/-home-claude/fd63146d-511d-5ae9-8dab-3100fb85f0db/scratchpad/mid/`:
`R1_findings.md`, `R2_findings.md`, `R3_findings.md`, `R4_findings.md`. Each finding gives file:line, the current text,
evidence, and paste-ready Korean replacement text. Read all four reports in full first.

Book root `/home/claude/qc-book/`; chapter text is `chNN/chNN_<slug>.md` (source of truth), script `chNN/chNN_<slug>.py`,
html is built with `cd chNN && python3 ../tools/md_to_html.py <slug-without-ext> "<title>" chNN_results.png [extra png]`
(title = the `<title>` of the existing html; ch17 also takes `ch17_checklist.png` as the extra png). Text is Korean;
write Korean literally (never unicode escapes). In md, `|` inside a table cell breaks the table, and `*`, `_`, `~`, `[`
are escaped as `\*`, `\_`, `\~`, `\[` in running text (follow the surrounding style of the file).

## What to apply

Apply every finding below unless you find, by your own check, that the proposed text is wrong — then fix the proposal
and say so in your report. Locate text by exact string (line numbers may have shifted); every replacement must match
exactly once. Keep a log.

### R1 (exercises ch01–09)
A1 (use the R3-A2 wording instead, see below), B1, B2 (both places), B3, B4 (both places), B5, B6 (both md places and the
two print labels in `ch04_universality.py`), B7, C1–C10 as proposed. C11: add one sentence to the problem statement of
ch05 문제 3 defining a channel by Kraus operators (ρ → Σ KₖρKₖ†, Σ Kₖ†Kₖ = I, "추적 보존"). C12: skip (backlog).

### R2 (exercises ch10–17, ch17 second review)
A1: all four places (ch11 md lines ~13, ~29, ~241 and the print string in `ch11_shor.py`). Additionally add to
`ch11_shor.py` a short block that enumerates, for N = 15, 21, 33, 57, 77, the fraction of a coprime to N with r even and
a^{r/2} ≢ −1 (mod N), prints it, and re-run the script to confirm the numbers quoted in the replacement text (15: 6/8,
others exactly 1/2); if the script's stdout is quoted elsewhere in the md, keep md and stdout consistent.
B1 as proposed. B2: do NOT add R2's note — instead unify ch17 to the ch10/ch12 definition as in R4-B2:
γ = |⟨ψ₀|φ⟩|² is a probability everywhere in ch17; repetitions 1/γ; with amplitude amplification 1/√γ. Check every
occurrence of γ in ch17 (text, tables, exercises, code comments) and make each consistent; list them in your report.
B3: all three parts. C1 (both parts), C2 (formula grouping only; leave the attribution), C3, C4, C5, C6 (all three
places 294, 48, 202 if the same wording occurs), C7, C8, C9, C10 as proposed.

### R3 (citations)
A1 (both places), A2 (this wording for ch09:232, replacing R1-A1), A3 (three parts incl. the new §10 entry), A4 (both),
B1, C1–C11 as proposed (C11: link text and URL), C12 (mechanical replace of `https://docs.quantum.ibm.com/` by
`https://quantum.cloud.ibm.com/docs/en/` in all chapters), C13 (add the missing bibliography entries to the §10 list of
the respective chapter, in that list's existing style; for "Cerezo 2021" write 2021a/2021b in text and §10).
B2 is already done (do not touch ch15:59).

### R4 (plan compliance / coherence)
- A1: ch12:32 as proposed. Also add to ch12 §10 the reference
  `- Pfeuty, "The one-dimensional Ising model with a transverse field", *Annals of Physics* 57, 79 (1970).`
- B1, B3 (table rows), B9's ch15:59 and ch15:110, C9: ALREADY DONE in ch15 and ch17 — do not redo. ch17 §5(A) now has
  a caveat paragraph and a second table; leave them.
- B2: see R2-B2 above.
- B4: the ch16:321 replacement, and the sentence to add after the ch17 §5(C) break-even table, as proposed.
- B5: apply the whole sentence table (ch04, ch07 ×2, ch08, ch09 ×2, ch10 ×2, ch13 ×2). Do not touch the plan.
- B6: add the proposed "게이트 비용 기록" paragraph to ch08 (end of §4 or where circuit costs are discussed). Nothing else.
- B7: before each first code cell that does not parse (ch06, ch10, ch12, ch13, ch14, ch15, ch16 — lines in R4) insert the
  sentence "아래는 `chNN_<slug>.py`에서 핵심만 간추린 발췌다(그대로 실행되지 않는다). 실행 가능한 전체 코드는 스크립트 파일에 있다."
  with the right file name. If a label like "numpy 구현(실행 완료)" precedes the cell keep it (it refers to the script).
- B8: skip (needs the author's decision).
- B9 remaining items — use these texts:
  - ch08:226 "Cirq의 `_decompose_` 패턴이 이 책 3부에서 가장 유용하다:" → "Cirq의 `_decompose_` 패턴은 오라클의 비용을 셀 때 유용하다:"
  - ch07:254 "…회로 전체의 오류 예산 1차 추정이고, 14장에서 \"이 회로가 돌기는 할까\"를 가늠하는 공식이다." →
    "…회로 전체의 결어긋남 몫의 1차 추정이다. 14장의 회로들은 짧아서 이 항보다 게이트 오류와 판독 오류가 지배한다(14장 5절)."
  - ch06:211 and ch11:227: leave (chapter 18 keeps those promises).
- B10: only the ch12:144 sentence. B11: the note at ch15:254 as proposed (the table cell is already done).
- B12: append "(공개 보정 데이터에서 흔히 보이는 범위이며, 이 책이 직접 잰 값은 아니다.)" after the first of the four
  empirical statements in ch14 (line ~19).
- C1: add seeds where the call is a plain constructor: `StatevectorSampler()` → `StatevectorSampler(seed=7)`;
  `cirq.Simulator()` → `cirq.Simulator(seed=7)`; `cirq.DensityMatrixSimulator(` add `seed=7` if no seed; before
  `cudaq.sample(` add `cudaq.set_random_seed(7)` once per snippet; Aer: add `seed_simulator=7` to the `run(`/`AerSimulator(`
  call. If a snippet prints "# → value" comments that depend on the sample, leave the comments.
- C2: `cuda-q==0.16` → `cudaq==0.16.x` everywhere. Nothing else.
- C11: ch14:33 "위 기여 표" → "5절의 기여 표"; ch12:227 "게이트 수가 O(αt + …)" → "오라클 호출 수가 O(αt + …)";
  ch08:236 "7장의 T2가 이 책에서 가장 중요한 하드웨어 숫자인 이유." → "깊은 위상 되돌림 회로에서 7장의 T2를 가장 먼저 봐야 하는 이유(14장의 얕은 회로에서는 게이트·판독 오류가 더 크다).";
  ch13:9 "오류 정정 전에는 돌릴 수 없다" → "의미 있는 크기에서는 오류 정정 전에 돌릴 수 없다(14장의 3비트 QPE 같은 작은 예는 예외)".
- C12: the pointer sentence in ch12:11 as proposed.
- Skip: C3, C4, C5, C6, C7, C8, C10 (backlog / author decisions).

## After editing
1. For ch04 and ch11 (scripts changed) re-run the script (`MPLBACKEND=Agg python3 <script>` inside the chapter dir) and
   check stdout still agrees with the md.
2. Rebuild the html of every chapter whose md changed and check each html contains a distinctive phrase of one edit.
3. Scan all 17 md files for garbled Hangul introduced by your edits (rare syllables) and for broken tables
   (a table row whose cell count differs from its header).
4. Write the log to
   `/tmp/claude-0/-home-claude/fd63146d-511d-5ae9-8dab-3100fb85f0db/scratchpad/mid/apply_log.json`: a list of
   {"chapter": "chNN", "finding": "R2-A1", "old": "<exact old text>", "new": "<exact new text>"} for every md edit
   (this log drives the later sync of the Claude Docs copies, so old/new must be the exact strings).
5. Final message: counts per chapter, anything you changed relative to the reviewers' proposals and why, anything
   you could not apply. Do not edit ch18, the plan, or the tools directory. Do not touch Claude Docs.
