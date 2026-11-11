# Brief: add the "three cost numbers" to ch10–13 (author decision Q3, 2026-11-11)

Book: Korean quantum-computing textbook at `/home/claude/qc-book/`. Chapter text `chNN/chNN_<slug>.md` is the
source of truth; `chNN/chNN_<slug>.py` is the numpy verification script (big-endian: qubit 0 is the most
significant bit) whose stdout the text quotes; html is rebuilt with
`cd chNN && python3 ../tools/md_to_html.py <slug-without-ext> "<title>" chNN_results.png` (title = `<title>` of
the existing html). No quantum SDK can be installed (PyPI blocked): everything is numpy/scipy.

## The decision
Plan decision of 2026-10-27: every circuit in part 2 records **"2큐비트 게이트 수 / 깊이 / T-count"**. Only
ch08 (`ch08/ch08_oracles.md`, paragraph starting "**게이트 비용 기록(이 책의 규칙).**") and ch09
(`ch09/ch09_grover.md`, same heading) do. The author decided: **add it to ch10, ch11, ch12, ch13.**
Amendment adopted with it: for circuits with arbitrary-angle rotations a T-count exists only relative to a
synthesis precision, so write **"회전 수 R"** and **"T ≈ 3·log₂(1/ε)·R"** (Ross & Selinger, "Optimal ancilla-free
Clifford+T approximation of z-rotations", Quantum Inf. Comput. 16, 901 (2016), arXiv:1403.2975); with
ε = 10⁻⁶ per rotation that is about 60 T per rotation.

## Counting convention (state it in one sentence per chapter and use it everywhere)
- Basis: CX (CNOT) plus arbitrary 1-qubit gates, all-to-all connectivity. Expansions: CZ = 1 two-qubit gate;
  SWAP = 3 CX; controlled-phase CP(θ) = 2 CX + 3 phase rotations (θ/2 on each qubit and −θ/2 between the CXs);
  controlled-U for a 1-qubit U = 2 CX; exp(−iθ Z⊗Z) = 2 CX + 1 Rz; exp(−iθ(XX+YY+ZZ)) = 3 CX; Toffoli/CCZ = 6 CX
  and 7 T (ch04, ch09); Fredkin = Toffoli + 2 CX.
- 2큐비트 게이트 수: number of CX/CZ after that expansion.
- 깊이: number of layers when every gate of the expanded circuit is placed as early as possible (ASAP), 1-qubit
  and 2-qubit gates each occupying one layer on the qubits they touch, measurements not counted. Check that this
  definition reproduces ch08's "깊이 ≈ 7" for BV (n = 6, s = 101101) and is compatible with ch09's "깊이 ≈ 30";
  if ch08/ch09 used a looser count, say so in your report (do not edit ch08/ch09).
- T-count: exact when every gate is Clifford or a multiple of π/4 (T = 1, Toffoli = 7); otherwise report the number
  R of rotations that are not multiples of π/4 and the estimate 3·log₂(1/ε)·R (plus the exact T from Toffolis, if any).
- Every number must come from code: build each circuit as an explicit gate list in numpy, count from the list, and
  **verify numerically that the gate list implements the same unitary/state as the chapter's circuit** (compare
  with the matrices the chapter script already builds; up to global phase) whenever the circuit has ≤ 12 qubits.
  A number you cannot back with a verified gate list must not be written; instead say in the text that the three
  numbers are not defined for that circuit and why (ch08 does this for Simon's oracle), or quote a published count
  from a source the chapter already cites, with the citation.

## What to build
1. `tools/gate_cost.py`: a small gate-list toolkit — gates as tuples (name, qubits, angle), `expand()` to CX + 1q,
   `two_qubit_count()`, `depth()` (ASAP), `t_count()` returning (exact_T, R), and `unitary(gates, n)` for
   verification. Include a self-test (`python3 tools/gate_cost.py`) that checks CP, Toffoli (6 CX, 7 T, correct
   unitary), SWAP and the BV example.
2. In each of `ch10/ch10_qpe.py`, `ch11/ch11_shor.py`, `ch12/ch12_hamiltonian.py`, `ch13/ch13_variational.py`:
   append a block at the very end, headed `# ── 게이트 비용 기록 (2큐비트 게이트 / 깊이 / T-count) ──`, that
   imports the toolkit (`sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tools"))`),
   builds the chapter's circuits as gate lists, verifies them, and prints one line per circuit. The block must not
   change any existing output line: use no global rng (or its own `np.random.default_rng(...)`), and confirm by
   diffing stdout before/after (timing lines aside) and checking the png is unchanged.
3. Which circuits: read the chapter and cover the circuits it actually builds and discusses — at least
   - ch10: textbook QPE with t = 4 (the chapter's main example), the φ = 1/3, t = 6 circuit, one round of iterative
     QPE, and the QFT on n qubits as a formula (n(n−1)/2 controlled phases → n(n−1) CX, with/without final swaps);
   - ch11: the order-finding circuit for N = 15 as the chapter constructs it. If the chapter uses a permutation
     matrix for modular multiplication, build an explicit gate version (controlled SWAP/X network for the chapter's
     a), verify it against the permutation matrix, and count that; say clearly it is this specific compiled
     construction. For N = 21 and the general case give only what you can verify or what a source cited in the
     chapter states (cite it);
   - ch12: the Trotter circuits the chapter runs (1st/2nd order, the chapter's n and step counts), and the
     energy-QPE / amplitude-estimation circuits if the chapter gives enough to define a gate list — otherwise say why not;
   - ch13: the HEA, HVA and QAOA circuits at the sizes/depths the chapter optimises.
4. In each md: add a paragraph (or a small table when there are ≥ 3 circuits) headed
   "**게이트 비용 기록(이 책의 규칙).**" at the place where the circuit's cost is discussed (section 4 or 5), in
   the style of ch08/ch09, quoting exactly the numbers your block prints, the one-sentence convention, and the
   rotation rule with its citation. Where the chapter already states partial counts (e.g. ch10 QPE t = 4 gives
   2q and depth, ch12 Trotter gives CNOT and depth, ch13 HEA gives 2q and depth), make old and new statements
   agree — if your verified count differs from the existing text, the verified count wins; report every such
   change with old/new values. Add the Ross–Selinger reference to "10. 출처와 환경" if not there. Mention the
   new block in the environment line if runtime changes noticeably. End the paragraph with
   "(규약은 책 앞의 「기호와 비용 단위」에도 적어 둔다.)".
5. Other chapters quoting a number you changed (grep `ch*/ch*.md`, especially ch14 which re-uses the part-2
   circuits with CZ counts on a line/heavy-hex layout — those are *routed* counts and need not equal yours; say so
   if a reader could be confused) — report, and fix only with a one-sentence edit.

Write Korean literally (never `\uXXXX` escapes). In md running text `~ _ [ * |` are written `\~ \_ \[ \* \|`
(follow the file's style); a literal `|` inside a table cell breaks the table. Keep a log of every md edit as
`/tmp/claude-0/-home-claude/fd63146d-511d-5ae9-8dab-3100fb85f0db/scratchpad/q3/edits.json`:
[{"chapter": "ch10", "old": "<exact old text or empty for an insertion>", "new": "<exact new text>", "after": "<for insertions: the exact text it follows>"}].

## Then, in the same chapters, apply decision Q2
Read `/home/claude/qc-book/tools/q2_one_sdk_brief.md` and apply it to ch10, ch11, ch12, ch13 (the mapping table
there says which SDK to keep). Another agent handles ch07–09 and ch14–17; do not touch those.

## Finish
Rebuild the four htmls. Sync the four Claude Docs copies as described in the Q2 brief ("Docs sync") until the
differ prints `hunks: 0` for each. Final message: per chapter the table of circuits with (2q, depth, T or R),
how each was verified, every pre-existing number you changed (old → new), the Q2 result, html and Docs status,
and anything left undone with the reason.
