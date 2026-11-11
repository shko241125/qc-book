# R9 findings — gate-cost records (Q3) and one-SDK cut (Q2) in ch10–13

Reviewer R9, 2026-11-10. Scope: today's edits to `ch10_qpe`, `ch11_shor`, `ch12_hamiltonian`, `ch13_variational` (md + appended script blocks), `tools/gate_cost.py`, checked against `ch00_notation.md` §3 and the older records in ch08/ch09.

Headline: **every number in the four new tables and "읽는 법" paragraphs was reproduced by an independent counter** (own gate representation, own expansion, own ASAP depth, own matrix-based Clifford/T/rotation classifier, own simulator, own dense references). No grade-A finding. Two grade-B findings (one unsourced/understated claim, one formula with a missing domain) and eight grade-C items.

Counts: **A 0 · B 2 · C 8**

Scratch code and outputs: `/tmp/claude-0/-home-claude/fd63146d-511d-5ae9-8dab-3100fb85f0db/scratchpad/mid/R9/` (`r9count.py`, `r9_ch10.py` … `r9_ch13.py`, `r9_toolcheck.py`, `r9_beauregard.py`, `out_ch10.txt` … `out_ch13.txt`, re-runs of the chapter scripts in `run/` and `run_strict/`).

Line numbers refer to the files as they stood at 17:40 UTC (ch12 md after the `2^{t\_b}` notation change, which is not reported here).

---

## Grade B

### [B1] ch11 §4 and §7 — `ch11/ch11_shor.md:99` and `:220`
Current text (line 99): `N = 15(n = 4)에서 11큐비트·수천 게이트, N = 21에서 13큐비트.`
Current text (line 220): `그런 회로는 N = 15에서도 수천 게이트라 현재 하드웨어에서 노이즈에 묻힌다.`

What is wrong.
- "11큐비트" and "13큐비트" are backed: Beauregard 2003 (quant-ph/0205095) §3 says "The circuit for factorization of N=15 uses eleven qubits", and 2n + 3 = 13 for n = 5.
- "수천 게이트" is **not** in the cited source. The paper gives orders only ("2n+3 qubits, O(n³ k_max) gates and a depth of O(n³)"), no count for N = 15. Checked at https://ar5iv.labs.arxiv.org/html/quant-ph/0205095 (abstract and §3).
- The sentence sits inside a paragraph headed "게이트 비용 기록(이 책의 규칙)", so a reader takes "게이트" in this book's unit. In that unit the figure is too low for the full circuit. I built the uncompiled Beauregard controlled multiplier for N = 15 on 2n + 3 = 11 qubits (Fourier adder, doubly controlled modular adder of Fig. 5, CMULT, controlled SWAP, inverse CMULT), verified it on all 2 × 15 basis inputs plus a control-superposition phase test, and counted with the book's rule (extra rule I had to add: doubly controlled phase = 3 CP + 2 CX; exact QFT without SWAPs; zero-angle phases dropped):

  | circuit (N = 15, uncompiled Beauregard) | 2-qubit | depth | T-count |
  | --- | --- | --- | --- |
  | controlled ×7, one multiplication | 1752 | 2148 | T 708 + R 1546 |
  | controlled ×4 (or ×1, ×2), one multiplication | 1512 | ≈ 1790 | T 708 + R 1276 |
  | all 2n = 8 controlled multiplications (a = 7) | 12336 | 14589 | T 5664 + R 10478 |
  | only the two with a^(2^j) ≠ 1 (needs 7⁴ ≡ 1, i.e. the answer) | 3264 | 3939 | T 1416 + R 2822 |

  So "수천" holds for one multiplication or for a circuit that already uses the answer, and is about 3–4× too low for the honest full circuit (≈ 1.2×10⁴ CX; by arithmetic only, not simulated, still ≈ 1.0×10⁴ if the doubly controlled phase is done with 6 CX instead of 8). A structure-only lower bound needs no simulation: 2n controlled-U × 2n modular adders = 4n² = 64 adders, each with four (n+1)-qubit QFTs of 10 CP → 2,560 CP = 5,120 CX from the adders' QFTs alone.
- The book's own rule (Q3 brief) is that a number without a verified gate list in the chapter script is not written. My count lives in reviewer scratch code (`R9/r9_beauregard.py`), so the primary replacement below uses only the cited qubit count and the structural arithmetic.

Proposed replacement, line 99 (replace the quoted sentence only):

> N = 15(n = 4)에서 11큐비트, N = 21에서 13큐비트다(2n + 3). 게이트 수는 논문이 차수만 주고, 이 장도 그 회로를 짜서 세지 않았다. 구조만 따지면 제어 곱셈 2n = 8번에 모듈러 덧셈기가 4n² = 64개 들어가고, 덧셈기마다 5큐비트 QFT를 네 번 쓰므로 QFT의 제어 위상만 2,560개, 이 책의 규칙으로 CX 5,120개다. 나머지 게이트는 여기에 더해진다.

Optional extra sentence, only if the editor ports the reviewer's construction into `ch11_shor.py` and re-verifies it there:

> 실제로 짜서 세면(이중 제어 위상을 CP 3개 + CX 2개로 풀 때) 제어 ×7 한 번이 CX 1,752개·깊이 2,148이고, 여덟 번을 모두 넣은 회로는 CX 12,336개다.

Proposed replacement, line 220 (replace the quoted clause only):

> 그런 회로는 N = 15에서도 2큐비트 게이트가 5천 개를 넘어(4절 게이트 비용 기록) 현재 하드웨어에서 노이즈에 묻힌다.

### [B2] ch13 §4 — `ch13/ch13_variational.md:100`
Current text: `깊이는 사슬 01·12·23이 큐비트를 공유해 계단으로 놓이므로 3L + 2다(일반 n은 n + 3L − 2)`

What is wrong. The general formula fails for n = 2. With one CX per layer there is no staircase: the circuit is Ry, CX, Ry, CX, …, Ry and the depth is 2L + 1 (my counter: L = 1, 2, 3 → 3, 5, 7), while n + 3L − 2 = 3L gives 3, 6, 9. For n ≥ 3 the formula is right (checked n = 3, 4, 5, 6, 8, 12 and L = 1, 2, 3, 5). The script's own assert covers only n ∈ {3, 4, 5, 6, 8}, so the gap is not caught there. Graded B, not A, because the chapter's circuits are all n = 4 and every tabulated value is correct.

Proposed replacement:

> 깊이는 사슬 01·12·23이 큐비트를 공유해 계단으로 놓이므로 3L + 2다(n ≥ 3인 일반 n은 n + 3L − 2, n = 2는 2L + 1)

---

## Grade C

### [C1] toolkit — `tools/gate_cost.py:204` and `:206`
Current code: `np.allclose(A, ph * B, atol=tol)` (and `np.allclose(A, B, atol=tol)` two lines above), with `tol=1e-9`.

`np.allclose` keeps its default `rtol=1e-5`, so the effective tolerance is about 1e-5 relative, not the 1e-9 the signature suggests. Evidence: `gc.same_up_to_phase(diag(1,1,1,e^{0.7i}), diag(1,1,1,e^{0.7i+1e-6i}))` returns True (`R9/r9_toolcheck.py`). The verifications are still real: a wrong gate list gives O(1) differences. I patched a scratch copy to `rtol=0` and re-ran all four chapter scripts: every assert passes and stdout is byte-identical (`R9/run_strict/`). So the fix is safe.

Proposed code (both lines): `np.allclose(A, ph * B, rtol=0, atol=tol)` and `np.allclose(A, B, rtol=0, atol=tol)`.

### [C2] ch12 §5 — `ch12/ch12_hamiltonian.md:156`
Current text: `2차가 1차보다 같은 스텝에서 10배 정확하고 게이트는 2배 미만 — 항상 2차를 쓴다.`

The new §4 table shows the ratio depends on merging: as written it is exactly 2× (1 step 6 → 12; 32 steps 192 → 384), 1.68× after merging adjacent halves (322), and nearly 1× after grouping the ZZ block (198 vs 192). "2배 미만" is true only for the merged variants, which this sentence does not say.

Proposed replacement:

> 2차가 1차보다 같은 스텝에서 10배 정확하다. 게이트는 그대로 두면 2배지만 이웃한 반쪽을 합치면 1.7배(10r + 2), ZZ 묶음까지 합치면 1차와 거의 같다(6r + 6 대 6r, 4절 표) — 항상 2차를 쓴다.

### [C3] ch11 §7 — `ch11/ch11_shor.md:220`
Current text: `×7 mod 15가 SWAP 세 개로 줄어드는 것은 답을 알기 때문이다`

The new §4 record and the snippet text both give SWAP 3 + X 4 (the X layer is what makes it ×7 and not ×8, and it is why y = 0 and 15 are exchanged). Pre-existing sentence, now inconsistent with the paragraph added today.

Proposed replacement:

> ×7 mod 15가 SWAP 세 개와 X 네 개로 줄어드는 것은 답을 알기 때문이다

### [C4] convention page vs chapters — `ch00/ch00_notation.md:45` and `:46`
Current text (line 45, rule cell): `환산: 제어 위상 2, exp(−iθZ⊗Z) 2, SWAP 3, 일반 2큐비트 게이트 3, Toffoli 6. 14장의 수치는 장치 연결에 맞춰 SWAP을 넣은 뒤의 CZ 수라 더 크다`
Current text (line 46, rule cell): `1큐비트 게이트도 한 층을 차지한다. 측정은 세지 않는다`

No contradiction with ch08–13, but the page omits rules that ch10–13 state and use: a controlled phase of angle π is one CZ (ch10:100; ch12 "제어 Z 1"), controlled rotation = 2 CX (ch12:88), CCZ = 6 CX (ch12:88), Fredkin = 8 CX (ch11:89), and "회전은 합치지 않는다" (ch10:100, while ch12 lists merged variants as separate rows). ch14 counts the controlled Z as 2 (CZ 12 for the 3-bit QPE before routing), which only ch10:112 points out. The depth of Toffoli/Fredkin rows also depends on the gate order inside the Toffoli, which is fixed in ch04:149 but named nowhere in the cost rule.

Proposed replacement, line 45 rule cell:

> 기준은 CX와 임의의 1큐비트 게이트, 모든 큐비트가 이어진 장치. 환산: 제어 위상 2(각도가 π인 제어 Z는 CZ 1개), 제어 회전 2, exp(−iθZ⊗Z) 2, SWAP 3, 일반 2큐비트 게이트 3, Toffoli·CCZ 6, Fredkin 8. 14장의 수치는 제어 Z도 2개로 풀고 장치 연결에 맞춰 SWAP을 넣은 뒤의 CZ 수라 더 크다

Proposed replacement, line 46 rule cell:

> 1큐비트 게이트도 한 층을 차지한다. 측정은 세지 않는다. 이어 붙은 회전은 합치지 않고 세며, 합친 값을 적을 때는 따로 밝힌다(12장). Toffoli는 4장의 15게이트 순서로 편다

### [C5] ch12 §4 — `ch12/ch12_hamiltonian.md:103`
Current text: `s = 8(확률 0.495·0.185, ‖U\_s − U‖ = 0.0034)에서는 40964개·깊이 49136 — "분해능 × Trotter 스텝 × 항 수"의 곱이 그대로 보인다.`

The table's "T ≈" uses ε = 10⁻⁶ per rotation for every row. ch00 line 47 now says the error accumulates to Rε. For the last QPE row R = 28679, so the whole-circuit bound is Rε ≈ 0.03 (about 1.1×10⁻⁴ per controlled U, still below the Trotter error 0.0034 per U, so the row is not wrong). One pointer sentence keeps the big-R rows from being read as "T for a fixed total accuracy".

Proposed replacement (append one sentence):

> s = 8(확률 0.495·0.185, ‖U\_s − U‖ = 0.0034)에서는 40964개·깊이 49136 — "분해능 × Trotter 스텝 × 항 수"의 곱이 그대로 보인다. 표의 T 어림은 회전 하나의 오차를 ε = 10⁻⁶으로 둔 값이어서, 회전이 28679개인 이 행에서는 회로 전체로 Rε ≈ 0.03까지 쌓일 수 있다(책 앞 「기호와 비용 단위」).

### [C6] ch13 §3 vs §6 — `ch13/ch13_variational.md:42` (pre-existing; now next to today's reason sentence)
Current text (line 42): `이것이 "미분 가능 회로"의 정확한 뜻이고 PennyLane의 기본 미분 모드다.`
Line 216 (moved into §6 today) says: `` `qp.grad`는 시뮬레이터(default.qubit)에서는 역전파(backprop), 실기에서는 파라미터 이동을 자동 선택한다 ``

The two sentences disagree about the default. Line 216 matches PennyLane's `diff_method="best"` behaviour as I know it; I could not run PennyLane here, so treat the direction of the fix as "make them agree".

Proposed replacement, line 42:

> 이것이 "미분 가능 회로"의 정확한 뜻이고 PennyLane이 실기에서 쓰는 미분 모드다(시뮬레이터의 기본은 역전파, 6절).

### [C7] wording of the breakdowns — `ch10/ch10_qpe.md:112`, `ch12/ch12_hamiltonian.md:103`
Current text (ch10): `제어 Z 1 + 제어 S 2 + 제어 T 2 + 역 QFT의 CP 6개 12 + SWAP 2개 6.`
Current text (ch12): `제어 Q 한 번 = 제어 Z 1 + 제어 Ry 2개 4 + CCZ 6 = 2큐비트 11`

Arithmetic is right (23 and 11), but "CP 6개 12", "SWAP 2개 6", "제어 Ry 2개 4" put a count and a cost side by side with no marker.

Proposed replacement (ch10):

> 제어 Z 1 + 제어 S 2 + 제어 T 2 + 역 QFT의 CP 6개(CX 12) + SWAP 2개(CX 6).

Proposed replacement (ch12):

> 제어 Q 한 번 = 제어 Z 1 + 제어 Ry 2개(CX 4) + CCZ(CX 6) = 2큐비트 11

### [C8] ch11 §4 — `ch11/ch11_shor.md:135`
Current text: `제어 버전(Fredkin 3 + CNOT 4)으로도 CZ 28(위 게이트 비용 기록)`

The record it points to counts CX, not CZ. The number is the same in either basis (one CX = one CZ plus two H), but the unit should match the table or say why it differs.

Proposed replacement:

> 제어 버전(Fredkin 3 + CNOT 4)으로도 CX 28(위 게이트 비용 기록; CZ 기저로 옮겨도 28)

---

## Summary

- All tabulated 2-qubit counts, depths, exact T-counts, rotation counts R and "T ≈" values in ch10–13 match my independent recount, and match the scripts' printed lines. The stated formulas hold: QPE 2t + t(t−1) + 3⌊t/2⌋, QFT n(n−1) (+3⌊n/2⌋), HEA 3L + 2 (n = 4), HVA/QAOA 7p + 1, Trotter 12r / 10r + 2 / 6r + 6, energy QPE 5100s + 164.
- The one substantive problem is the uncited "수천 게이트" for the uncompiled N = 15 circuit (B1): the qubit count is sourced, the gate count is not, and under this book's rule the full circuit is about 1.2×10⁴ CX.
- The convention is applied the same way in all four chapters; ch00 §3 does not contradict it but omits several rules the chapters rely on (C4).
- The one-SDK cut is clean in all four chapters.

## Checked and found correct

1. **Independent recount** (own counter, never importing `gate_cost.py`): ch10 QPE t = 3, 4 (U = T) and t = 4, 6, 8 (φ = 1/3), QFT n = 3, 4, 6, 8 with and without SWAPs, IPE rounds; ch11 ×7 and ×4 networks, controlled ×7 / ×4, order finding t = 4, 8, 9; ch12 first/second-order steps, 32-step and 24-step circuits in as-is / merged / brick / grouped variants, controlled U_s, energy QPE s = 1 and 8 (12 qubits), amplitude estimation t_b = 3, 5, 7; ch13 HEA L = 1–4 (CZ and CX), HVA p = 1–4, QAOA p = 1–3. Zero disagreements with md or stdout.
2. **Depth robustness**: unchanged when the two +θ/2 phases of CP are placed after the CX pair instead of before, when the eigenstate-preparation X is dropped, and when the Fredkin is mirrored.
3. **Convention consistency**: CP = 2 CX + 3 phase rotations everywhere; CP(π) = one CZ and CP(0) dropped wherever they occur (ch10, ch12); π/4 multiples counted as exact T, other angles as R, in all four chapters; no rotation merging except in the explicitly labelled ch12 rows; ch08 (BV depth 7) and ch09 (Grover 12 / 26 / T 14, full 24 / 53 / T 28) reproduce under the same rule.
4. **Verifications in the script blocks are real**: full-unitary or full-state comparison against the chapter's own matrices and functions, global phase only. The toolkit's simulator, expansion, depth and classification agree with mine on 300 random circuits over every gate name (one apparent mismatch was my own classifier's looser tolerance near 3π/4, not a tool error).
5. **Re-runs** from a scratch copy: all four scripts exit 0 with empty stderr; all pre-existing stdout lines are byte-identical to `q3/base/` (34, 52, 27, 43 lines); the four pngs have the same md5 as the originals; each new block's numbers equal the md tables.
6. **Text claims**: Ross–Selinger (arXiv:1403.2975 abstract: "T-count 3log_2(1/ε) + O(log(log(1/ε)))" in the typical case; 3·log₂(10⁶) = 59.8); the ch11 caveat (the SWAP + X network equals U_7 on y = 1…14 and exchanges 0 and 15; the controlled version carries no relative phase); ch12 energy-QPE substitute values (‖U_s − U‖ = 0.2450 and 0.0034; modes −4.569 at 0.878, and −4.735 / −4.79 at 0.495 / 0.185; exact-U reference 0.428 / 0.238), reproduced by eigen-decomposition rather than repeated squaring; ch10's remark about ch14's CZ count (12 before routing = controlled Z counted as 2).
7. **Section 6**: one SDK each (ch10 PennyLane, ch11 Cirq, ch12 QuTiP, ch13 PennyLane); headings and reason sentences present; tables have consistent cell counts (22 tables, 0 bad rows); no remaining mention of a removed SDK or its API in any of the four chapters, and no other chapter points at a removed section-6 item; the QuTiP and PennyLane snippets are byte-identical to the pre-edit blocks after the move; code fences balanced; no unicode escapes, replacement characters or stray jamo.
8. Every md change versus `q3/base/` corresponds to an entry in `edits.json`; each `.py` is the original plus an appended block.

## Could not check

- SDK snippets (Qiskit, PennyLane, Cirq, QuTiP) cannot be run here. The PennyLane docs page confirms the `QuantumPhaseEstimation(unitary, target_wires, estimation_wires)` signature and `argmax / 2**n` read-out but does not state the wire order in words; the big-endian claim in ch10's reason sentence rests on PennyLane's standard `probs` ordering. The Cirq `ArithmeticGate` reference page was truncated on fetch.
- "수백 CZ" for the transpiled `UnitaryGate` (ch11:135) is unverified; it is plausible (generic 5-qubit lower bound 252 CX) but needs a Qiskit run.
- Stated runtimes are machine-dependent: here ch10 took 43 s (text: 38초) and ch13 took 58 s before and 60 s after the new block (text: 약 25초). The new block adds about 2.5 s to ch13, so no text change is needed on that account.
- The Beauregard count in B1 comes from my own construction with one extra expansion rule; it has not been reviewed by anyone else.
- html builds and the Claude Docs copies were not examined.

Sources: [arXiv:1403.2975](https://arxiv.org/abs/1403.2975), [Beauregard, quant-ph/0205095 (ar5iv)](https://ar5iv.labs.arxiv.org/html/quant-ph/0205095), [PennyLane QuantumPhaseEstimation](https://docs.pennylane.ai/en/stable/code/api/pennylane.QuantumPhaseEstimation.html)
