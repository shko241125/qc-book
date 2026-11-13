# R14 — Independent review of the 2026-11-12 changes (side-by-side tables, BV "0000" confirmation design, ch13/ch14/ch18 additions)

Reviewer R14 · tree reviewed: `/home/claude/qc-book` (nothing under it was edited; ch21 not read). Scripts were rerun only in a scratch copy (`.../scratchpad/r14/copy/`) with `PYTHONDONTWRITEBYTECODE=1`, so no repo file (including `appA_posthoc.json`, `appA_bv_plan.json`, `__pycache__`) was touched. Scratch code and outputs: `/tmp/claude-0/-home-claude/fd63146d-511d-5ae9-8dab-3100fb85f0db/scratchpad/r14/`.

## 0. Summary

**Counts: A = 2, B = 15, C = 14.**

The arithmetic is clean again: every number in the new text and tables matches script output, and my own statevector and density-matrix code reproduces the design check and the whole prediction table.

The problem is the experiment's central premise. The design assumes that the shots of one PUB run back to back ("한 PUB의 샷은 연달아 돈다고 본다"). IBM's own guide says the opposite: circuits of a primitive job are executed **column-wise — for each shot, the QPU iterates over all circuits of the job** — and it adds that neither the PUB order nor "one batch per job" is guaranteed. Under that documented order

- the prediction table for A1/A2 is wrong (P1 would *not* lose its 0000; P0 and P0b would differ under A2 and be flagged as "drift"),
- the fixed decision rule returns "no hypothesis matches" for a true A1 and a true A2 (I ran the author's `verdict()` on those cases),
- the 2026-11-11 argument "data qubits end in |1⟩ every shot yet show ≤ 1 %p" no longer bounds anything (≥ 4 repetition periods ≈ 1 ms ≈ 6 T1 lie between two BV shots), and `f = e/0.5 = 0.104` has no meaning,
- and carry-over from BV's own previous shot becomes physically implausible for 2026-11-11 (qubit 3 is used by no other PUB of that job; exp(−1 ms/175 µs) = 0.3 %), which by itself favours "history-independent" (B) or an in-circuit cause.

I could not test the device's real order (no SDK, no submission). The point is that the design must not depend on it. A two-PUB job `[P0, P1]` (≈ 4 s) separates "carry-over + consecutive", "carry-over + alternating", "history-independent" and "gone" with one comparison and tells you the execution order as a by-product (A1 below).

The second grade-A item is a wrong causal sentence in the new §4.1: the median model was pessimistic on GHZ/BV because of **readout error** (model 1.3 %, that day's measured qubits ≈ 0.5 %), not because of CZ error or the double count (the double count is in both the prediction and 사후 ①, so it cannot explain their difference).

Things that are fine: no code path submits without `real --yes`; nothing secret is printed; the side-by-side tables label post hoc values unmistakably; html was regenerated after the md edits.

## 1. What was run

| Step | Result |
| --- | --- |
| `appA_qpu_runs.py` in scratch copy | rc 0, stdout byte-identical to `appA/appA_stdout.txt`; `appA_posthoc.json` identical |
| `appA_bv_confirm.py` in scratch copy | rc 0, stdout byte-identical to `appA/appA_bv_stdout.txt`; `appA_bv_plan.json` identical |
| `ch14_hardware.py` in scratch copy | rc 0; `ch14_pred.json` identical; 7) block as quoted in ch14 |
| `ch13_variational.py` in scratch copy (1 m 43 s) | HVA p=1: E − E₀ = 1.13e-01, as in the ch13 table |
| `r14_sv.py` — own noiseless statevector simulator (two rz conventions) | task 2 (a)–(f), all confirmed (§2) |
| `r14_dm.py` — own noisy density-matrix replay of P0/P1/P2/P4 | 0.00060, e = 0.0517 ± 0.0036, 0.9115, f = 0.1035, e₁ = 0.00055, the whole §5 table, 1.00/0.60/0.18 %p, 10.7σ, 6.1 s / 9.2 s — all reproduced |
| `r14_hva.py` — own HVA p=1 optimisation | optimum −4.646249, E − E₀ = 0.11252; E(γ=0.2064, β=0.3927) = −4.6462 |
| `r14_extra.py`, `r14_two_pub.py` — supplementary numbers (first one uses a scratch copy of the author's `Replay`; stated where used) | HVA under each hypothesis, column-wise predictions, rule behaviour on mixed cases, alternative design |
| synthetic records → `appA_bv_confirm.py <record>` (scratch) | analysis path works; exposes B1 |
| `appA_bv_submit.py` | read only (no SDK here; not installed, as instructed) |

## 2. Independent recomputation

### 2.1 Noiseless design check (own code, `r14_sv.py`)

Qubits [2, 3, 4, 16, 114]; measured c0←2, c1←16, c2←114, c3←4; string = c3c2c1c0. Results identical for rz = diag(e^{−iθ/2}, e^{iθ/2}) and rz = diag(1, e^{iθ}) (global phase only). The ancilla tail rz(π/2)·sx·rz(π/2) equals e^{−iπ/4}·H.

| Check | Result |
| --- | --- |
| (a) P0 | `1011` with probability 1; ancilla ends with P(1) = 0.5 |
| (b) P1 (tail → x) | `1011` with probability 1; ancilla ends in \|0⟩ (P(1) = 0) |
| (c) ancilla starts \|1⟩ | P0: `0000` prob 1 (ancilla ends 0.5). **P1: `0000` prob 1, and the ancilla ends in \|1⟩ with probability 1** |
| (d) P2, perfect reset first | `1011` prob 1 for ancilla starting \|0⟩ or \|1⟩ |
| (e) data qubit starts \|1⟩ | q4 → `0011`, q16 → `1001`, q2 → `1010`, each prob 1 (same in P1). Extra: idle q114 → `1111`; ancilla **and** q4 → `1000`; ancilla X skipped → `1000` |
| (f) `carry_fixed_point` | Correct *under the script's own premise*. In P1 a \|1⟩ start reproduces itself (end P(1) = 0.9955 with noise, 0.00475 for a \|0⟩ start), so e ← f·[0.9955 e + 0.00475 (1 − e)]; with f = 0.1035 this contracts to e₁ = 0.00055 (mine and the script's). The script's `run(g, e)[1]` is linear in e and captures exactly this. Note the consequence: in P1 the 0000 events would come in runs (P(0000 at n+1 \| 0000 at n) ≈ f), which is a testable signature if per-shot data are kept (B6). |

### 2.2 Numbers in the new text vs script output

| Location | Value in text | Script output | Verdict |
| --- | --- | --- | --- |
| appA §4.1 table, ch14 §5 new table (all 8 rows, identical in both) | 0.918/0.920±0.004/0.955/0.959 … −4.433/−4.322±0.039/−4.493/−4.502/−4.377; diffs; 3.4/3.8/3.5 %p | A6 block | ✓ all cells |
| appA §4.1 text | e = 0.0517; 4.4σ → 1.4σ; +3.5 → +0.5; −1.0 → −2.3 | A4, A5, A6 | ✓ |
| appA §5 | 1.00·0.60·0.18 %p; e = 5.17 %; f = 0.104; table (5.15/0.11/0.06 %); p₀ < 1 %; 10.7σ; 200 sets; 6.1 s; "약 9초" | `appA_bv_stdout.txt` | ✓ (f = 0.10351; 9.2 s printed as "9") |
| appA §5 ② (older) vs §4.1 | success with e: 0.9115 vs 0.912 | 0.911516 | ✓ consistent (4 vs 3 decimals) |
| appA §6 / ch13 | 50 % ± 6 % | A5 | ✓ |
| ch13 new block | −4.646; E − E₀ = 0.113; −4.433 (−4.555 without readout); −4.322 ± 0.039; −4.493 (4.4σ); −4.377 (1.4σ); 0.646; 0.322; chain 2–3–16–23, CZ 6; γ, β | A2/A5/A6, `ch14_pred.json`, own `r14_hva.py` (0.1125) | ✓ |
| ch14 §5 (older sentence, same section as the new table) | "Grover는 0.779" | 0.778493 → table says 0.778 | ✗ see B11 |
| hand-derived, not printed by any script | "2 % 안팎" (§5, prepare-and-measure), "QPU 약 10초" (§9, submit docstring, CLAUDE.md) vs "약 9초" (§5) | — | C1 |

No other mismatch. Rounding ties that a reader can trip over: B11.

## 3. Grade A

### A1. The design assumes consecutive shots per PUB; IBM documents column-wise (interleaved) execution — predictions, rules and one supporting argument fall with it

- Where: `appA/appA_qpu_runs.md:190` ("한 PUB의 샷은 연달아 돈다고 본다"), `:151` (data-qubit argument), `:174` (f = 0.104, "다음 샷"), `:176-184` (table), `:186` (rules, P0b as drift check), `:170` ("P0b … 작업 안에서 상태가 변했는가"); `appA_bv_confirm.py:109-119` (`f_carry`, `carry_fixed_point`), `:122-129` (E table), `:163,176` (drift); `appA_bv_submit.py:147` comment; also the older (R11-derived) sentences `:145` ("다음 샷은 그 들뜸을 … 지운 뒤 시작한다") and `:147` (가).
- Evidence (IBM guide "Repetition rate execution", section on qubit initialization; fetched 2026-11-11 UTC): "Circuits you submit in a primitive job are batched together for execution on the QPU. These circuits are executed by iterating over the circuits for each shot requested; the execution is column-wise over a matrix of circuits and shots". Figure alt text: "The first column represents shot 0. The circuits are run in order from 0 through 3. The second column represents shot 1. …". And: "there is no guarantee on the order the circuits from PUBs are executed. For example, if you submit pubs=[pub1, pub2], circuits from pub1 might not run before those from pub2. There is also no guarantee that circuits from the same job would run as a single batch on the QPU." The executor options page describes `rep_delay` as "the delay between the end of one circuit and the start of the next within a shot loop". (https://quantum.cloud.ibm.com/docs/en/guides/repetition-rate-execution , https://quantum.cloud.ibm.com/docs/en/api/qiskit-ibm-runtime/ibm-quantum-schemas-executor-version-0-1-options-model.md)
- Not verified: what `ibm_fez` actually did on 2026-11-11 or does for the client-side Sampler (the records hold counts only, no shot order). R11's reasoning ("the next shot inherits…") made the same unverified assumption; I am correcting that too.
- Consequences, computed:
  1. **2026-11-11 data.** Job PUB order was GHZ, BV, Grover, QPE. Qubit 3 appears only in BV (QPE uses 16, 22, 23, 24; Grover 22, 23, 24; GHZ 121…143). If execution was column-wise, two BV shots are ≥ 4 repetition periods ≈ 1 ms apart; exp(−1000 µs / T1) = 0.0033 for qubit 3 (T1 175 µs), so a 50 % end-of-shot excitation leaves < 0.2 % even with no reset. A 5 % start excitation can then not be BV's own carry-over: it is history-independent (B), or caused inside the circuit, or the job was not run column-wise. Likewise the data qubits had ≈ 1 ms, so "초과가 1.0 %p 이하 → generic 10 % carry-over is ruled out" proves nothing in that case. (If shots were consecutive the argument is valid: predicted ≈ 10 %, observed total frequency ≤ 2.4 %, see B7.)
  2. **Proposed job, column-wise in PUB order P0, P1, P2, P3, P4, P0b** (author's model, carry-over from the immediately preceding circuit, f = 0.1035; `r14_extra.py`): A1 → P0 5.15, **P1 5.15**, P2 0.06, P3 0.06, P4 0.59, P0b 5.15 %; A2 → the same except **P0b 0.06 %**. The author's `verdict()` returns for both `표의 어느 가설과도 맞지 않는다 ('그대로', '사라짐', '사라짐')`, and for A2 additionally `P0 과 P0b 가 +14.5σ 다르다(작업 안의 드리프트)` with p₀ averaged to 2.6 %. B still classifies as B. So a true carry-over cause can never be recognised, and P1 = "그대로" (the B signature) is also what A1/A2 give.
  3. With column-wise execution P0 and P0b are sampled over the same seconds, so their difference is never drift; it is an order effect.
- Fix. Do not assume an order; make the first job decide it. Smallest decisive experiment — one job with exactly two PUBs `[P0, P1]`, 4,000 shots each, per-shot bitstrings saved. Predictions (own replay, `r14_two_pub.py`):

| Case | P0 | P1 |
| --- | --- | --- |
| carry-over, a PUB's shots consecutive | 5.15 % | 0.11 % |
| carry-over, column-wise alternation (P0 follows P1 and vice versa) | 0.63 % | 5.15 % |
| history-independent (B) or in-circuit cause | 5.15 % | 5.15 % |
| gone (T) | 0.06 % | 0.06 % |

  The pairs differ by ≥ 12σ at 4,000 shots. With two circuits "which PUB runs first" does not matter. QPU ≈ 4 s by the usage fit in C1. P2/P3/P4 then go into a second job whose arrangement is chosen knowing the order (or as single-PUB jobs, where the previous shot is the same circuit by construction; ≈ 2.8 s each).
- Proposed replacement for `appA_qpu_runs.md:190`:
  `*이 설계의 전제와 그 점검.* A1·A2의 예측은 "직전에 돈 회로가 무엇인가"에 달려 있다. IBM 문서("Repetition rate execution")는 한 작업의 회로들을 샷마다 한 번씩 돌아가며 실행한다고 적는다(샷 0에서 회로 0…k, 샷 1에서 다시 회로 0…k). PUB의 순서와 한 묶음 실행은 보장하지 않는다고도 적는다. 그렇다면 직전 샷은 같은 회로가 아니라 작업 안의 다른 회로다. 2026-11-11의 작업(GHZ·BV·Grover·QPE)에서 3번 큐비트는 BV에만 쓰였으므로, 이 순서였다면 BV의 두 샷 사이는 반복 지연 네 번(약 1 ms, T1의 6배)이고 BV 자신의 들뜸이 5% 남을 수 없다. 실제 순서는 기록에 없다. 그래서 첫 작업은 순서에 기대지 않는 회로 둘(P0, P1)만 넣는다. 잔류이고 같은 회로가 연달아 돌면 (P0, P1) = (5.15%, 0.11%), 잔류이고 번갈아 돌면 (0.63%, 5.15%), 직전 회로와 무관하면 (5.15%, 5.15%), 지금은 없으면 (0.06%, 0.06%)다. P2·P3·P4는 이 결과로 실행 순서를 안 뒤에 둘째 작업으로 낸다.`
  and `:151`: after "…1.00·0.60·0.18%p뿐이다." add `다만 이 비교는 같은 회로의 샷이 연달아 돌았을 때만 뜻이 있다. 회로를 번갈아 돌았다면 데이터 큐비트도 약 1 ms를 쉰 뒤라 잔류가 작은 것이 당연하다.` The `:174-186` table and rules must be regenerated after `appA_bv_confirm.py` models "the preceding circuit" explicitly (a `PREV = {...}` map per assumed order) instead of `carry_fixed_point(g)` on the circuit itself. CLAUDE.md lines 65–67 (order ①–④, "PUB 6개 × 4000샷") need the same update.

### A2. §4.1 names the wrong cause for the median model's pessimism

- Where: `appA/appA_qpu_runs.md:129`: "중앙값 모형은 비관적이었고(그날 쓰인 연결의 CZ 오류가 중앙값보다 작았고, 게이트 중 열완화를 두 번 셌다), 모형에 없는 손실(GHZ의 q4, BV의 보조 큐비트)이 그만큼을 채웠다."
- Evidence (ch14 output 2) "판독 오류 끔" column; appA A3 "판독 제외"; `r14_extra.py`):

| Circuit | 예측 | 사후 ① | ① − 예측 | readout off: 예측 / ① | readout loss: 예측 / ① |
| --- | --- | --- | --- | --- | --- |
| GHZ | 0.9181 | 0.9552 | +3.70 %p | 0.980 / 0.9797 | 6.2 / 2.45 %p |
| BV | 0.9245 | 0.9612 | +3.67 %p | 0.980 / 0.9890 | 5.5 / 2.78 %p |
| Grover | 0.7590 | 0.7785 | +1.95 %p | 0.791 / 0.7913 | 3.2 / 1.29 %p |
| QPE | 0.8809 | 0.9321 | +5.12 %p | 0.913 / 0.9392 | 3.2 / 0.71 %p |

  With readout switched off the prediction and 사후 ① coincide for GHZ and Grover. The model assumes 1.3 % readout error; the measured qubits that day averaged 0.51 % (GHZ) and 0.46 % (BV). For GHZ/BV (3–4 CZ) the CZ error rate is irrelevant; for Grover the lower CZ error was cancelled by 39 instead of 32 CZ. The gate-time double count is present in both the prediction and ①, so it cannot be a reason why the prediction is more pessimistic than ①.
- Also "예측이 가까웠던 것은 두 가지가 상쇄됐기 때문이다" is true for GHZ only (3.7 against 3.5 %p); in BV the missing loss (7.2 %p) is twice the pessimism; Grover and QPE have no cancellation (+5.5, +4.1 %p). The same sentence in `ch14_hardware.md:138` ("실행 전 예측이 가까웠던 것은 … 상쇄된 결과다") needs the same scope.
- Proposed replacement (after the script prints the readout-off values and the day's mean readout error):
  `상쇄가 일어난 곳은 GHZ 하나이고 BV에서는 절반만 상쇄됐다. 두 회로에서 중앙값 모형이 사후 ①보다 3.7%p 비관적이었던 까닭은 판독 오류다. 모형은 1.3%를 넣었고 그날 측정에 쓰인 큐비트는 평균 0.5%였다(판독 오류를 끄면 예측과 ①은 GHZ에서 0.980 대 0.9797로 같다). 그 3.7%p를 모형에 없는 손실이 GHZ에서는 거의 그대로(q4, 3.5%p), BV에서는 두 배로(보조 큐비트, 7.2%p) 채웠다. Grover와 QPE에서는 상쇄가 없었고 예측이 비관적인 채로 남았다(+5.5, +4.1%p). 게이트 중 열완화의 이중 계산은 예측과 ①에 똑같이 들어 있으므로 둘의 차와는 무관하다.`

## 4. Grade B

### B1. Decision rules: most realistic outcomes map to no hypothesis, and some labels are wrong
`appA_bv_confirm.py:149-178`, `appA_qpu_runs.md:186`. Checked with the author's `verdict()`:
- 27 (P1, P3, P4) patterns, 3 mapped. Unmapped but plausible: reset misfire on the worst-readout qubit (P1 ≈ P(1\|0) ≈ 2 %) → ('부분','사라짐','그대로'); half carry-over, half static → same; A2 with weaker carry-over on the other hub (2 %) → ('사라짐','부분','사라짐'). The text says of A2's P3 "0이 아니다만이 예측" but the rule demands "그대로" (> 2p₀/3).
- P3 skipped (no hub) → `('사라짐', None, '그대로')` → "맞지 않는다", although P1 and P4 alone identify A1/A2/B.
- p₀ = 0.6–0.9 % (10–15× the 0.06 % baseline, > 5σ above it) → "T — 0000 이 지금은 없다". It is present, only weaker.
- "그대로" has no test and no upper bound: p_k = 0.7 p₀ is 4.0σ below p₀ at these sizes and is still "그대로"; so is p_k = 3 p₀.
- P2 = "부분" is ignored silently. It is the expected result even under A/B: an explicit reset on qubit 3 can hardly do better than its readout error (0.026, the worst of the 13).
- Fix: report the ratio p_k/p₀ with a 95 % interval for every PUB and state the reading per hypothesis as an interval ("≤ 0.2", "0.2–0.8 = 혼합", "≥ 0.8 and consistent with 1"); give P1-and-P4-only readings; rename the T outcome "약해졌다(p₀ < 1%) — 판정하지 않는다".

### B2. "P2 그대로 → C(누설)" does not follow; add the free controls that separate the alternatives
`appA_bv_confirm.py:174-175`, md `:167`, `:180`. P2 unchanged is also what you get if (i) the initialization reset itself is the problem (an explicit reset is the same operation), (ii) no pulses were applied in those shots (0000 in every variant), (iii) the excitation happens inside the circuit before the first CZ (B3). Free additions: an X on the gate-free measured qubit 114 as a marker (R11's P3; a no-pulse shot then reads 0 there) — it was dropped without explanation; and a pre-registered reading of the P4 quantity the script already prints (`0000 인 샷 가운데 보조를 1 로 읽은 비율`): ≈ 0.5 if the ancilla started in \|1⟩ (X then H), ≈ 0 for no-pulse shots, ≈ 1 for \|2⟩ (read as 1). Change the verdict text to `P2 가 그대로다: 맨 앞의 reset 으로도 없어지지 않는다(초기화 reset 자체, 게이트가 걸리지 않은 샷, 회로 안에서 생기는 들뜸, 누설 가운데 하나)`.

### B3. A missing hypothesis and a sharper measurement: read the ancilla at the start of the BV circuit
All hypotheses share the unverified premise "0000 ⇔ ancilla \|1⟩ at the first CZ, and it got there before the circuit began". Not covered: D — the ancilla is excited inside the circuit (its three neighbours 4, 16, 2 all receive sx at t = 0 while it idles); pattern (그대로, 사라짐, 그대로) = verdict "B". Sharper circuit Qm: P0 with a mid-circuit `measure` of the ancilla into an extra bit before the first CZ (no feed-forward needed). Per shot, P(0000 \| start bit = 1) ≈ 1 and P(start bit = 1 \| 0000) ≈ 1 − P(0\|1) ≈ 0.97 if the premise holds; ≈ P(1\|0) ≈ 0.02 if 0000 arises later (D) or elsewhere. The md rejects "준비 후 즉시 측정" because it "보정 실험과 같은 양을 재므로 2 % 안팎" (`:149`); that applies to a stand-alone circuit, not to a measurement embedded in the BV circuit in the same job, where the shot-by-shot correlation removes the readout floor. Caveat: I could not check that the Open plan/`ibm_fez` accepts a measurement followed by gates without changing timing of the other qubits.

### B4. The HVA "independent confirmation" is a hypothesis-B calculation, and it already ranks the hypotheses
`appA_qpu_runs.md:127`, `:147` ("서로 다른 두 작업이 같은 크기를 가리킨다"), `:203`; ch13 `:175`. In the HVA job qubit 3 is **measured** every shot, and its end state differs by basis (P(1) = 0.498 in the Z circuit, 0.087 in the X circuit). Putting the same e in both bases is what B predicts. With the author's replay (scratch copy, `r14_extra.py`):

| Hypothesis for HVA | e (Z, X) | ⟨H⟩ | 실측 − 값 |
| --- | --- | --- | --- |
| B (what A5 computes) | 0.0517, 0.0517 | −4.3774 | +0.055 (1.4σ) |
| A1, previous circuit = same basis | 0.0516, 0.0090 | −4.4444 | +0.122 (3.1σ) |
| A1, bases alternate | 0.0090, 0.0516 | −4.4064 | +0.084 (2.2σ) |
| A2 (qubit 3 is measured → no carry-over) | 0, 0 | −4.4930 | +0.171 (4.4σ) |

σ is shot noise only, so this is weak evidence, but it costs nothing and should be stated: `같은 e를 두 기저에 똑같이 넣는 것은 "직전 회로와 무관하다"(B)는 가정이다. HVA에서 3번 큐비트는 매 샷 측정되므로 A2라면 e = 0이어야 하고 차는 4.4σ로 남는다. A1이라면 기저마다 e가 달라 2.2~3.1σ다. 샷 잡음만 센 값이라 결정적이지는 않지만 이미 있는 자료는 B > A1 > A2 순으로 가리킨다.` (print the three values in A5 first).

### B5. P3: choosing the hub with the lowest readout error biases the test against A2; the plan is not pinned
`appA_bv_submit.py:39-52`. A2 claims "any qubit". Picking the best-read (and hence probably best-reset) hub makes "사라짐" likelier, which the rule reads as "qubit 3 specific". Use a hub whose readout error and T1 are closest to qubit 3's, or two hubs (median and worst), and ignore neighbours' quality only knowingly (the script never looks at the neighbours' or the edges' calibration, nor at `props.faulty_gates()`). Also `real` (plan) and `real --yes` each call `pick_hub`; if the calibration updates in between, the submitted hub differs from the one the user approved. Add `--hub N` (or write the chosen map into a plan file that `--yes` must read).

### B6. `appA_bv_submit.py`: keep what the analysis will need
1. `:152` saves counts only. Save per-shot data too (`pub.data.c.get_bitstrings()`, 4,000 short strings per PUB): it is the only way to see execution-order effects, runs of 0000 in P1 (§2.1 f), bursts in time, and the Qm correlation (B3). Unverified: that the array order equals acquisition order.
2. Options actually used are not recorded (`rep_delay_default` only). Set them explicitly and dump them: `sampler.options.execution.init_qubits = True`, `rep_delay = backend.default_rep_delay`, and store `dataclasses.asdict(sampler.options)` if available. R11 already noted that the 2026-11-11 records cannot say which options ran.
3. `:150` `job.result()` — if it raises or is interrupted the record keeps the job id but there is no path back: `fetch_ibm_job.py` stores counts in another shape and `appA_bv_confirm.py` cannot read it. Add `real --resume <job id> <record>`.
4. Local mode prints "(1.0000 이어야 한다)" but does not assert; `env/test_all.sh --sdk` would report 통과 for a wrong circuit. Add `assert good == 1.0`.
5. The P3 path (relabelling) is never simulated. `appA_bv_confirm.py:87-88` asserts only that gate *names* are unchanged, which is true by construction, and `relabel()` is otherwise dead code (the submit script remaps inside `build`). Simulate the relabelled list noiselessly instead: `assert abs(p_data(R2({…}, CAL, noisy=False).run(relabel(P0, _mp)).dist(…), GOODS) − 1) < 1e-9`.
6. `has_reset = "reset" in tgt.operation_names` → `tgt.instruction_supported("reset", (anc,))`; record `tgt["reset"][(anc,)].duration` (the reset delays `cz[4, 3]`, so qubit 4 idles in superposition; this changes single-bit errors, not 0000).

### B7. "초과 = 그 큐비트의 초기 들뜸의 상한" is not an upper bound
`appA_bv_confirm.py:102`, md `:151`. The excess over the replay is an estimate that can be too low if the replay overstates other routes to the same string (e.g. calibrated P(0\|1) too high that day). The bound is the observed frequency itself: 2.03 %, 2.38 %, 1.00 % for qubits 4, 16, 2. The conclusion (≪ 10 %) survives. Proposed: `…인데, 그 오답의 실측 빈도는 2.03·2.38·1.00%이고(이것이 초기 들뜸의 상한이다) 재현보다 많은 몫은 1.00·0.60·0.18%p다.`

### B8. "더 많이 아는 모형이 평균으로는 더 가깝지 않았다 (3.4 → 3.8 → 3.5 %p)" — true as arithmetic, fragile as a reading
`appA_qpu_runs.md:129`, `ch14_hardware.md:138`.
- One circuit decides it. Without BV: 3.28 → 2.66 → 2.20 %p (monotonically better). With BV's ③ value in ①: 2.56 %p. Shot noise on a mean of four |differences| is ≈ 0.25 %p, so 3.4 vs 3.5 is no difference.
- The signed means are +1.6, −2.0, −3.3 %p: the two informed models are optimistic in six of eight differences (only Grover is positive in each), which is what models that omit loss channels should be. So "②는 QPE에서 더 멀어진다" is not evidence against removing the double count; ①'s agreement with QPE (−1.0) may be the same kind of cancellation as GHZ's. Two statements are being mixed: that r already contains gate-time relaxation (true by the definition of RB, no data needed) and that removing it improves the fit (only Grover).
- A pre-run value is missing from the table: fake (0.9685/0.9782/0.8742/0.9530, run 15:37, one minute before the device) — mean |실측 − fake| = 5.7 %p, the worst. "예측" in the side-by-side table is therefore the better of two pre-run values. Add a fake column or a sentence.
- Proposed: `성공 확률 넷의 평균 차는 예측 3.4, 사후 ① 3.8, 사후 ② 3.5%p로 사실상 같다(샷 잡음만으로 ±0.25%p). 이 순서는 BV 하나가 정한다. BV를 빼면 3.3, 2.7, 2.2%p로 정보가 많을수록 가깝다. 부호를 보면 ①과 ②는 여덟 값 중 여섯에서 실측보다 낙관적이다(Grover만 반대다). 손실 기제를 빠뜨린 모형이 보일 모습이다. 그러므로 ②가 QPE에서 더 멀어진 것은 이중 계산을 빼는 것이 틀렸다는 근거가 아니다. r에 게이트 중 열완화가 들어 있다는 것은 RB의 정의에서 나오고, 실측이 보태는 근거는 Grover 하나다. 실행 전에 낸 값은 fake도 있었고 그 평균 차는 5.7%p였다.`

### B9. Sentences elsewhere that now contradict the new caution or the new naming rule
- `ch14_hardware.md:123` (2): "실측과의 나머지 차는 이 장의 모형이 게이트 중 열완화를 두 번 세는 데서 온다" — stated as fact; the new block 15 lines below says the evidence is one circuit. → `…나머지 차는 게이트 중 열완화의 이중 계산을 빼면 0.5%p로 줄어든다(사후 ②). 같은 수정이 다른 세 회로에서는 차를 키운다(아래 병기 표).`
- `appA_qpu_runs.md:106`: "Grover는 … 이것은 큐비트가 아니라 모형 쪽의 문제다(3.4절)" → `…모형 쪽에서 설명될 수 있다(3.4절, 4.1절)`.
- `ch14_hardware.md:189`: "이 수정은 Grover에서만 실측에 가까워지고 QPE에서는 멀어진다" — GHZ (−3.5 → −3.9) and BV (−7.2 → −7.5) also move away. → `Grover에서만 가까워지고 나머지 셋(GHZ·BV·QPE)에서는 멀어진다`.
- `appA_qpu_runs.md:90`: "같은 회로를 세 가지 모형으로 예측해 실측과 나란히 놓으면" — the third is 사후 ①; CLAUDE.md now says "사후 값을 예측이라 부르지 않는다". → `같은 회로에 대한 세 가지 모형의 값(실행 전의 둘과 실행 뒤의 하나)을 실측과 나란히 놓으면`.
- `ch13_variational.md:263`: "그래서 14장의 예측·제출 목록에 들어간다(실기 실행은 아직 못 했다)" contradicts the new §1 sentence and §5 block. → `(HVA p = 1은 2026-11-11에 돌았다 — 5절)`.
- (가)(나)(다) vs A1/A2/B/C: the mapping is stated (`:151`, `:157-158`) and consistent. OK.

### B10. `appA_stdout.txt` A6 last line contradicts its own table
`appA_qpu_runs.py:332`: "뒤로 갈수록 실측에 가깝지만 그만큼 실측을 보고 만든 값이다." The line above prints 3.4 → 3.8 → 3.5 and the md says "더 가깝지 않았다". Also "② 는 모형의 이중 계산을 뺀 값". → `읽는 법: ① 은 '어느 큐비트·어느 날'을 넣은 값, ② 는 ① 에서 게이트 중 열완화를 뺀 값, ③ 은 맞춘 매개변수 하나(e)를 쓴 값이다. 뒤로 갈수록 실측을 보고 만든 값이고, 평균 차가 줄지는 않았다(③ 을 쓴 BV 와 HVA 만 줄었다).` Pin the reading: `assert _ad('pred') < _ad('post_no_gate_relax') < _ad('post_calib')`.

### B11. Rounding: "0.779" vs "0.778" in the same ch14 section, and two visible ties
`ch14_hardware.md:123` "Grover는 0.779" vs the new table `:131` "0.778" (value 0.778493; 0.779 is a double rounding of 0.7785). The appendix prints 0.7785 (§3.4, §4) and 0.778 (§4.1), and 0.9635 (§3.4) and 0.963 (§4.1; value 0.963455) — a reader rounding the four-decimal figures gets 0.779 and 0.964. Fix: use four decimals for the three 사후 columns in both side-by-side tables (as the §4 table does), and change `:123` to 0.7785.

### B12. ch13: "그 손실의 큰 몫이 큐비트 하나(물리 3번)에 걸려 있다는 것이 부록 A 6절의 분석이다" overstates §6
§6 says the readout share inside the *model* is 0.068 of 0.089 and that with e the gap shrinks, "다만 HVA는 에너지 하나만 주므로 다른 원인과 가를 수는 없다". Of the 0.324 lost, 0.068 is modelled readout on qubit 3 and 0.116 is the fitted, unconfirmed e. → `그 손실 0.324 가운데 0.068은 재현 모형 안에서 3번 큐비트의 판독 오류 몫이고, BV에서 맞춘 초기 들뜸 e를 같은 큐비트에 넣으면 0.116이 더 설명된다(사후 ③). 뒤의 것은 확인되지 않은 가설이다(부록 A 5·6절).` (0.116 = 4.4930 − 4.3774; have A5 print it.)

### B13. ch18 paragraph: "같게 놓는다" is not what ch00/ch19 define, and it covers only one of the four uses of 6′
`ch18_search_tools.md:146`.
- ch00 TTS₉₉ and the ch19 pre-registration compare cost-to-solution in a common **unit** under a common cap ("단위는 목적함수 평가 1회(회로 실행 1샷 = 담금질 1걸음 = 무작위 비트열 1개) … 예산 10⁷ 평가"). The counts are not set equal; TTS is the thing compared. "평가" alone is also the word ch19 exercise 3 uses for the *wrong* unit (one expectation estimate = 1,000 shots).
- ch19 uses 6′ for kernels ("같은 특성의 고전 커널과의 성능 차", `:50`), VQE (`:290`, exact diagonalisation) and Max-Cut (`:296`, "보장값 비교"). None of these is a count of objective evaluations.
- "이전 이름 '6′ 실측 비교'를 바꿨고" is revision history in book text; "21장이 이 단계를 실제로 돌린다" sits one sentence after "실측" and can be read as a device run (I did not read ch21).
- Proposed: `6′의 "같은 예산"은 양쪽의 비용을 같은 단위로 세고 같은 상한 안에서 견준다는 뜻이다. 조합 최적화에서는 단위가 목적함수를 한 번 계산하는 일이다(회로 실행 1샷 = 담금질 1걸음 = 무작위 비트열 1개). 변분 알고리즘이 훈련에 쓴 샷도 같은 단위로 세어 TTS₉₉에 한 번 더한다(0장, 19장 3절). 커널에서는 같은 특성과 같은 훈련 데이터를 준 고전 커널이 상대다. 어느 경우든 장치에서 잰 값끼리의 비교라는 뜻은 아니고, 이 책에서 QPU 값에만 쓰는 "실측"(0장)과 무관하다. 21장이 이 단계를 numpy 계산으로 수행한다.`
- `grep "6′"` over ch17–ch20: 12 hits, all "같은 예산 비교" or bare "6′"; no "실측 비교" left. CLAUDE.md notes the Claude Docs copies of ch17–19 still carry the old name.

### B14. Step (1), the free checks, has no pre-registered reading, and one free check is missing
`appA_qpu_runs.md:188`. (i) `init_error`: under A-type carry-over with consecutive shots e ≈ 0.5 × init_error(q3), so init_error(q3) ≈ 0.10 and ≲ 0.024 on 4/16/2; a value ≈ 0.01–0.02 on qubit 3 rules that reading out before any QPU time is spent. Write these numbers down now. (ii) Under B, a static 5 % also contradicts the calibration's own P(1\|0 prep) = 0.0195 ± 0.0043 (which contains the initialization error) unless it varies in time; `calib_now` records today's value, so state in advance what "P0 ≈ 5 % with P(1\|0 prep) ≈ 2 %" would mean (context dependence — B3/D). (iii) Missing: the per-shot bitstrings of the finished job `dbq0qmoqtr0c738f53u0` can still be fetched with no QPU use (`service.job(id).result()[1].data.c.get_bitstrings()`); whether the 206 `0000` shots are spread evenly or bunched in time separates a steady 5 % from an episode, and the same holds for GHZ's q4. `fetch_ibm_job.py` saves counts only.

### B15. The self-check is circular and does not cover "every hypothesis"
`appA_bv_confirm.py:189-200`, md `:186` ("가설마다 … 규칙은 매번 그 가설을 골랐다"). Samples are drawn from the same PRED table that defines the mapping, at 5.15 % vs ≤ 0.11 %; it can only fail through shot noise, which at 10σ separation it cannot. It says nothing about model misspecification (A1, B1). C is not tested (PRED\[P1\]\[C\] is None) and "C" never appears as a verdict prefix. → `자체 점검: A1·A2·B·T의 예측값 그대로에서 뽑은 합성 표본 200묶음에서 규칙은 매번 그 가설을 골랐다. 이것은 샷 잡음이 판정을 바꾸지 않는다는 확인일 뿐이고, 예측 모형이 틀렸을 때의 동작은 시험하지 않는다(C는 대상이 아니다).`

## 5. Grade C

- **C1 QPU time.** "약 9초" (§5) vs "약 10초" (§9, submit docstring, CLAUDE.md). The "×1.5" scaling is not the right model: the three records fit usage ≈ 1.7 s + 267 µs × shots (Bell 1,024 → 2 s; Sampler 16,000 → 6 s; the Estimator's 3 s then corresponds to ≈ 4,800 shots). That gives 8.1 s for 24,000 shots, 3.9 s for the two-PUB job, ≈ 17 s for six single-PUB jobs. Integer-second records make all of this ± 0.5 s. `rep_delay` 250 µs is consistent with the 267 µs slope; the script records `default_rep_delay` at submission.
- **C2 runtimes.** "약 3초"/"약 5초": here 5.6 s and 6.6 s (the second imports the first).
- **C3** `appA_bv_stdout.txt` says "C 의 P1·P3·P4 는 모형이 없어 적지 않았다" while the table prints C's P3 (0.06 %) and P4 (5.15 %); the md says only P1 has no model.
- **C4** `appA_bv_confirm.py` imports `appA_qpu_runs`, which recomputes everything and rewrites `appA_posthoc.json` at import. Guard the write (and ideally the prints) with `if __name__ == "__main__":`. `env/test_all.sh` runs `appA_bv_confirm.py`, (`--sdk`: `appA_bv_submit.py`), `appA_qpu_runs.py` in that order, sequentially, each from inside `appA/`: no race, the plan file exists before the submit script needs it, the JSON is rewritten twice with identical content. Outputs go to `runs/regress/`, so nothing checks that `appA_stdout.txt`/`appA_bv_stdout.txt` are current.
- **C5** `test_all.sh` SDK detection: `from qiskit import …` (line 15) matches, so without `--sdk` the submit script is skipped and with `--sdk` it runs with no arguments → local Aer mode → no device contact. It tests `import qiskit` only; with qiskit but without qiskit-aer the run is reported as 실패 instead of 건너뜀.
- **C6** `appA_bv_submit.py:69,160` write raw exception text into the committed record; `get_service` truncates to 200 characters for the same reason. Truncate here too.
- **C7** "P4 … 검토 뒤에 더한 회로다" (`:169`) reads as if the review proposed it; R11's list had a different P4 (measure–measure–X). Say `독립 검토 뒤에 저자가 더한 회로다` and list what was dropped from R11's proposal and why (marker on 114, double measurement, rep_delay scan).
- **C8** By the appendix's own new rule ("다음 실행에서는 두 모형의 예측을 둘 다 실행 전에 적어 두고 견준다"), the confirmation run should pre-register the success probabilities of its PUBs under both models; the plan holds P(0000) only.
- **C9** `appA_bv_confirm.py:135-136`: two comments say the same thing in different words; `:64` `[m for m in A.idle]`; `E["P3"]["C"] = 0.0` is an assumption presented as a prediction.
- **C10** ch13 table lists 사후 ① and ③ and skips ② (−4.502) without a word; ch13 `:13` gives device and date but not the job id (it is in §5); `:173` "resilience 0" rests on the script source, not on the record (the records hold no options — R11 §7).
- **C11** Glossary (§8) lacks the terms the new §5 uses: 사후 확인값, 반복 지연(rep\_delay), reset, 허브 큐비트(차수 3), 드리프트, 합성 표본. A reader new to hardware meets them undefined.
- **C12** md escaping: mechanical check of the four md files — all tables have consistent column counts, `|` inside cells is `\|`, no paragraph starts with `|`, `p\_k`, `rep\_delay`, `\~` are escaped; `<` in "p₀ < 1%" renders as `&lt;` in the html. OK. The four html files are newer than their md and contain the new phrases.
- **C13** `appA_qpu_runs.md:149`: "2% 안팎이 나올 가능성이 높다" is a hand estimate; say `보정값 P(1\|0 준비) = 0.0195와 같은 양이므로 그 근처가 나올 것이다`.
- **C14** `appA_bv_submit.py` exits if `cz[4, 3]` is not in the target "방향이 바뀌었을 수 있다"; CZ is symmetric, so flipping the pair is equivalent and would be the friendlier behaviour. (On 2026-11-11 the ISA contained both `cz[22, 23]` and `cz[23, 22]`, so the target lists both directions.)

## 6. Answers to the review questions, briefly

1. **Numbers** — all match; one stale rounding (B11); two hand-derived estimates (C1, C13).
2. **Recomputation** — (a)–(e) confirmed; (f) the ancilla ends in \|1⟩ and `carry_fixed_point` handles it correctly under its own premise.
3. **Design logic** — hypotheses are well defined only under an execution order that IBM's documentation contradicts (A1). Beyond that: rules leave realistic outcomes unmapped (B1), P2 cannot identify leakage (B2), an in-circuit cause and the direct start-of-circuit measurement are missing (B3), existing HVA data already discriminate (B4), P3's hub choice is biased (B5), the data-qubit "upper bound" is not one (B7), the self-check is circular (B15). P1's X-for-H swap and P4's extra readout are not confounds for 0000 by themselves. QPU estimate credible (C1).
4. **Submit script** — API usage agrees with the documentation I could reach: `qubit_property(q, name)` → (value, datetime); `faulty_qubits()`, `readout_error(q)` exist; `backend.default_rep_delay` is the documented accessor; `reset` is listed among Heron supported instructions (shown for `ibm_marrakesh`); Sampler twirling and DD default to off; `shots=` keyword on `run`. Bit handling `s[-4:]`/`s[0]` is right. Safety: `sampler.run` is reachable only with both `real` and `--yes`; local mode imports Aer only. Refuses to overwrite (name is time-stamped to the second anyway). Could not verify: that the 0.52 top-level `Sampler` is the executor one and accepts a 156-qubit hand-built circuit with `reset` unchanged; `job.metrics()` layout; BitArray shot order; Open-plan acceptance of `reset`/mid-circuit measurement. The 0.52/0.53 release notes were not reachable in this session.
5. **Honesty of the side-by-side presentation** — labelling is unmistakable in all three places and in the JSON note. Weak spots: the reading of the mean (B8), the wrong cause (A2), contradictory older sentences (B9, B10), ch13 overstatement (B12). No "실측" on a non-QPU value; one post hoc value still called 예측 (`appA:90`, B9).
6. **ch18** — B13.
7. **Conventions / code** — C3–C12, B6.

## 7. Verified OK

- Stdout of both appendix scripts, `appA_posthoc.json`, `appA_bv_plan.json`, `ch14_pred.json`: reproducible byte for byte.
- Every cell of the two side-by-side tables (appA §4.1, ch14 §5) and of the ch13 table.
- §5 prediction table, f, e₁, 10.7σ, 1.00/0.60/0.18 %p, 6.1 s: reproduced with independent code.
- Noiseless outputs of P0, P1, P2, P4 and the three data-qubit strings; rz convention immaterial.
- HVA p = 1: (γ, β) = (0.2064, 0.3927) is the optimum, E = −4.6462, E − E₀ = 0.1125 → "0.113" consistent with the ch13 table (1.13×10⁻¹).
- ch13 statement that the p = 4 training result (10⁻¹⁰) is not observable on this device: fair (shot error alone is 0.039; half the p = 1 descent was lost with 6 CZ).
- (가)(나)(다) ↔ A1/A2/B/C mapping is explicit and consistent.
- No submission path without `real --yes`; no key/CRN in code, arguments or output; record name collision is refused.
- `env/test_all.sh`: correct treatment of both new scripts, no ordering problem.
- md escaping and table structure; html regenerated.
- Nothing has been submitted: no `*_bv_confirm_record.json` in `runs/`.

## 8. Could not verify

- The device's actual shot order (2026-11-11 or now), and whether the client-side Sampler changes it. A1 rests on IBM's written description.
- Anything requiring the SDK or the account (see 4 above).
- That calibration experiments for `prob_meas1_prep0` use the same initialization and repetition delay as user jobs.
- ch21 (out of scope by instruction), so the ch18 sentence "21장이 이 단계를 실제로 돌린다" was checked for wording only.
- My proposed code changes for `appA_bv_submit.py` are untested.

## Sources

- https://quantum.cloud.ibm.com/docs/en/guides/repetition-rate-execution (also the `.md` and raw `.ipynb` renderings of the same guide: column-wise execution, no order/batch guarantee, implicit reset, `init_qubits`, `backend.default_rep_delay`)
- https://quantum.cloud.ibm.com/docs/en/api/qiskit-ibm-runtime/ibm-quantum-schemas-executor-version-0-1-options-model.md (`init_qubits`, `rep_delay` "within a shot loop")
- https://quantum.cloud.ibm.com/docs/guides/qpu-information (`init_error`, `prob_meas1_prep0`, readout error; supported instructions incl. `reset`, `measure_2`, `reset_2`, `if_else`)
- https://quantum.cloud.ibm.com/docs/en/api/qiskit-ibm-runtime/models-backend-properties (`qubit_property`, `faulty_qubits`, `readout_error`)
- https://quantum.cloud.ibm.com/docs/en/guides/sampler-options (execution options, twirling/DD defaults, shots precedence, client-side import path)
- https://quantum.cloud.ibm.com/docs/en/guides/migrate-from-sampler-to-executor.md , https://quantum.cloud.ibm.com/docs/guides/get-started-with-executor (Executor: ISA requirement, no built-in suppression)
