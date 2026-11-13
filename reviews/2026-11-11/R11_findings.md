# R11 — Independent review of Appendix A (first QPU runs) and the ch14 changes

Reviewer R11 · 2026-11-11 · tree reviewed: `/home/claude/qc-book` (files under review byte-identical at start and end of the review; other chapters were being edited concurrently by someone else, not by me). Nothing under `/home/claude/qc-book` was edited. All scratch code and outputs are in `/tmp/claude-0/-home-claude/fd63146d-511d-5ae9-8dab-3100fb85f0db/scratchpad/mid/R11/`.

## 0. Verdict

The arithmetic is clean: every number in `appA_stdout.txt`, in `appA_qpu_runs.md`, in the ch14 section-5 table and in report §8.1 was reproduced by independent code from the frozen JSON records, including the calibration replay (to 4 decimals), the 345/12 fault enumeration, e = 0.0517 and the HVA replay. The bit-order handling is correct everywhere I checked (GHZ reversed exactly once; BV/Grover already book order; QPE is "x in binary").

The interpretation is where the problems are. Seven claims are wrong or rest on a wrong premise (grade A), and the most consequential ones all concern the BV/"qubit 3" story and the meaning of the replay:

1. The calibration values were **not** "measured 37 min 20 s before the job". `last_update_date` is the last calibration event; the property values come from (daily) benchmarking and carry their own timestamps, which `fetch_ibm_job.py` throws away.
2. `e / P(1|0 prep) = 2.65` does **not** compare like with like, and `P(1|0)` is **not** an upper bound on the initial excitation in this circuit. IBM publishes a separate `init_error` property for exactly the situation in BV (qubit left excited by the previous shot). The ancilla ends every BV shot in |−⟩, unmeasured.
3. The Grover conclusion ("RB → depolarizing is pessimistic for deep circuits", three candidate causes) misses the most direct cause: the replay model (and ch14's) adds gate-time T1/T2 **on top of** RB error rates that already contain it. Removing that double count moves the Grover replay from 0.7785 to 0.8090 (measured 0.8137; +0.5 %p, 0.8σ).
4. Two plain factual errors: "BV의 손실 대부분" (0000 is 46 % of BV's loss) and "실측은 BV를 빼면 둘 사이에 있었다" (HVA is not between either).
5. "유일한 단일 결함 기제" is not true at circuit level (a single faulty `cz[4, 3]` also gives 0000 — that is where the replay's own 0.0006 comes from), and "12 of 345" counts list positions, not distinct faults.

One positive result the appendix does not use: putting the BV-fitted e = 0.0517 on qubit 3 into the HVA replay (one line with the existing `init_exc` hook) moves measured−replay from +4.4σ to +1.4σ. And a 4-parameter "initial |1⟩ population per qubit" model fits the whole BV histogram (χ² 17417 → 0.3 over 8 bins) and the GHZ histogram (725 → 26 over 16 bins). The "qubit starts in |1⟩" reading is therefore better supported than the appendix says — but its *cause* is more likely circuit-dependent carry-over than "qubit 3 got worse within 37 minutes".

**Counts: A = 7, B = 16, C = 18.**

## 1. What was run

| Step | Result |
| --- | --- |
| `cd appA && python3 appA_qpu_runs.py` | rc 0, output byte-identical to `appA_stdout.txt` |
| `ch14_hardware.py` (run in a scratch copy so `ch14_results.png` in the repo is not rewritten) | rc 0; blocks 1–6 identical to the pre-merge version's output; block 7 as quoted in ch14 |
| syntax check (`ast.parse`) of `ch14_submit.py`, `ch14_hardware.py`, `fetch_ibm_job.py`, `run_ibm.py`, `appA_qpu_runs.py` | all OK (not `py_compile`, to avoid writing `__pycache__` into the repo) |
| `r11_counts.py` — success, SE, TVD, σ, Bell, QPU seconds, wrong bitstrings in book order | all match (§2) |
| `r11_replay.py` — own kron-based density-matrix replay (GHZ, BV, Grover, QPE) | 0.9552 / 0.9612 / 0.7785 / 0.9321, readout-off 0.9797 / 0.9890 / 0.7913 / 0.9392, GHZ q4-only 0.0055, BV P(0000) 0.00060 — all identical to the appendix |
| `r11_more.py` — model-choice sensitivity, calibration sampling error, own statevector fault enumeration (1-qubit and circuit-level), e fit | 345/12 reproduced; 84 circuit-level faults / 6 hits; e = 0.0517 ± 0.0036 |
| `r11_hva.py`, `r11_hva_isa.py` — HVA replay (author-style reconstruction and the *actual* ISA gate list with fitted basis blocks), init error on q3, estimator std | −4.4930 / −4.5816 / −4.5610 reproduced; actual-ISA replay −4.4967 (Δ −0.004) |
| `r11_carryover.py` — per-qubit initial-excitation fits to BV and GHZ | see A2 |
| `appA_qpu_runs_patched.py` (scratch) — the proposed additions applied to a copy of the author's script | runs; prints the new numbers quoted below; diff in `r11_appA_patch.diff` |

No SDK was available: `ch14_submit.py`, `fetch_ibm_job.py`, `run_ibm.py` were reviewed by reading only. The code replacements I propose for those three are untested.

## 2. Independent recomputation — result

Bit order of the records: `bell2.json` has both orders (book = reversed Qiskit, verified). The two `*_ch14_record.json` files have **no** `bit_order` key: GHZ is `measure_all` → Qiskit order (reversed for book order), BV/Grover reversed classical mapping → already book order, QPE direct mapping → the string is x in binary.

| Quantity | Mine | appA stdout / md | ch14 | report §8.1 |
| --- | --- | --- | --- | --- |
| Bell TVD, leak, shot-noise mean, 00:11, σ | 0.02832, 0.02637, 0.01247 (exact binomial 0.01246), 0.4844, −0.98σ | 0.0283, 0.0264, 0.0125, 0.484, −1.0σ ✓ | — | 0.0283 ✓ |
| GHZ success ± SE, TVD | 0.92050 ± 0.00428, 0.07950 | ✓ | 0.920 ± 0.004 ✓, TVD "0.079" (see C1) | ✓ |
| BV | 0.88875 ± 0.00497, 0.11125 | ✓ | ✓ | ✓ |
| Grover | 0.81375 ± 0.00616, 0.13156 | ✓ | ✓ | ✓ |
| QPE | 0.92225 ± 0.00423, 0.07775 | ✓ | ✓ | ✓ |
| measured − ch14 prediction (%p) | +0.25, −3.63, +5.47, +4.13 | ✓ | ✓ | ✓ |
| fake success ± SE | 0.9685 ± 0.0028, 0.9782 ± 0.0023, 0.8742 ± 0.0052, 0.9530 ± 0.0033 | ✓ | 0.969/0.978/0.874/0.953 ✓ | ✓ |
| GHZ wrong, book order | 11110:74, 00001:47, 01111:31 | ✓ | — | ✓ ("q4") |
| BV wrong | 0000:206, 1001:95, 0011:81 | ✓ | — | ✓ |
| Grover wrong | 111:153, 001:130, 100:99 | ✓ | — | — |
| QPE wrong | 000:90, 011:55, 101:54 | ✓ | — | ✓ |
| HVA | −4.3224 ± 0.03895; −pred(−4.4329) = +0.1105 = +2.84σ | "+0.111 (2.8σ)" | "**+0.110**(2.8σ)" (B11) | "0.11 (2.8σ)" |
| QPU seconds | 2 + 6 + 3 = 11 = 1.83 % of 600 | ✓ | 6 s, 3 s ✓ | ✓ |
| e fit, success with e, residual | 0.0517 ± 0.0036, 0.9115, −2.28 %p | 0.0517, 0.9115, −2.3 %p ✓ (no ± given) | — | "2.65배" ✓ |
| HVA descent fractions | 100 / 90 / 76 / 50 % (measured 49.9 ± 6.0 %) | ✓ (no ± given, C4) | — | — |

Only one numeric mismatch between documents (B11) and three rounding-tie cosmetics (C1). Everything else agrees.

Observation worth knowing (C2): in all four circuits the "TVD" column is *exactly* (ideal success − measured success) — 0.0795 = 1 − 0.9205, 0.1112 = 1 − 0.8888, 0.1316 = 0.9453 − 0.8137, 0.0777 = 1 − 0.9223 — because every wrong outcome is above its ideal weight. The column carries no information beyond the success column in this data.

## 3. Grade A findings

### A1. "보정값은 작업 37분 20초 전에 잰 값" — wrong premise
- `appA/appA_qpu_runs.md:104`, `:122` ("그 37분 사이에"), `:161` (glossary), `:193`; `ch14/ch14_hardware.md:178` ("보정값은 몇십 분 전의 평균이고"); `appA/appA_qpu_runs.py:170`, `:189` (labels "보정 시각"); `env/fetch_ibm_job.py:39` (takes `[0]`, discards the timestamp).
- Evidence (IBM docs, "Monitoring, calibrations, and benchmarking"): each property in `backend.properties()` is "paired with the timestamp of when it was measured"; `last_update_date` "records when the underlying device calibrations were last updated" and "is updated whenever a calibration occurs, independently of the per-property benchmark timestamps". T1/T2, readout fidelity, 1q RB and CZ error come from **daily** benchmarking ("Daily benchmarking gives a comprehensive view…"); only short monitoring jobs run "approximately once an hour". `BackendProperties.qubit_property` returns "a tuple of the value and the time it was measured" — the script keeps only the value. So 37 min 20 s is the age of the last *calibration event*; the numbers used in the replay may be many hours old. The record cannot say.
- Consequence: the dichotomy in the verdict ("did qubit 3's initialization get worse in those 37 minutes?") is built on this.
- Replace `appA_qpu_runs.md:104`:
  - Current: `3. **당일 보정값도 "이 값"을 주지는 못한다.** 재현과 실측의 차는 GHZ −3.5, BV −7.2, Grover +3.5, QPE −1.0%p다. 보정값은 작업 37분 20초 전(보정 15:00:59, 작업 15:38:19 KST)에 잰 값이라 한 큐비트의 순간 상태는 담지 못한다. 차가 큰 두 회로는 모두 특정 큐비트 하나로 설명된다(GHZ의 q4, BV의 보조 큐비트).`
  - Proposed: `3. **당일 보정값도 "이 값"을 주지는 못한다.** 재현과 실측의 차는 GHZ −3.5, BV −7.2, Grover +3.5, QPE −1.0%p다. 재현에 쓴 값은 작업 시각에 조회되는 보정 데이터다. 기록한 것은 그 데이터의 마지막 보정 갱신 시각(`last_update_date`)이 작업 37분 20초 전(15:00:59, 작업 15:38:19 KST)이라는 것뿐이다. T1·T2·판독 오류·게이트 오류는 항목마다 측정 시각이 따로 있고(IBM은 이 값들을 하루 단위의 벤치마킹에서 잰다), 그 시각은 받아 두지 않았다. 값이 몇 시간 묵었을 수 있고, 어느 쪽이든 한 큐비트의 순간 상태는 담지 못한다. 실측이 재현보다 나쁜 두 회로(GHZ, BV)는 손실이 특정 큐비트 하나에 몰려 있다(GHZ의 q4, BV의 보조 큐비트). Grover는 반대로 실측이 3.5%p 좋고, 이것은 큐비트가 아니라 모형 쪽의 문제다(3.4절).`
  (The last two sentences also fix B1.)
- Replace `appA_qpu_runs.md:161`: current `| 보정 데이터(calibration) | IBM이 주기적으로 재서 공개하는 큐비트·게이트별 오류율과 T1·T2 | 작업 37분 전 값을 받아 재현에 썼다 |` → `| 보정 데이터(calibration) | IBM이 주기적으로 재서 공개하는 큐비트·게이트별 오류율과 T1·T2. 항목마다 측정 시각이 따로 있다 | 작업 시각에 조회되는 값을 재현에 썼다. 마지막 보정 갱신은 작업 37분 전이고, 항목별 측정 시각은 받지 않았다 |`
- Replace `appA_qpu_runs.md:193`: `(보정 시각 2026-11-11 15:00:59 KST)` → `(마지막 보정 갱신 `last_update_date` 2026-11-11 15:00:59 KST. 항목별 측정 시각은 기록하지 않았다)`
- Replace in `ch14_hardware.md:178`: `보정값은 몇십 분 전의 평균이고, 한 큐비트의 순간 상태(초기화 실패, TLS 결함)는 담지 못한다.` → `보정값은 하루 단위의 벤치마킹에서 잰 값이고(이 실행에서는 항목별 측정 시각을 받아 두지 않았다), 한 큐비트의 순간 상태나 회로에 따라 달라지는 초기화 오류는 담지 못한다.`
- Script labels: `appA_qpu_runs.py:170` `보정 시각` → `마지막 보정 갱신(last_update_date)`; `:189` `보정 시각과 작업 생성 시각의 차` → `마지막 보정 갱신과 작업 생성 시각의 차`. Fetch fix: see B13.

### A2. "e는 보정값 상한의 2.65배" — not like with like; `P(1|0)` is not an upper bound here
- `appA/appA_qpu_runs.md:120`, `:122`; `CLAUDE.md:57`; `docs/저장소 분석·검토 보고서.md:178`; `ch14/ch14_hardware.md:123` (indirectly).
- What `prob_meas1_prep0` is (IBM docs, "View backend details"): "the probability of measuring a qubit in the |1⟩ state when it was intended to be prepared in the |0⟩ state … it reflects SPAM errors, with measurement errors being the predominant contributor". So it is a SPAM sum, as the appendix says.
- What the appendix misses: the same page defines a **separate** per-qubit property `init_error`: "The initialization error quantifies the percentage of the population in the |1⟩ state after the default rep-delay and qubit initialization procedure when the prior experiment prepared the qubit in the |1⟩ state. Note that this is subtly different than the quantity prob_meas1_prep0…" (accessed as `backend.properties().qubit_property(q, 'init_error')`). IBM's repetition-rate guide adds: "imperfections in the reset process can introduce state-preparation errors".
- Why it matters for this circuit (from the recorded ISA gate list): the ancilla's last gates are `cz[2,3]`, `rz`, `sx`, `rz` = H on |1⟩ → it ends **every shot in |−⟩ (50 % |1⟩) and is never measured**. The next shot inherits whatever the rep-delay + initialization leaves. In the calibration's "prepare 0" experiments the qubit is not left excited. So `P(1|0)` bounds the initial excitation *for the calibration's history*, not for BV. The like-for-like comparison is e ≈ 0.5 × `init_error`(q3) (+ thermal floor) — a number that can be fetched for the job's timestamp **without any QPU time**, and was not.
- Additional points: (i) all readout values in the JSON are multiples of 1/1024 (20/1024, 34/1024 for q3) → 1024-shot estimates; `P(1|0)` = 0.0195 ± 0.0044 (binomial), so the ratio 2.65 is 1.8–4.8 at ±2σ of the denominator. The difference e − P(1|0) is still 5.8σ (shot noise of both), so "larger" stands; "2.65×" is not a robust number. (ii) e itself is 0.0517 ± 0.0036; no uncertainty is quoted.
- Supporting evidence I computed from the existing records (hypothesis, not proof — `r11_carryover.py`): the data qubits 2, 4, 16 end every BV shot in |1⟩ (answer bits). An initial |1⟩ on phys 4 / 16 / 2 gives exactly 0011 / 1001 / 1010 — the three single-bit wrong answers that the appendix leaves as "데이터 큐비트 쪽의 별도 손실". A maximum-likelihood fit of four initial-|1⟩ probabilities gives e₃ = 0.053, e₄ = 0.012, e₁₆ = 0.008, e₂ = 0.0025, success 0.8904 (measured 0.8888), χ² 0.3 over 8 bins (replay alone 17417; e₃ only 63). The same model on GHZ (targets q1..q4): e₁₂₁ = 0.027, e₁₂₃ = 0.009, e₁₃₆ = 0.005, e₁₂₂ = 0, success 0.9175 (measured 0.9205), χ² 725 → 26. Caveat: for a *data* qubit an initial flip is indistinguishable in Z-basis counts from extra readout P(0|1) or a Z error between its two H's, so only the ancilla's e is identified uniquely.
- Replace `appA_qpu_runs.md:120`:
  - Current: `**③ 보정값과의 대조.** 물리 3번은 이 작업에 쓰인 13개 큐비트 중 판독 오류가 가장 크다(0.0264). P(1\|0 준비)도 0.0195로 높다. 이 값에는 초기화 오류와 판독 오류가 함께 들어 있으므로, 그 시각의 초기화 오류는 많아야 약 2%였다. 맞춘 e = 5.17%는 그 2.65배다.`
  - Proposed: `**③ 보정값과의 대조.** 물리 3번은 이 작업에 쓰인 13개 큐비트 중 판독 오류가 가장 크다(0.0264). P(1\|0 준비)도 0.0195로 높다. 맞춘 e = 5.17%는 그 2.65배다. 하지만 두 값은 같은 양이 아니다. P(1\|0 준비)는 큐비트를 \|0⟩에 두는 보정 실험에서 1로 읽힌 비율이고, 초기화 오류와 판독 오류가 섞여 있다. 값이 1/1024의 배수(20/1024)라 1,024샷 추정으로 보이고, 그렇다면 통계 오차만 ±0.0044다. 반면 BV의 보조 큐비트는 매 샷이 \|−⟩(절반이 \|1⟩)로 끝나고 측정되지 않는다. 다음 샷은 그 들뜸을 반복 지연(rep\_delay)과 초기화 절차로 지운 뒤 시작한다. 직전 실험이 \|1⟩이었을 때 초기화 뒤에 남는 \|1⟩의 비율을 IBM은 `init_error`라는 별도 항목으로 공개하며, `prob_meas1_prep0`과 다른 양이라고 문서에 적는다. e와 견줄 값은 `init_error`의 절반쯤인데 이 실행에서는 받아 두지 않았다. 그러므로 e가 "보정값의 상한을 넘었다"고는 말할 수 없다.`
- Replace `appA_qpu_runs.md:122`:
  - Current: `**판정.** "보조 큐비트가 처음부터 |1⟩이었다"는 0000을 만드는 **유일한 단일 결함 기제**이고, 그 크기가 재현 대비 손실 7.2%p 중 5.1%p를 설명한다. 이 기제가 맞다는 근거는 강하다. 다만 크기가 보정값의 상한보다 크다. 그래서 그 37분 사이에 3번 큐비트의 초기화가 나빠졌는지(열적 들뜸, TLS 결함), 아니면 누설처럼 단일 파울리 결함 모형에 없는 기제인지는 이 자료로 가를 수 없다.`
  - Proposed: `**판정.** 0000은 "첫 CZ 때 보조 큐비트가 |1⟩이었다"는 것을 가리킨다. 회로 수준의 단일 결함 가운데 0000을 내는 것은 보조의 시작 상태 오류와 `cz[4, 3]` 하나의 상관 오류뿐이고, 뒤쪽은 보정된 CZ 오류율로는 0.06%밖에 만들지 못한다(①·②). 크기 5.1%p는 재현 대비 손실 7.2%p의 대부분이다. 남은 질문은 왜 |1⟩이었는가다. 후보는 셋이다. (가) 직전 샷의 들뜸이 덜 지워졌다. 보조는 매 샷 절반이 |1⟩인 채 측정 없이 끝나므로 보정 실험의 P(1\|0 준비)보다 커도 이상하지 않다. (나) 3번 큐비트의 열적 들뜸이나 리셋 품질이 그 시각에 나빴다. (다) 누설처럼 파울리 모형에 없는 기제다. 이 자료로는 셋을 가르지 못한다. 다만 같은 e를 3번 큐비트에 넣으면 곧이어 돈 HVA 작업의 실측−재현도 4.4σ에서 1.4σ로 줄어든다(6절). 서로 다른 두 작업이 같은 크기를 가리킨다.`
  ("4.4σ → 1.4σ" needs the one-line script addition in B4.)
- `CLAUDE.md:57`: `…설명되지만 크기가 보정값 상한의 2.65배라 확정하지 못했다(부록 A 5절).` → `…설명된다. 왜 |1⟩이었는지(직전 샷 들뜸의 잔류, 그 시각의 초기화 품질, 누설)는 가르지 못했다. 견줄 보정 항목은 P(1|0 준비)가 아니라 init_error 이고 아직 받지 않았다(부록 A 5절).`
- Report `§8.1:178`: `크기는 보정값 상한의 2.65배라 기제는 강하게 지지되지만 확정은 아니다.` → `크기 e = 5.17%는 보정값 P(1\|0 준비)의 2.65배지만 둘은 같은 양이 아니다(견줄 항목은 init\_error, 미수집). "첫 CZ 때 보조가 \|1⟩"이라는 기제는 강하게 지지되고, 그 원인은 가르지 못했다.`

### A3. Grover: the "pessimistic RB → depolarizing" conclusion is explained by a double count in the model
- `appA/appA_qpu_runs.md:72` (and `:78` "네 회로 중 재현이 가장 잘 맞았다"); model description `:98`; same construction in `ch14/ch14_hardware.md:29`.
- An RB error rate is the *total* average infidelity of the physical gate; decoherence during the gate is inside it. Both ch14's model and `Replay.gate` apply depolarizing(p from r) **and then** thermal relaxation for the gate duration (`appA_qpu_runs.py:135-137`). Per-gate size of the extra term for the Grover qubits (32 ns / 68 ns, my computation): thermal-relaxation infidelity 1.0–2.2×10⁻⁴ per 1q gate against calibrated sx r ≈ 2.1–2.2×10⁻⁴ (qubits 22, 23) — i.e. the 1q error is roughly counted 1.5×; qubit 24 has T1 = 68 µs, T2 = 78 µs.
- Numbers (my `r11_more.py`; identical when the same flag is added to the author's script, see patch): keeping idle relaxation but dropping gate-time relaxation gives replay GHZ 0.9591, BV 0.9635, **Grover 0.8090**, QPE 0.9449. Measured − replay becomes −3.9, −7.5, **+0.5 (0.8σ)**, **−2.3** %p. ALAP instead of ASAP changes little (0.9595 / 0.9627 / 0.7798 / 0.9330).
- So (a) the Grover "+3.5 %p, 5.7σ" discrepancy is a model artefact, not evidence about coherent cancellation or deep circuits; (b) "QPE matched best" is model-dependent; (c) the replay has a systematic uncertainty of 1–3 %p from modelling conventions alone, larger than the shot σ it is compared with.
- I could not open Qiskit Aer's source in this session; from memory, `NoiseModel.from_backend` chooses the depolarizing strength so that depolarizing ∘ thermal relaxation reproduces the reported gate error (the Aer API page I could open only says errors consist "of a depolarizing_error() followed by a thermal_relaxation_error()" with parameters "tuned … based on … the gate error and gate time"). If my recollection is right, this also explains part of "fake is the most optimistic" (B3).
- "CZ가 더 많은데도 좋았던 것은 그날 쓰인 연결이 좋았기 때문이다" is contradicted two sentences later by the script's own replay (good edges give only 0.7785).
- Replace in `appA_qpu_runs.md:72`:
  - Current: `CZ가 더 많은데도 좋았던 것은 그날 쓰인 연결이 좋았기 때문이다. 이 작업이 쓴 연결의 CZ 오류는 1.58×10⁻³\~2.46×10⁻³으로, 14장이 가정한 중앙값 3×10⁻³보다 작다. 그런데 당일 보정값으로 다시 돌려도 0.7785로, 실측이 3.5%p(5.7σ) 더 좋다. 깊은 회로에서는 "RB 오류율 → 탈분극" 모형이 비관적이라는 뜻이다. 원인 후보는 셋이다. RB 오류율에 섞인 결맞은 성분은 회로에 따라 상쇄될 수 있다. 게이트 시간(14장 가정 32·68 ns)과 실제 일정이 다를 수 있다. 보정 이후 오류율이 바뀌었을 수도 있다. 이 자료로는 가르지 못한다.`
  - Proposed: `CZ가 더 많은데도 좋았던 이유의 일부는 그날 쓰인 연결이 좋았던 데 있다. 이 작업이 쓴 연결의 CZ 오류는 1.58×10⁻³\~2.46×10⁻³으로, 14장이 가정한 중앙값 3×10⁻³보다 작다. 하지만 그 값으로 다시 돌린 재현은 0.7785로 14장 예측 0.759와 크게 다르지 않고, 실측은 그보다 3.5%p 더 좋다. 이 차의 가장 그럴듯한 원인은 재현 모형의 이중 계산이다. RB 오류율 r은 게이트 한 번의 오류 전체이고, 게이트가 걸리는 동안의 T1·T2 손실도 그 안에 들어 있다. 14장과 이 부록의 모형은 r로 만든 탈분극 위에 게이트 시간만큼의 열완화를 한 번 더 얹는다. 게이트 중 열완화만 빼고(유휴 시간의 열완화는 그대로) 다시 돌리면 재현이 0.8090이 되어 실측과 0.5%p 차이다(0.8σ). 같은 변형에서 QPE는 0.9449로 올라가 실측이 오히려 2.3%p 낮다. 재현값은 모형을 어떻게 짜느냐에 따라 1\~3%p 움직인다. 이 크기의 차를 장치의 성질로 읽어서는 안 된다. 그 밖의 후보(결맞은 오류의 상쇄, 가정한 게이트 시간 32·68 ns, 보정 뒤의 변화)는 이 자료로 가르지 못한다.`
- Replace in `appA_qpu_runs.md:78`: `네 회로 중 재현이 가장 잘 맞았다.` → `이 모형에서는 네 회로 중 차가 가장 작다. 게이트 중 열완화를 뺀 변형에서는 차가 −2.3%p로 벌어진다(3.4절).`
- Add to `appA_qpu_runs.md:98` after `…큐비트별 비대칭 판독 오류 P(1\|0 준비)·P(0\|1 준비)를 넣는다.`: `이 모형의 근사는 넷이다. 게이트 중 열완화는 r에 이미 들어 있는 몫을 한 번 더 센다(3.4절). 유휴 위상 손실에 쓴 T2는 IBM이 Hahn 에코로 잰 값이라 에코 없는 회로에서는 낙관적이다. 게이트 시간은 14장의 가정(32·68 ns)이다. 크로스토크·누설·초기화 오류는 없다.` (IBM docs: T2 "is reported from a Hahn echo sequence".)
- Suggested addition to `ch14_hardware.md` §7 (new paragraph after `:172`): `**게이트 중 열완화는 이중 계산이다.** RB 오류율 r에는 게이트 동안의 T1·T2 손실이 이미 들어 있다. 이 장의 모형은 r로 만든 탈분극에 게이트 시간만큼의 열완화를 더 얹으므로 CZ가 많은 회로일수록 비관적이다. 실제 ISA 회로와 당일 보정값으로 한 재현에서 Grover는 이 몫을 빼면 0.7785가 0.8090이 되어 실측 0.8137과 맞는다(부록 A 3.4절).`
- Script addition (verified in a scratch copy; full diff in `r11_appA_patch.diff`): add `gate_relax=True` to `Replay.__init__`, guard `for p in phys: self.relax(dt, p)` with `if self.gate_relax:`, and in the A3 loop print `sum(replay(g, cal, gate_relax=False).get(b, 0) for b in GOOD[k])`. Output: `0.9591 / 0.9635 / 0.8090 / 0.9449`, `실측−재현 -3.9 / -7.5 / +0.5 / -2.3 %p`.

### A4. ch14 §7: "실측은 BV를 빼면 둘 사이에 있었다" is false for HVA
- `ch14/ch14_hardware.md:178`. HVA: ch14 prediction −4.433, fake −4.581, measured −4.322 — worse than both. Also "다섯 회로 모두 이 장의 예측보다 좋았고" lists four numbers and compares fake HVA (no readout error) with the readout-included prediction.
- Current: `fake(`fake_marrakesh`)는 다섯 회로 모두 이 장의 예측보다 좋았고(GHZ 0.969, BV 0.978, Grover 0.874, QPE 0.953), 실측은 BV를 빼면 둘 사이에 있었다.`
- Proposed: `fake(`fake_marrakesh`)는 성공 확률을 잰 네 회로 모두 이 장의 예측보다 좋았고(GHZ 0.969, BV 0.978, Grover 0.874, QPE 0.953), 실측은 GHZ·Grover·QPE에서 둘 사이에 있었다. BV(0.889)와 HVA(−4.322)의 실측은 두 예측보다 나빴다. fake의 HVA 값 −4.581은 판독 오류가 없는 값이라 "판독 오류 제외" 예측 −4.555와 견줘야 한다.`

### A5. "BV의 손실 대부분" — 0000 is 46 % of BV's loss
- `appA/appA_qpu_runs.md:11`; `ch14/ch14_hardware.md:123` item (3). BV loss = 1 − 0.8888 = 0.1112 (445 shots); 0000 = 206 shots = 0.0515 = 46 %. "Most" is true only of the shortfall against the replay (5.1 of 7.2 %p = 71 %), which is what `:122` says correctly.
- `appA_qpu_runs.md:11` current: `BV의 손실 대부분은 측정하지 않는 보조 큐비트 하나가 |1⟩에서 시작한 탓으로 설명된다. 그 큐비트는 HVA 에너지 측정에도 들어갔다.`
- Proposed: `BV에서 가장 큰 오답 0000(전체 손실 0.1112 중 0.0515, 재현보다 더 잃은 7.2%p 중 5.1%p)은 측정하지 않는 보조 큐비트 하나가 |1⟩에서 시작했을 때의 모양과 같다. 그 큐비트는 HVA 에너지 측정에도 들어갔다.`
- `ch14_hardware.md:123` current: `(3) BV의 손실은 대부분 오답 `0000`(5.2%)이고, 이것은 측정하지 않는 보조 큐비트가 |1⟩에서 시작했을 때만 나오는 지문이다.`
- Proposed: `(3) BV 손실 11.1%p 가운데 가장 큰 몫은 오답 `0000`(5.15%p)이고, 이것은 측정하지 않는 보조 큐비트가 |1⟩에서 시작했을 때 나오는 지문이다(실제 ISA 회로에 단일 결함을 넣어 확인 — 부록 A 5절).`

### A6. "유일한 단일 결함 기제" / "12 of 345" — the enumeration does not show uniqueness
- `appA/appA_qpu_runs.md:116`, `:122`, `:169` (glossary "BV 345가지 중 12가지"); `ch14_hardware.md:123` ("…했을 때만 나오는 지문"); report `§8.1:178` ("…경우만 만든다"); `appA_qpu_runs.py:205-219`.
- I reproduced 345 / 12 exactly. But: (i) the script inserts a Pauli on *every* qubit after *every entry of the serialized gate list* (23 positions × 5 qubits × 3). The 12 hits are 6 list positions (start + the five entries that happen to precede `cz[4, 3]` in the list, all on other qubits) × {X, Y}. The ancilla has no gate before its first CZ, so these are **one** physical fault (X and Y act identically on |0⟩ up to phase). "12 of 345" has no physical meaning. (ii) The enumeration contains no two-qubit faults. A standard circuit-level count (initial state 3 per qubit; 3 after each sx/x; 15 after each cz; rz virtual) has 84 faults, of which **6** give 0000 with certainty: initial X/Y on the ancilla and the four weight-2 Paulis after `cz[4, 3]` (X or Y on the ancilla with X or Z on phys 4). Those four are exactly the source of the replay's own P(0000) = 0.00060 (4/15 × 1.25 × 1.76×10⁻³ = 5.9×10⁻⁴) — so the script's replay number contradicts "only". (iii) "No gates applied at all" also gives 0000. The existing data argues against a job-wide version of that: QPE's 000 is 0.0225 measured vs 0.0223 in my replay, and GHZ's 00000 − 11111 imbalance is 76 counts against 44 expected (0.5σ excess; a 5 % no-pulse fraction would add ≈ 200).
- Useful physics the text could state (verified by simulation): in the ISA the first CZ acts on an ancilla in |0⟩, i.e. it is the identity in the ideal circuit; q0's phase flip was compiled into a virtual Z. That is why an ancilla starting in |1⟩ un-does q0's flip *and* disables the other two. Skipping the ancilla's X gives 1000, not 0000 (observed 0.27 %), so a failed X is not the cause.
- Replace `appA_qpu_runs.md:116`:
  - Current: `**① 단일 결함 열거.** 모든 위치(시작 직후, 각 게이트 뒤)와 모든 큐비트에 파울리 X·Y·Z 하나를 넣는다. 모두 345가지다. 이 가운데 0000을 만드는 결함은 12가지뿐이다. 모두 **보조 큐비트 3번의 X 또는 Y이고, 첫 CZ `cz[4, 3]` 이전**이며, 각각 확률 1로 0000을 낸다. 이것은 "보조가 |1⟩에서 시작한 것"과 같다. 그러면 첫 CZ가 q0에 위상을 주고, X가 보조를 |0⟩으로 돌려 나머지 둘은 위상을 받지 못한다. 결과적으로 세 비트가 한꺼번에 뒤집힌다. 판독 오류나 다른 큐비트의 오류로는 0000이 나오지 않는다.`
  - Proposed: `**① 단일 결함 열거.** 게이트 목록의 모든 자리(시작 직후, 각 게이트 뒤)에서 모든 큐비트에 파울리 X·Y·Z 하나를 넣는다. 모두 345가지다. 이 가운데 0000을 만드는 것은 12가지이고, 모두 **보조 큐비트 3번의 X 또는 Y이며 첫 CZ `cz[4, 3]` 이전**이다. 각각 확률 1로 0000을 낸다. 12는 목록에서 첫 CZ 앞에 놓인 자리 6곳에 X·Y를 넣은 수일 뿐이다. 보조에는 그 앞에 게이트가 없으므로 물리적으로는 모두 "보조가 |1⟩에서 시작했다"는 결함 하나다. 이상적인 ISA 회로에서 첫 CZ는 보조가 |0⟩이라 아무 일도 하지 않고, q0의 위상은 가상 Z로 미리 들어가 있다. 보조가 |1⟩에서 시작하면 첫 CZ가 그 위상을 되돌리고, X가 보조를 |0⟩으로 돌려 나머지 둘도 위상을 받지 못한다. 그래서 세 비트가 한꺼번에 틀린다. 2큐비트 게이트의 상관 오류까지 넣어 회로 수준으로 세면 결함은 84가지이고 0000을 내는 것은 6가지다. 보조의 시작 상태 X·Y와, `cz[4, 3]` 직후 보조의 X·Y가 q0 쪽의 X·Z와 함께 생기는 4가지다. 뒤의 4가지가 ②의 재현값 0.00060을 만든다. 1큐비트 게이트 오류, 판독 오류, 다른 CZ의 오류 하나로는 0000이 나오지 않는다. 보조의 X가 빠지면 0000이 아니라 1000이 나온다.`
- `appA_qpu_runs.md:169` glossary: `| 단일 결함 열거 | 모든 위치에 오류 하나씩만 넣어 어떤 오답이 나오는지 전부 시험하는 방법 | BV 345가지 중 12가지 |` → `| 단일 결함 열거 | 모든 위치에 오류 하나씩만 넣어 어떤 오답이 나오는지 전부 시험하는 방법 | BV: 0000을 내는 것은 보조의 시작 상태 오류와 첫 CZ의 상관 오류뿐 |`
- Report `§8.1:178`: `단일 결함 345가지 중 보조 큐비트(물리 3번)가 \|1⟩ 에서 시작한 경우만 만든다.` → `단일 파울리 결함 가운데 보조 큐비트(물리 3번)가 \|1⟩ 에서 시작한 경우만 만든다(회로 수준으로는 첫 CZ 의 상관 오류도 만들지만 보정된 오류율로는 0.06%).`
- Script: the circuit-level enumeration block is in `r11_appA_patch.diff` (prints `84가지 중 … 6가지` and the list).

### A7. The one-line summary attributes all deviations to "one specific qubit"
- `appA/appA_qpu_runs.md:11`: `**14장의 예측은 회로마다 −3.6\~+5.5%p 안에서 맞았고, 어긋난 곳은 "사양서의 평균"이 아니라 "특정 큐비트 하나"였다.**` The two halves refer to different comparisons. The −3.6…+5.5 %p deviations are against the *ch14 prediction*, and the two largest (Grover +5.5, QPE +4.1) are not single-qubit effects at all (different CZ counts, better readout that day, and the A3 double count). "One specific qubit" describes only the GHZ/BV shortfall against the *replay*.
- Proposed: `**14장의 예측은 회로마다 −3.6\~+5.5%p 벗어났다. 예측보다 좋았던 Grover·QPE는 트랜스파일 결과와 그날의 보정값이 가정과 달랐던 탓이고, 당일 보정값으로도 설명되지 않는 손실(GHZ 3.5%p, BV 7.2%p)은 "사양서의 평균"이 아니라 큐비트 하나씩에 몰려 있었다.**`

## 4. Grade B findings

### B1. "차가 큰 두 회로는 모두 특정 큐비트 하나로 설명된다" (`appA_qpu_runs.md:104`)
Grover's |measured − replay| (3.53 %p) is larger than GHZ's (3.47 %p). And "설명된다" for GHZ means "localized", not "explained" (the same paragraph in §3.2 says the cause cannot be determined). Fixed by the A1 replacement text.

### B2. σ against a model counts shot noise only; say so, and give the model's own uncertainty
- `appA_qpu_runs.md:38` (2.8σ), `:72` (5.7σ), `:78` (2.3σ), `:137` (4.4σ), `:153` (glossary); stdout `-8.1σ`, `-14.6σ` (`appA_qpu_runs.py:183`); `ch14_hardware.md:123` (2.8σ).
- All formulas and N are right (binomial SE with N = 4000 on the summed success outcomes; Bell uses N = 997 conditional; HVA uses the Estimator's std). But the denominator is the measurement's shot noise alone. Against a noise-model value this tells you only that the gap is not shot noise; it is not a significance of anything, because the model side has (a) calibration sampling error — propagating the 1024-shot readout values alone gives ±0.35 / ±0.47 / ±0.34 / ±0.26 %p on the four replay values (Monte Carlo, `r11_more.py`), which already turns −8.1 / −14.6 / +5.7 / −2.3 into −6.3 / −10.6 / +5.0 / −2.0; and (b) modelling-convention shifts of 1–3 %p (A3). Conversely a small σ does not validate the model (GHZ vs ch14 prediction is +0.6σ by cancellation, as §4 item 1 says).
- `appA_qpu_runs.md:153` current: `| σ (표준편차 단위) | 차이 ÷ 표준오차. 2σ 안은 우연과 구분이 어렵고, 3σ 넘으면 우연으로 보기 어렵다 | HVA 실측 − 예측 = 2.8σ |`
- Proposed: `| σ (표준편차 단위) | 차이 ÷ 표준오차. 2σ 안은 우연과 구분이 어렵고, 3σ 넘으면 우연으로 보기 어렵다. 이 부록의 σ는 실측의 샷 잡음만 센다. 모형 쪽의 불확도(보정값의 통계 오차, 모형을 짜는 방식에 따른 1\~3%p)는 들어 있지 않으므로, 모형과의 차에 붙은 σ는 "샷 잡음은 아니다"라는 뜻일 뿐이다 | HVA 실측 − 예측 = 2.8σ |`
- Script: change the label at `appA_qpu_runs.py:183` from `σ` to `σ(샷 잡음만)`.

### B3. fake HVA is a readout-free number shown unlabelled next to readout-included ones; "fake snapshot was better than the real device" is an untested inference
- `appA_qpu_runs.md:38`, `:96`, `:133` (tables), `:103` (text), `:171` (glossary "다섯 회로 모두 가장 낙관적").
- fake HVA −4.5807 (±0.02 sampling) should be compared with the readout-excluded values: ch14 −4.5551, replay −4.5816 — it is *not* an outlier there. In the §6 table "fake 90 %" vs "replay 76 %" vs "measured 50 %" mixes the two kinds.
- "가짜 장치의 보정 스냅숏이 그날의 실제 장치보다 좋았다" was not checked (the snapshot values were never read; different device, different layout, and a differently built noise model — see A3).
- `:103` current: `2. **fake는 다섯 회로 모두 가장 낙관적이었다.** 가짜 장치의 보정 스냅숏이 그날의 실제 장치보다 좋았다. fake의 HVA 값은 Aer가 판독 오류 없는 정확값에 N(0, 0.02) 표본 잡음을 더한 것이다. 그래서 다시 돌릴 때마다 ±0.02쯤 움직인다.`
- Proposed: `2. **fake는 성공 확률을 잰 네 회로 모두에서 가장 낙관적이었다.** 원인은 확인하지 않았다. 가짜 장치(다른 장치인 marrakesh)의 보정 스냅숏이 더 좋았을 수도 있고, Aer의 잡음 모형이 이 부록의 재현과 다르게 짜여 있기 때문일 수도 있다(3.4절). fake의 HVA 값은 판독 오류를 넣지 않은 기댓값에 Aer가 N(0, 0.02) 표본 잡음을 더한 것이다. 그래서 다시 돌릴 때마다 ±0.02쯤 움직이고, 견줄 상대는 판독 오류를 뺀 값(14장 −4.5551, 재현 −4.5816)이다.`
- Tables: in `:38`, `:96`, `:133` write the fake HVA cell as `−4.5807 (판독 오류 없음)`; in `:133` the "이상 대비" cell as `90% (판독 오류 없음)`. Glossary `:171`: `다섯 회로 모두 가장 낙관적` → `네 회로 모두 가장 낙관적(HVA는 판독 오류 없는 값)`.

### B4. HVA §6: compute the "speculation" — it takes one line
- `appA_qpu_runs.md:137`: `실측은 재현보다도 0.171(4.4σ) 높다. BV에서 본 3번 큐비트의 초기화 문제가 여기에도 있었다면 이 방향의 차를 만든다. 다만 이것은 추측이고 검증하지 않았다.`
- With the existing hook: `energy(calh, init_exc={anc[0]: e_fit})` = −4.3774 → measured − replay = +0.055 = **+1.4σ** (my independent code: −4.3775; on the actual ISA gate list: −4.3811, +1.5σ). e on q3 alone that would close the gap completely: 0.076. Sensitivity dE/de ≈ 2.2 per qubit.
- Proposed: `실측은 재현보다도 0.171(4.4σ) 높다. BV에서 맞춘 e = 0.0517을 3번 큐비트의 시작 상태에 넣으면 재현이 −4.3774가 되어 차가 0.055(1.4σ)로 줄어든다. 서로 다른 작업(Sampler와 Estimator)이 같은 크기의 시작 오류를 가리키는 셈이다. 다만 HVA는 에너지 하나만 주므로 다른 원인과 가를 수는 없다.`
- Code (after `en_no3 = …`, verified): `en_e = energy(calh, init_exc={anc[0]: e_fit})` and a print of `en_e`, `ev - en_e`, `(ev - en_e) / sd`.
- Also checked: the author's "logical reconstruction" shortcut is fine — replaying the actual ISA prefix with fitted rz·sx·rz·sx·rz basis blocks gives −4.4967 against −4.4930 (0.1σ). The claim in the code comment ("에너지에 주는 차는 작다") holds.

### B5. The proposed confirmation experiment would not discriminate the alternatives
- `appA_qpu_runs.md:124`; report `§8:162`; `CLAUDE.md:58`.
- (a) "보조만 다른 큐비트에" is not realizable on heavy-hex without SWAPs (phys 2, 4, 16 are neighbours of 3 only); the transpiler would produce a different circuit. "0000이 사라지면 3번 큐비트 문제로 확정" does not follow: without the original layout as a same-job control, a disappearance can equally be time (the condition ended). And if 0000 *persists* on another ancilla — the carry-over outcome — the text has no reading for it. It also cannot separate initialization from leakage or from a correlated `cz[4, 3]` fault.
- (b) "prepare-and-measure on qubit 3" measures the same SPAM sum as the calibration, in a sequence where the qubit is never left excited — it would most likely return ≈ 2 % and look like a refutation.
- Proposed replacement for `:124`:

```
**확정하려면.** 먼저 QPU를 쓰지 않는 확인이다. `env/fetch_ibm_job.py`가 작업 시각의 `init_error`(직전 실험이 |1⟩이었을 때 초기화 뒤에 남는 |1⟩의 비율)와 항목별 측정 시각, 작업 옵션(`rep_delay`, `init_qubits`)까지 받게 고쳐 새 파일로 받는다. 3번 큐비트의 `init_error`가 10% 안팎이면 e ≈ 0.5 × `init_error`로 (가)가 설명된다.

그다음이 QPU 약 15초가 드는 실험이다(승인 후). 2026-11-11에 돈 ISA 게이트 목록을 같은 물리 큐비트에 그대로 올리고(최적화 수준 0, 다시 합성하지 않는다), 한 작업의 PUB로 묶어 같은 시각에 잰다. 각 4,000샷.

- P0 대조: 그때의 BV 그대로. 0000이 다시 5% 안팎인지 본다. 재현 수준(0.06%)으로 내려가 있으면 그날의 일시적 상태였고 아래 비교는 뜻이 없다.
- P1: P0에서 보조의 마지막 H를 X로 바꿔 보조가 |0⟩으로 끝나게 한다. 0000이 사라지면 (가) 직전 샷의 잔류이고, 그대로면 (나)나 (다)다.
- P2: P0의 맨 앞에 보조 `reset`을 두 번 넣는다. 0000이 줄면 |1⟩ 들뜸이고((가)·(나)), 줄지 않으면 (다) 누설 쪽이다.
- P3: P0에서 게이트가 하나도 없던 물리 114번에 X를 넣어 표지로 쓴다. 게이트가 걸리지 않은 샷이라면 0000과 함께 표지도 0이고, 보조의 시작 오류라면 표지는 1이다.
- P4: 3번 큐비트만으로 "측정 → 측정 → X"와 "측정 → 측정" 두 회로. 첫 측정의 1 비율이 직전 샷이 |1⟩로 끝났을 때와 |0⟩으로 끝났을 때의 초기화 오류다. 연달아 두 번 재므로 판독 오류와 나눌 수 있다.
- P5: 같은 BV를 다른 허브 큐비트와 그 이웃 셋에 올린다. heavy-hex에서는 보조만 옮기면 SWAP이 들어가 다른 회로가 되므로 네 큐비트를 함께 옮긴다. 0000이 여기서도 나오면 3번 큐비트가 아니라 "들뜬 채 측정 없이 끝나는 보조"라는 회로 구조의 문제다.
- 별도 작업 하나: P0를 `rep_delay`를 장치가 허용하는 최댓값으로 늘려 다시 잰다. 0000이 지연만큼 줄면 (가), 그대로면 (나)다.
```
- "약 15초" is my estimate from the 2026-11-11 rate (16,000 shots in 6 s; 7 PUBs + 1 = 32,000 shots ≈ 12 s plus a second job's overhead; a longer `rep_delay` costs proportionally more). `reset` is listed among Heron's supported instructions in IBM's docs; I could not check that the Open plan / this backend accepts mid-circuit measurement followed by gates (P4).
- Report `§8:162` and `CLAUDE.md:58` should point to the new list (`보조만 다른 큐비트에 두고 재실행, 3번 큐비트의 준비-측정 오류 직접 측정` → `init_error 조회(QPU 불필요) 뒤, 같은 ISA 회로의 변형 여섯 개를 한 작업으로`).

### B6. "버전 주의" / Runtime interface — what could and could not be confirmed
- `ch14_hardware.md:107`; `appA_qpu_runs.md:143`, `:177`, `:194`; `ch14_submit.py:5-6`; `run_ibm.py:15-16`; report `§6-2`; `requirements/sdk.txt:5`.
- **Confirmed** (IBM release notes page): 0.50.0 is dated 2026-09-24; its deprecation note reads "EstimatorV2 and SamplerV2 are deprecated. Please use the new client-side implementations of these primitive interfaces: Estimator and Sampler instead. … At the end of the deprecation period, the legacy server-side implementations will be removed, and the top-level imports (qiskit_ibm_runtime.Sampler[V2] and qiskit_ibm_runtime.Estimator[V2]) will be updated to point to the client-side implementations."; new features: "Added the client-side primitives qiskit_ibm_runtime.executor_sampler.Sampler and qiskit_ibm_runtime.executor_estimator.Estimator"; they "perform pre- and post-processing on the client side and delegate execution to the … Executor". 0.51.0 is dated 2026-10-22.
- **Not confirmed**: that in 0.52.0 the top-level `Sampler`/`Estimator`/`SamplerV2`/`EstimatorV2` already point to the client-side classes, and that a module named `qiskit_ibm_runtime.server_side` exists. The release-notes page's newest entry today is 0.51.0 (no 0.52 notes; the string "server_side" does not occur); the 0.50 note puts the re-pointing "at the end of the deprecation period"; and IBM's current options guides still label `from qiskit_ibm_runtime import Sampler` as "Legacy server-side primitive before qiskit-ibm-runtime 0.50.0" and `from qiskit_ibm_runtime.executor_sampler import Sampler` as the client-side import. The author's claim rests on reading the installed 0.52.0 `__init__.py`, which I could not inspect (no SDK here). Indirect support in the repo: the frozen job inputs are executor-style (`inputs['quantum_program'].items[i].circuit`, parametrised basis block with separately stored values), and `samplomatic==0.21.0` is in the lock.
- Two wording problems: (i) `appA_qpu_runs.md:143` says the earlier sentence "0.50부터 SamplerV2·EstimatorV2가 폐기 예정" "은 정확하지 않았다" — but that is literally the release note's first sentence. (ii) "사용 시간의 키 이름은 0.51에서 바뀌었으므로" (`ch14_hardware.md:107`, snippet `:102`, `ch14_submit.py:10`): the 0.51.0 note is a *bug fix for local testing mode* — "In local testing mode, job.metrics() now returns the same layout as jobs run on IBM Quantum Platform: … usage contains qpu_charge_time_seconds and status instead of quantum_seconds and seconds." For real jobs the key comes from the platform, not from the client version.
- Recommendation: import the client-side classes by their explicit, documented paths (`from qiskit_ibm_runtime.executor_sampler import Sampler`, `from qiskit_ibm_runtime.executor_estimator import Estimator`) in `ch14_submit.py`, `run_ibm.py` and the ch14 snippet, and write `type(sampler).__module__` into the run record. Then the statement "this run used the client-side implementation" is evidence in the record instead of a reading of `__init__.py`, and it holds for 0.50–0.53 whatever the top-level names point to.
- `ch14_hardware.md:107` current: `**버전 주의.** qiskit-ibm-runtime 0.50(2026-09-24)부터 폐기 예정(deprecated)인 것은 서버 쪽 `SamplerV2`·`EstimatorV2`(`qiskit_ibm_runtime.server_side`)다. 0.52의 최상위 `Sampler`·`Estimator`는 이미 클라이언트 쪽 구현(`executor_sampler.Sampler`, `executor_estimator.Estimator`)을 가리키고, 최상위 `SamplerV2`·`EstimatorV2`도 같은 클래스의 별칭이다(설치된 패키지의 `__init__.py`로 확인). 그래서 이 장의 코드는 최상위 `Sampler`·`Estimator`를 쓰고 runtime ≥ 0.52를 전제로 한다. 사용 시간의 키 이름은 0.51에서 바뀌었으므로(`qpu_charge_time_seconds`) 키 대신 `job.usage()`를 쓴다.`
- Proposed: `**버전 주의.** qiskit-ibm-runtime 0.50(2026-09-24) 릴리스 노트는 "EstimatorV2 and SamplerV2 are deprecated"라고 적는다. 폐기되는 것은 서버 쪽 구현이고, 대체는 클라이언트 쪽 구현 `qiskit_ibm_runtime.executor_sampler.Sampler`·`executor_estimator.Estimator`다. 노트는 폐기 기간이 끝나면 최상위 이름이 클라이언트 쪽 구현을 가리키게 된다고 적는다. 이 책이 설치한 0.52.0에서는 최상위 `Sampler`·`Estimator`(와 `SamplerV2`·`EstimatorV2`)가 이미 그 클래스였다(설치된 패키지의 `__init__.py`로 확인. IBM의 0.52 릴리스 노트로는 확인하지 못했다). 최상위 이름이 무엇을 가리키는지는 버전에 따라 다르므로, 구현을 못 박으려면 `executor_sampler`·`executor_estimator` 경로에서 직접 가져온다. 사용 시간은 `job.usage()`로 읽는다. `job.metrics()['usage']`의 키는 플랫폼 작업에서 `qpu_charge_time_seconds`이고, 로컬 시험 모드는 0.51부터 같은 꼴을 낸다.`
- `appA_qpu_runs.md:143` current first sentences: `"0.50부터 `SamplerV2`·`EstimatorV2`가 폐기 예정"은 정확하지 않았다. 폐기 예정인 것은 서버 쪽 클래스(`qiskit_ibm_runtime.server_side`)다.` → `"0.50부터 `SamplerV2`·`EstimatorV2`가 폐기 예정"은 릴리스 노트의 문장 그대로지만, 폐기되는 것은 그 이름이 아니라 서버 쪽 구현이다(설치된 0.52.0에서는 `qiskit_ibm_runtime.server_side`).`

### B7. "실행 완료(2026-11-11, fake·real)" overstates what the current `ch14_submit.py` has run
- `ch14_submit.py:1`; `ch14_hardware.md:77` ("**Qiskit Runtime 제출 스크립트(실행 완료, 2026-11-11).**"), `:80` ("같은 호출을 … 실행").
- The 2026-11-11 hardware run used the pre-merge script (names `SamplerV2`/`EstimatorV2`, GHZ `measure_all`, no early dump, no stats). The current `real --yes` path — `record["job_ids"]`, `pub.data.c` for the re-mapped GHZ, `getattr(d, "stds")`, the `metrics()` lookup, the frozen-copy write — has only been exercised up to the plan (`real` without `--yes`) and in `fake` mode with Aer primitives.
- `ch14_submit.py:1` proposed: `"""14장 — IBM Open Plan 제출 스크립트. qiskit==2.6.2, qiskit-ibm-runtime==0.52.0, qiskit-aer==0.18.2 — 실행 완료(2026-11-11: fake 전체, real 은 계획 출력까지). 실기 제출은 이전 판(SamplerV2·EstimatorV2 이름, GHZ measure_all)으로 2026-11-11 에 했고, 지금의 제출·기록 부분은 실기에서 아직 돌리지 않았다`
- `ch14_hardware.md:77` proposed heading: `**Qiskit Runtime 제출 스크립트(2026-11-11 실기 실행은 이전 판으로, 지금 판은 `fake`와 `real` 계획 출력까지 실행).**`; and at the end of that paragraph replace `기록(`ch14_record_<모드>.json`)에는 두 작업 id, 회로별 2큐비트 수·깊이·물리 큐비트, HVA 표준오차가 남는다.` with `지금 판의 기록(`ch14_record_<모드>.json`)에는 두 작업 id, 회로별 2큐비트 수·깊이·물리 큐비트, HVA 표준오차가 남는다. 2026-11-11 기록은 이전 판이 만든 것이라 이 항목이 없고, 표준오차와 실기 CZ 수는 `runs/20261111-063819Z_ibm_fez_ch14_isa_calib.json`에서 읽는다.`
- `ch14_hardware.md:80` snippet comment: `(같은 호출을 ch14_submit.py real --yes 로 2026-11-11 실행)` → `(2026-11-11 실기 실행은 같은 호출을 SamplerV2·EstimatorV2 이름으로 했다)`.

### B8. `ch14_hardware.py` block 7 reads a mutable file, skips silently, and has no check
- `ch14_hardware.py:407-408` reads `ch14/ch14_record_real.json` (CLAUDE.md: "다음 실행이 덮어쓴다") and does nothing if it is absent. The section-5 table therefore (a) silently changes source at the next `real --yes`, (b) has no source and still "passes" `test_all.sh` if the file is missing, and (c) crashes with `KeyError: 'ghz5'` while a new job is pending, because `ch14_submit.py:122` rewrites the file with an empty `jobs` before any result exists. The task's question "do they read only frozen records in runs/" → **appA yes, ch14 block 7 no.**
- Replacement for lines 406-409 (then dedent the block by one level):
```python
_rec = _here.parent / "runs" / "20261111-063819Z_ibm_fez_ch14_record.json"   # 고정 사본. ch14/ch14_record_real.json 은 다음 real --yes 가 덮어쓴다
rec = json.loads(_rec.read_text(encoding="utf-8"))                            # 없으면 여기서 멈춘다 — 5절 표의 출처가 없다는 뜻이다
assert (rec["backend"], rec["shots"]) == ("ibm_fez", 4000)
assert rec["jobs"]["ghz5"]["job_id"] == "dbq0qmoqtr0c738f53u0" and rec["jobs"]["hva"]["job_id"] == "dbq0qpk4148s739vbffg"   # 본문 5절이 적은 작업
```
- Also add after the loop: `assert abs(s_all["GHZ"] - 0.9205) < 5e-5` style checks for the five table values (ch14 currently has no assert at all — report §6-3), or at least for the four success counts (3682, 3555, 3255, 3689 of 4000).
- `ch14_hardware.md:113` then reads correctly if `ch14_record_real.json` is replaced by `runs/20261111-063819Z_ibm_fez_ch14_record.json`; also add `(HVA의 표준오차는 `runs/20261111-063819Z_ibm_fez_ch14_isa_calib.json`)` — today the caption attributes ±0.039 to a file that does not contain it.

### B9. `ch14_submit.py`: ids and results can be left only in memory / on the console
Safety is fine: `real` without `--yes` returns at line 106 before any `.run()`; a typo in the mode falls back to `fake`; no path spends QPU time unintentionally. The bit-order conversion is right and applied once (`measure(range(n), range(n)[::-1])` for GHZ/BV/Grover, direct for QPE, `"bit_order": "book"` written, block 7 flips GHZ only for records without the key). Depth uses the shared `gate_depth`. Remaining defects (untested replacements — no SDK here):
1. `:114-122` — the sampler job is submitted, then `estimator.run` is called, and only then is anything written. If `estimator.run` raises, a job that costs QPU time exists and its id is only on the console.
2. `:124-134` — Sampler counts are fetched, then `ejob.result()` is awaited, and only then `dump()`. If the Estimator job fails or the wait is interrupted, the counts already received are not on disk, and nothing in the kit can rebuild them (`fetch_ibm_job.py` saves no counts; `run_ibm.py --fetch` handles only Bell/GHZ).
3. `:111`, `:122` — `dump()` writes `ch14_record_real.json` relative to the current directory (run from the repo root → the file lands in the wrong place) and overwrites the previous record before any result exists (see B8).
4. `:147-149` — the frozen copy in `runs/` is written only at the very end.
5. `:138-141` — no retry; `usage()` can be 0 before aggregation (the script's own comment says so) and 0 is recorded silently. `run_ibm.final_usage` already solves this.
6. `:108` — `"t"` is local time without zone ("2026-11-11 15:38"); the frozen-copy name uses UTC at *completion*, while the existing 2026-11-11 file is named by job creation time.

Replacement for lines 108-149:
```python
    here = pathlib.Path(__file__).resolve().parent
    stamp = time.strftime("%Y%m%d-%H%M%SZ", time.gmtime())           # 제출 시각(UTC) — 고정 사본의 이름
    record = {"backend": backend.name, "shots": shots, "t_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
              "bit_order": "book", "primitives": [type(sampler).__module__, type(estimator).__module__], "jobs": {}}

    def dump(final=False):
        text = json.dumps(record, indent=1, ensure_ascii=False)
        if mode == "fake":
            (here / "ch14_record_fake.json").write_text(text, encoding="utf-8")
            return
        (here.parent / "runs" / f"{stamp}_{backend.name}_ch14_record.json").write_text(text, encoding="utf-8")   # 처음부터 고정 사본에 쓴다
        if final:                                         # ch14/ 의 기록은 결과가 다 온 뒤에만 바꾼다
            (here / "ch14_record_real.json").write_text(text, encoding="utf-8")

    job = sampler.run(list(isa.values()), shots=shots)
    jid = job.job_id() if hasattr(job, "job_id") else "fake"
    print("sampler job id:", jid, flush=True)
    record["job_ids"] = {"sampler": jid}; dump()          # Estimator 제출이 실패해도 이 id 는 파일에 있다
    ejob = estimator.run([(hva_isa, H_isa)], precision=0.02)
    ejid = ejob.job_id() if hasattr(ejob, "job_id") else "fake"
    print("estimator job id:", ejid, flush=True)
    record["job_ids"]["estimator"] = ejid; dump()

    for k, pub in zip(isa.keys(), job.result()):
        counts = pub.data.c.get_counts()
        record["jobs"][k] = {"job_id": jid, "counts": counts, **st[k]}
        print(k, sorted(counts.items(), key=lambda kv: -kv[1])[:3])
    dump()                                                # Estimator 결과를 기다리다 끊겨도 counts 는 남는다
    d = ejob.result()[0].data
    ev, sd = float(d.evs), float(np.asarray(getattr(d, "stds", np.nan)))
    print(f"HVA <H> = {ev:.4f} ± {sd:.4f} (이상 -4.6462 | 노이즈 예측: 판독 오류 포함 -4.433, 제외 -4.555)")
    record["jobs"]["hva"] = {"job_id": ejid, "evs": ev, "stds": sd, **st["hva"]}
    dump()
    if mode != "fake":
        from run_ibm import final_usage                   # 집계 전의 0 을 기록하지 않는다(None = 미확정)
        record["qpu_seconds"] = [final_usage(j) for j in (job, ejob)]
        print("qpu_seconds:", record["qpu_seconds"])
    dump(final=True)
    return record
```

### B10. `fetch_ibm_job.py`: hangs or aborts on an unfinished job, can overwrite a frozen record, and drops what the analysis needs
It uses no QPU time (only `service.job`, `job.inputs`, `job.status`, `job.result`, `backend.properties(datetime=…)`) and writes no token/CRN into the JSON or stdout (the only account-related output path is `run_ibm.get_service`'s error message, which prints the first 200 characters of the exception — I could not check whether an "invalid instance" error echoes the CRN). Defects:
1. `:68` `job.result()` is called unconditionally: for a queued/running job it **blocks** until the job ends; for a failed/cancelled job it raises, and because the file is written only after the loop (`:72-73`), everything already fetched for the other ids is lost.
2. `:72-73` overwrites an existing `--out` file. The appendix's own §9 command would silently replace a frozen record ("`runs/` 의 기록은 고치지 않는다").
3. `:39` discards each property's measurement time (A1); `init_error`, `readout_length`, gate lengths, job options (`rep_delay`, `init_qubits`, resilience) and job timestamps are not fetched (A2, A3); Sampler counts are not saved (B9-2), so appA cannot assert that the record's counts are the job's counts.
4. `:62` assumes executor-style inputs (`quantum_program`); a legacy job gives a bare `KeyError`.
5. `:74` `p.relative_to(ROOT)` raises after a successful write if `--out` is absolute.

Replacement sketch (untested), keeping every existing key so `appA_qpu_runs.py` still reads the file:
```python
def calib(props, qubits, gates):
    q = {}
    for k in sorted(qubits):
        def qp(name):
            try:
                v, t = props.qubit_property(k, name); return v, str(t)
            except Exception:
                return None, None
        vals = {n: qp(n) for n in ("T1", "T2", "readout_error", "prob_meas1_prep0", "prob_meas0_prep1", "init_error", "readout_length")}
        q[k] = {"T1_us": props.t1(k) * 1e6, "T2_us": props.t2(k) * 1e6, "readout_error": props.readout_error(k),
                "p_meas1_prep0": vals["prob_meas1_prep0"][0], "p_meas0_prep1": vals["prob_meas0_prep1"][0],
                "init_error": vals["init_error"][0], "readout_length_s": vals["readout_length"][0],
                "measured_at": {n: t for n, (_, t) in vals.items()}}          # 항목별 측정 시각
    g = {}; glen = {}
    for name, qs in sorted({(n, tuple(qs)) for n, qs, _, _ in gates if n in ("sx", "x", "cz")}):
        key = f"{name}{list(qs)}"
        try:
            g[key] = props.gate_error(name, list(qs)); glen[key] = props.gate_length(name, list(qs))
        except Exception:
            g[key] = glen[key] = None
    return {"last_update_date": str(props.last_update_date), "qubits": q, "gates": g, "gate_length_s": glen}
```
and in `main()`:
```python
    p = ROOT / a.out
    if p.exists():
        sys.exit(f"{a.out} 이 이미 있다. runs/ 의 기록은 고치지 않는다 — 다른 이름을 준다.")
    ...
        st = str(job.status())
        if "DONE" not in st.upper():          # 끝나지 않은 작업은 기다리지 않는다(result() 는 끝날 때까지 막는다)
            out[jid] = {"backend": be.name, "created": str(job.creation_date), "status": st}
            print(f"{jid}: {st} — 결과와 회로는 받지 않았다"); continue
        ...
        res = job.result(); data = res[0].data
        if hasattr(data, "evs"):
            out[jid]["estimator"] = {"evs": float(data.evs), "stds": float(data.stds)}
        else:                                  # Sampler: PUB 별 counts 도 남긴다(기록과 대조용)
            out[jid]["counts"] = [pub.join_data().get_counts() for pub in res]
    p.write_text(...); print("저장:", p)
```
Run it to a **new** file name; do not regenerate `20261111-063819Z_ibm_fez_ch14_isa_calib.json`.

### B11. HVA difference: "+0.111" (appendix) vs "+0.110" (ch14) — hard-coded rounded predictions in `appA_qpu_runs.py`
- `appA_qpu_runs.md:38`, `:84` vs `ch14_hardware.md:123`. Cause: `appA_qpu_runs.py:33-34` hard-codes `PRED = {0.918, 0.925, 0.759, 0.881}` and `PRED_HVA = -4.4329`; −4.3224 − (−4.4329) = 0.1105 prints as 0.111, while ch14 uses the unrounded prediction and prints 0.110. Nothing ties these constants to `ch14_hardware.py`; a change in the ch14 model would leave the appendix silently stale.
- Fix: let `ch14_hardware.py` write its predictions to a small JSON (e.g. `ch14/ch14_pred.json`: four success values, HVA with/without readout, ideal) and have appA read it; until then use `0.110` in both places. Other hard-coded constants that should come from a record or be asserted: `T_1Q, T_CZ` (assumed gate times — record real ones, B10), `"3"` at `:270` (use `anc[0]`), `GAMMA, BETA` (duplicated from `ch14_submit.hva_tfim`), the `2026-11-11 … QPU 6+3초` string in `ch14_submit.py:105`.

### B12. Bit-order statement is not true for QPE
- `appA_qpu_runs.md:30`: `비트열은 모두 책의 빅 엔디안(왼쪽이 q0)으로 적는다.`; `:141`: `다른 세 회로는 측정 매핑으로 처음부터 책 순서였다.`; `ch14_submit.py:108` labels every record `"bit_order": "book"`.
- In the submitted QPE circuit `QFTGate(t).inverse()` gives x = Σ q_j 2^j and the mapping is direct, so the printed string is c2c1c0 = q2q1q0: the left character is q2, the 4's place. (ch14's numpy version puts the MSB on q0 by dropping the SWAPs, so the *strings* agree, but the *qubit labels* do not.) The §3.5 reading ("000 = 1의 자리, 011 = 2의 자리") is correct as arithmetic. Anyone mapping a wrong QPE bit to a physical qubit with "left is q0" would repeat the GHZ mistake of §7.
- `:30` proposed: `비트열은 책의 빅 엔디안(왼쪽이 q0)으로 적는다. QPE만 예외로, 읽은 정수 x의 이진 표기다(왼쪽이 4의 자리이고, 제출한 회로에서는 q2다).`
- `:141` proposed: `BV와 Grover는 측정 매핑으로 처음부터 책 순서였고, QPE는 x의 이진 표기다.`
- `ch14_submit.py:108`: `"bit_order": "book"` → `"bit_order": "book", "bit_order_note": "qpe_t 는 x 의 이진 표기(c2c1c0 = q2q1q0)"`.

### B13. Cross-references "부록 A3" do not point to an appendix section
- `ch14_hardware.md:123` (`보정값 1.6\~2.5×10⁻³(부록 A3)`), `:178` (`…설명하지 못한다(부록 A3)`); report `§8.1:179`. "A3" is a label in the script's stdout; the appendix's sections are 1–10. → `(부록 A 3.4절)` at `:123`, `(부록 A 4절)` at `:178` and in the report. No stale sentence saying the QPU columns are empty or the script is unexecuted remains in ch14 (checked by diff against the pre-merge tree and by grep); section numbering 1–10 is intact.

### B14. ch14 §1 message is now contradicted by its own measured result
- `ch14_hardware.md:15`: `예측 안에 들어오면 노이즈 모델이 맞는 것이고, 벗어나면 사양서가 담지 않는 것(…)이 있다는 신호다.` The appendix shows GHZ landing on the prediction (+0.2 %p) for the wrong reason (§3.2, §4 item 1).
- Proposed: `예측 안에 들어왔다고 노이즈 모델이 맞았다는 보장은 없고(실측의 GHZ가 그랬다 — 부록 A 4절), 벗어나면 사양서가 담지 않는 것(크로스토크, 특정 큐비트의 나쁜 보정, 라우팅으로 늘어난 깊이)이 있다는 신호다.`

### B15. Single-run generalizations in ch14 §5
- `ch14_hardware.md:123`: `CZ가 더 많은 Grover가 예측보다 좋았으므로 그날 쓰인 연결은 중앙값(3×10⁻³)보다 좋았다 — 보정값 1.6\~2.5×10⁻³(부록 A3).` The inference runs backwards (the calibration values are the evidence, not the Grover result), the quoted range is job-wide (Grover's own two edges are 1.73×10⁻³ and 1.80×10⁻³), and better edges explain only 0.759 → 0.7785 (A3).
  - Proposed: `그날 이 작업이 쓴 연결의 CZ 오류는 1.6\~2.5×10⁻³으로 중앙값(3×10⁻³)보다 작았다. 다만 그 값으로 다시 계산해도 Grover는 0.779이고, 실측과의 나머지 차는 이 장 모형이 게이트 중 열완화를 두 번 세는 데서 온다(7절, 부록 A 3.4절).`
- Same line, last sentence: `사양서 숫자 네 개가 보지 못하는 "특정 큐비트 하나"가 얕은 회로의 결과를 정한다.` → `이 실행에서는 사양서 숫자 네 개가 보지 못하는 "특정 큐비트 하나"가 얕은 회로 둘(GHZ, BV)의 결과를 정했다.`
- `ch14_hardware.md:13`: `BV만 나쁘다(보조 큐비트 하나의 시작 오류로 설명된다 — 부록 A).` → `BV만 나쁘다(가장 큰 오답은 보조 큐비트 하나가 |1⟩에서 시작한 모양이다 — 부록 A 5절).` Also add the job ids here or `(작업 id는 5절)` — CLAUDE.md asks for device, date, job id and shots wherever a measured value is quoted; this line has three of the four.

### B16. Glossary: "깊이(게이트만)" is not the book's depth
- `appA_qpu_runs.md:159`: `| 깊이(게이트만) | 동시에 할 수 없는 게이트 층의 수. barrier·측정은 세지 않는다(`gate_cost` 규약) | 키트 전체가 이 규약으로 통일 |`. The kit's `gate_depth` shares only the barrier/measure exclusion with `tools/gate_cost.py`. It is taken on the ISA circuit in the native basis and counts virtual `rz` layers (Bell: depth 7 with rz 6, sx 3, cz 1), whereas ch00 §3 defines depth on CX + 1-qubit gates with all-to-all connectivity (Bell: 2). The two numbers are not comparable.
- Proposed: `| 깊이(게이트만) | 동시에 할 수 없는 게이트 층의 수. barrier·측정은 세지 않는다. 실행 기록의 깊이는 ISA 회로(rz·sx·x·cz, 가상 게이트 rz도 한 층)에서 센 값이라 본문의 깊이(CX와 1큐비트 게이트, 전결합 — 0장 3절)와 다른 수다 | Bell: 기록 7, 본문 규약으로는 2 |`

## 5. Grade C findings (polish)

- **C1 rounding ties.** `ch14_hardware.md:123` "GHZ 0.079" is 318/4000 = 0.0795 exactly; the appendix prints 0.0795 and plain rounding gives 0.080 (the script's 0.079 is a floating-point artefact; with 0.920 it does not even sum to 1). Print four decimals in block 7 as the appendix does. Likewise `0000` is written "5.15 %" in the appendix (`:64`) and "5.2 %" in `ch14_hardware.md:123`, `CLAUDE.md:57`, report `§8.1:178` (Python's own `f"{206/4000:.1%}"` gives 5.1 %). Use 5.15 % or 206/4000 everywhere.
- **C2 TVD column.** Add one sentence at `appA_qpu_runs.md:30`: `이 네 회로에서는 오답이 모두 이상값보다 많이 나와서 TVD가 (이상 성공 확률 − 실측 성공 확률)과 정확히 같다.`
- **C3 e without uncertainty.** `appA_qpu_runs.md:118`: `e = 0.0517에서 실측과 맞는다` → `e = 0.0517 ± 0.0036(0000 206회의 샷 잡음)에서 실측과 맞는다` (print `np.sqrt(p_obs*(1-p_obs)/n)/slope` in the script; slope dP/de = 0.984).
- **C4 "절반".** `appA_qpu_runs.md:135`, `:137`: the measured fraction is 50 ± 6 % (0.039 / 0.646). Write `50 ± 6%` and `그 이득의 **절반쯤**만 남았다`.
- **C5 precision vs stds.** `appA_qpu_runs.md:84` (`precision = 0.02는 목표였을 뿐이고 실제 오차는 그 두 배였다`), `:174`, `ch14_hardware.md:105`. It is not a missed target. From my replay distribution the single-shot variance of H is 2.43 (Z-basis part) + 1.87 (X-basis part); with 1/0.02² = 2500 shots per basis that gives 0.041, i.e. the reported 0.039 is what 2500–2800 shots per basis yield for this observable. Proposed `:84`: `precision = 0.02는 기저 하나에 쓸 샷 수(≈ 1/0.02² = 2500)를 정하는 값이다. H는 두 기저에서 잰 항 일곱 개의 합이라 한 샷의 분산이 1이 아니라 4쯤이고, 그래서 표준오차가 0.02가 아니라 0.039다.` Glossary `:174`: `Estimator의 목표 표준오차 · 실제 표준오차` → `Estimator에 주는 목표 정밀도(샷 수를 정한다) · 결과와 함께 돌아오는 표준오차 추정값`. (Print the two variances in appA before quoting "4쯤".)
- **C6 GHZ q4.** `appA_qpu_runs.md:54` last sentence `원인은 이 자료로 가를 수 없다.` → `q4는 사슬의 끝이라 q4의 시작 상태 오류, 마지막 CX 뒤의 X 오류, 판독 오류가 모두 같은 꼴을 낸다. 11110(74회)이 00001(47회)보다 많으므로 일부는 1 → 0 방향(이완이나 P(0\|1))이다. 이 자료로는 더 가르지 못한다.` (A symmetric initial flip of 2.6 % reproduces the total but predicts 64 : 57, observed 47 : 74, χ² 9.6.)
- **C7 "판독 오류로는 나올 수 없는 꼴"** (`appA_qpu_runs.md:64`) → `판독 오류로는 사실상 나올 수 없는 꼴이다` (three independent P(0|1) errors ≈ 6×10⁻⁷, not zero).
- **C8 glossary rows.** PUB (`:173`): `(회로 + 매개변수 + 샷/정밀도)` → `(회로 + [Estimator는 관측량] + 매개변수 값 + 샷 또는 정밀도)`. 판독 오류 (`:164`): add `IBM의 P(1\|0 준비)·P(0\|1 준비)는 준비 오류까지 섞인 값(SPAM)이다`. 초기화 오류 (`:165`): add `IBM은 직전 실험이 \|1⟩이었을 때의 값을 `init_error`로 따로 공개한다`. T1·T2 (`:163`): add `IBM의 T2는 Hahn 에코 값`. Clifford (`:170`): `BV는 Clifford라 트랜스파일러가 크게 다시 썼다` is an unverified causal claim (the rewrite seen in the ISA — H·H cancellation, X moved through a CZ — is ordinary commutation/resynthesis); → `BV는 Clifford라 파울리 결함 하나가 오답 하나로 정해진다(5절의 "확률 1")`. 클라이언트 쪽 프리미티브 (`:177`): `0.52 최상위 `Sampler`·`Estimator`` → `executor_sampler.Sampler`·`executor_estimator.Estimator`(설치된 0.52.0에서는 최상위 이름도 이것)`. The table has 30 rows (the task brief says 31 terms).
- **C9 CZ glossary row** (`:158`) "1큐비트 게이트보다 오류가 약 10배 크다" — that day 5–14× (qubit 24's sx is 1.0×10⁻³). Fine as "약 10배"; add "그날".
- **C10 ch14 time budget.** `ch14_hardware.md:140` predicts ≈ 5.4 s and says billed time "수십 초가 될 수 있지만"; the measured usage (6 + 3 = 9 s) is in the table but not here. Add `실측: Sampler 6초 + Estimator 3초(5절 표).`
- **C11 ch14 §7 line 174** `5절의 BV 오답이 "한 비트만 다른 이웃"에 모여 있는 것이 그 지문이다.` now reads as if about the measured BV, whose top wrong answer is 0000. → `5절의 BV 예측 히스토그램에서 오답이…`.
- **C12 ch14 line 170** "(하루 1회 보정)" vs the appendix's "37분 전" — after A1 both should say the same thing (daily benchmarking, calibrations as needed).
- **C13 ch00:13** defines "작성됨, 미실행" as "SDK를 설치할 수 없어 돌려 보지 못했다"; the ch14 snippet now carries that label for a different reason. Add `또는 발췌라 그대로는 돌지 않는다`.
- **C14 md escaping.** No `\uXXXX`, no unescaped `~ _ [` outside code, all table rows have consistent column counts and escape `|` as `\|` (mechanical check, both files). Bare `|` in running text (`|1⟩`, 11 lines in the appendix, e.g. `:11`, `:42`, `:122`, `:137`) violates the letter of CLAUDE.md's rule but matches the practice of every existing chapter including ch00 and ch14; either relax the rule's wording to "표 칸 안에서" or escape them all.
- **C15 `appA_qpu_runs.py` robustness.** `:226` `np.interp` clamps silently if `p_obs` leaves the grid → `assert f[0] < p_obs < f[-1]`. `:188-189` `_gap.seconds` drops whole days → `int(_gap.total_seconds())`. `:197` assumes one unmeasured qubit → `assert len(anc) == 1`. `:69`, `:184-185` reverse GHZ unconditionally → `assert "bit_order" not in real`. HVA value is read from the ISA file while the record has its own → `assert isa[EJ]["estimator"]["evs"] == real["jobs"]["hva"]["evs"]`; also `assert real["backend"] == isa[SJ]["backend"]`. (All in `r11_appA_patch.diff`, verified to pass.)
- **C16 `appA_qpu_runs.md:26`** `그 시각의 보정값` → `그 시각에 조회되는 보정값` (A1).
- **C17 Open plan wording.** IBM now states the allowance as "the RU equivalent of up to 10 minutes per 28-day rolling window" (1 RU = 1 minute of QPU time for standard workloads). "600초" is still right; no change needed unless the glossary wants to mention RU.
- **C18 report §8** mixes UTC ("06:21Z") and KST ("15:38") for the two runs; add the zone to the second.

## 6. Answers to the eight review questions, briefly

1. **Recomputation** — everything reproduces; one cross-document mismatch (B11) and tie-rounding cosmetics (C1).
2. **Replay (A3/A5)** — reproduced to four decimals with independent code (GHZ, BV, Grover, QPE, HVA; HVA also on the real ISA gate list). The approximations are only partly stated: missing are the gate-time double count (A3), echo-T2, the calibration's own sampling error and age (A1, B2). "실측−재현 −8.1σ" appears only in stdout, but 5.7σ / 2.3σ / 4.4σ are in the text with no caveat (B2).
3. **BV hypothesis (A4)** — enumeration reproduced; the mechanism "ancilla is |1⟩ at the first CZ" is well supported (and gains support from HVA, B4, and from the full-histogram fit, A2). Not supported: "only single-fault mechanism" (A6), "exceeds the calibration's upper bound" (A2), "within those 37 minutes" (A1). Mechanisms the enumeration cannot see and how they fare: two independent faults and triple readout flips — negligible (≈10⁻⁶); correlated `cz[4,3]` fault — real but 0.06 % at the calibrated rate; ancilla X failure / over-rotation — gives 1000, observed 0.27 %; no-pulse shots — give 0000, argued against by QPE/GHZ in the same job; leakage and reset/carry-over — open, and carry-over is the natural reading given the ancilla ends each shot 50 % excited and unmeasured; crosstalk in the H layer — cannot produce a deterministic triple flip; wrong layout — excluded by the noise-free check on the recorded ISA; stale calibration — likely in general (A1) but irrelevant to the mechanism. The proposed experiment does not discriminate (B5); a sharper one is given.
4. **Statistics** — formulas and N correct; σ against models needs the caveat (B2); generalizations flagged in A7, B15, A3.
5. **Glossary** — mostly correct and consistent with ch00/ch14 (r → p conversion, TVD, ISA, heavy-hex, TREX at resilience level 1, DD/twirling defaults all check out against IBM docs). Corrections: calibration row (A1), single-fault row (A6), depth row (B16), σ row (B2), and the C8 list. EPLG is not used anywhere in the repo.
6. **ch14** — still reads as a prediction chapter; measured values are labelled with device/date/shots (job ids missing only at `:13`). New columns agree with the records and with the appendix. Wrong or overstated: `:178` (A4, A1), `:123` (A5, B15, B13), `:15` (B14), `:77`/`:80` (B7), `:107` (B6).
7. **Scripts** — B8, B9, B10, B11, C15. No unintended QPU spend; no secret leakage found.
8. **Rules** — numbers come from scripts except hand-rounded "5.2 %" and fake values in `ch14_hardware.md:178` (taken from appA's four-decimal output); "실측" is used only for QPU values in the reviewed files; escaping per C14; html for both md files is current (same strings present).

## 7. What I could not verify

- Anything that needs the SDK or the IBM account: the contents of qiskit-ibm-runtime 0.52.0 (`__init__.py`, `server_side`, `finalize_options.py`, `options_models/`), qiskit-aer 0.18.2's `estimator_v2.py` (the N(expectation, precision) claim) and Aer's `from_backend` construction (A3 relies on my recollection there); that `init_error` is actually populated for ibm_fez; the true gate durations on ibm_fez (the replay assumes 32/68 ns); the options the two jobs actually ran with (`resilience_level`, twirling, DD, `rep_delay`, `init_qubits` — none is in the frozen records; IBM's docs do confirm the *defaults* claimed: Sampler DD/twirling off, Estimator `resilience_level` default 1 and "turn off all error mitigation and suppression" at 0); how IBM initializes qubits between shots on Heron (the docs say only that "an implicit reset is typically inserted at the beginning of the circuit" and that reset imperfections cause state-preparation errors).
- IBM's release notes for 0.52 (the page ends at 0.51.0 today) — see B6.
- That calibration readout values are 1024-shot estimates: inferred from every value being a multiple of 1/1024, not documented.
- That ibm_kingston / ibm_fez are "Heron r2": not in the records; IBM's processor-types page confirms Heron r2 = 156 qubits but does not list backends.
- My replacement code for `ch14_submit.py` and `fetch_ibm_job.py` is untested. The appA patch and all numbers quoted from it were run.
- The carry-over explanation (A2) is a hypothesis that fits the existing counts; it is not established. For data qubits it is not separable from extra readout error with this data.

## 8. Scratch files

`/tmp/claude-0/-home-claude/fd63146d-511d-5ae9-8dab-3100fb85f0db/scratchpad/mid/R11/`: `r11_counts.py`, `r11_replay.py`, `r11_more.py`, `r11_hva.py`, `r11_hva_isa.py`, `r11_carryover.py` (independent code), `appA_qpu_runs_patched.py` + `r11_appA_patch.diff` + `appA_patched_out.txt` (proposed additions applied to a copy of the author's script and run), `appA_rerun.txt`, `ch14_rerun.txt`, `ch14_before.txt` (outputs), `copy/`, `before/` (scratch copies used to run scripts without writing into the repo).

## Sources

- https://quantum.cloud.ibm.com/docs/en/api/qiskit-ibm-runtime/release-notes (first 100,000 of 125,709 characters read; covers 0.51.0 and 0.50.0)
- https://quantum.cloud.ibm.com/docs/guides/qpu-information (prob_meas1_prep0, prob_meas0_prep1, readout error, init_error, T2 Hahn echo, EPLG, CZ/SX error definitions)
- https://quantum.cloud.ibm.com/docs/guides/calibration-jobs (per-property timestamps, last_update_date, hourly monitoring, daily benchmarking)
- https://quantum.cloud.ibm.com/docs/en/guides/repetition-rate-execution (rep_delay, implicit reset, init_qubits)
- https://quantum.cloud.ibm.com/docs/en/guides/estimator-options (resilience_level default 1, level 0, default_precision, twirling/DD defaults, init_qubits, rep_delay; legacy vs client-side import comments)
- https://quantum.cloud.ibm.com/docs/en/guides/sampler-options (DD/twirling defaults False; legacy vs client-side import comments)
- https://quantum.cloud.ibm.com/docs/en/guides/configure-error-mitigation (TREX/ZNE and resilience level; link to server-side → client-side migration)
- https://quantum.cloud.ibm.com/docs/en/guides/migrate-from-sampler-to-executor (Executor has no built-in suppression/mitigation)
- https://quantum.cloud.ibm.com/docs/en/api/qiskit-ibm-runtime/models-backend-properties (qubit_property returns value and measurement time)
- https://quantum.cloud.ibm.com/docs/en/api/qiskit-ibm-runtime/dev/options (options classes)
- https://quantum.cloud.ibm.com/docs/guides/plans-overview (Open plan: RU equivalent of 10 minutes per 28-day rolling window)
- https://quantum.cloud.ibm.com/docs/en/guides/processor-types (Heron r1/r2/r3)
- https://qiskit.github.io/qiskit-aer/stubs/qiskit_aer.noise.NoiseModel.html (from_backend: depolarizing followed by thermal relaxation)
- https://docs.quantum.ibm.com/api/qiskit/0.39/qiskit_aer.noise.device.basic_device_gate_errors (parameter descriptions only)
- Not retrievable in this session: the qiskit-aer source for `noise/device/models.py` (two mirrors refused or returned no content).
