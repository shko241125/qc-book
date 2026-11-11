# R4 findings — plan compliance and whole-book coherence (ch01–ch17)

Reviewer R4. Scope: (1) compliance of the 17 chapters with the plan document (`plan.txt` §2, §4, §5, §6, §8),
(2) read-through for cross-chapter coherence, (3) what ch18 must deliver, (4) five structural weaknesses.
Nothing under `/home/claude/qc-book/` was edited. Scratch code: `mid/R4/ff_check.py`, `mid/R4/syn.py`.

How I read: ch01, ch12, ch14, ch15, ch17 in full; ch13 and ch16 sections 1, 3, 5 (partly), 7, 8, 9; for the
other chapters sections 2, 6, 7 and the header/goal lines, plus book-wide greps for every cross-reference,
symbol and repeated number listed in the task. Per-chapter numeric agreement with the scripts was not re-checked
(out of scope per the brief).

---

## Part 1 — Compliance matrix

Legend: P = pass, p = partial, F = fail, – = not applicable. Columns follow the task letters.

- **a** goal of plan §4 delivered · **b** template / 교과서 비교 / Quirk · **c** §6 checklist (sources · one other SDK + reason · §7 classical simulability · 3 exercises · 학습 노트)
- **d** "2-qubit gates / depth / T-count" (decision 2026-10-27, part 2 only) · **e** code rules · **f** 교정 규칙 · **g** simulator AND QPU, honesty (part 2 only)

| ch | a goal | b template | c checklist | d 3 numbers | e code rules | f 교정 | g hardware |
|---|---|---|---|---|---|---|---|
| 01 | P | P | P (QuTiP, reason given) | – | p (no pin in first cell) | p (no English gloss for 중첩/측정/간섭) | – |
| 02 | P | P | P (Cirq, reason) | – | p (same) | p | – |
| 03 | P | P | P (Cirq, reason) | – | p (same) | p | – |
| 04 | p (no transpile output exists to "read"; §5 table is hand-derived; "(c)(d) 결과는 실행 후 수치로") | P | p (§7 has no classical-simulability item) | – | p (`target = ...` l.70) | p | – |
| 05 | P (numpy); Qiskit `if_test` unexecuted | P | p (§7 no classical-simulability item) | – | p | p | – |
| 06 | P (numpy) | P | P (PennyLane, reason) | – | **F** (first cell does not parse, l.83) | p | – |
| 07 | p (goal reworded "설명"; plan's main tool QuTiP unexecuted, numpy RK4 substitutes) | P (reason for limited Quirk given) | **F**/p (5 SDKs, no single choice; §7 no classical-sim item) | – | p (`...` l.84; unseeded Aer/Cirq/CUDA-Q sampling) | P | – |
| 08 | P | P | p (2 other SDKs, no reason) | **F** (no numbers at all) | p (`...; return` l.98; `StatevectorSampler()` unseeded) | p | p (phase-damping sim only; "14장에서 실제 장비에 올릴 때") |
| 09 | P | P | p (3 other SDKs) | **P** (n=3: 2q 12, T 14, depth ≈30 per iteration) | p (`...` l.85; unseeded) | p (analogy without break line, l.262) | p ("14장에서 실제 올리는 회로") |
| 10 | P | P | p (2 other SDKs) | p (2q ≈15, depth ≈12; no T-count; only t=4, U=T) | **F** (first cell does not parse, l.94) | p | p ("14장에서 시도") |
| 11 | p ("회로를 끝까지" = permutation applied to a state array; no gate-level modular arithmetic) | p (no Quirk link; reason given: URL too long) | p (2 other SDKs; one B-grade news source, labelled) | **F** (asymptotics and literature numbers only) | p (unseeded) | p | p (honest: "14장은 Shor 회로를 실행하지 않는다"; no noise sim) |
| 12 | p (`PauliEvolutionGate` snippet present, unexecuted) | P | **F**/p (4 other SDKs; §10 has no links and no grade labels) | p (CNOT, depth; no T-count) | **F** (first cell: body-less `def`s) | p (t used for time and for counting bits) | **F** (no noise sim, no QPU, not in ch14's list, and no statement about it) |
| 13 | P (numpy); PennyLane unexecuted | P | **F**/p (Cirq + CUDA-Q snippets, 5-column table; §10 no links/grades) | p (2q only; no depth for HVA/QAOA, no T) | **F** (pseudo-code `psi = |0...0>`) | p | p (HVA p=1 predicted in ch14; QAOA not; "14장은 … 올리므로") |
| 14 | **F** vs plan ("실제 QPU에서 돌리고" → prediction only; header goal rewritten, no §8 record) | P | p (2 snippets + 5-column table; checklist item "QPU 결과와 사용 분" empty) | p (CZ/SWAP/1q; no depth, no T) | **F** (body-less `def`s; `>=0.40` is not a pin) | p ("이 세션", "이전 판" narration ×10) | P on honesty (clear, consistent, columns "—") |
| 15 | P | P | **F**/p (4 SDK columns + 3 snippets; §10 no links/grades) | – | **F** (body-less `def`s) | p (analogy without break line, l.253) | – |
| 16 | P (all five primitives, numpy matrices) | P | **F**/p (Qualtran + pyqsp, unpinned; §10 no links/grades) | – | **F** (pseudo-code with `e^{iφ₀}`, U+2212) | p (two analogies without break line, §2) | – |
| 17 | p (checklist run on 8 candidates; plan's main tool `pennylane.estimator` unexecuted) | P | p (Qualtran snippet, no "why" sentence, unpinned) | – | p (excerpt, no imports) | p (analogy without break line, l.249) | – |

All 17 chapters: three exercises with answers directly below each (P), 학습 노트 present with 5–6 bullets (P),
교과서 비교 box at the end of section 3 comparing N&C and KLM (P).

Evidence for column e: `mid/R4/syn.py` runs `ast.parse` on all 55 ```python blocks. All 38 SDK blocks parse; 7 of
the 17 first (numpy) cells do not: ch06 l.83, ch10 l.94, ch12 l.66, ch13 l.58, ch14 l.56, ch15 l.74, ch16 l.108.
Five more parse but are excerpts with no `import` (ch07, 08, 09, 11, 17).

---

## Findings (sorted by severity)

### A — factually wrong

**[A1] ch12 §3 — `ch12/ch12_hamiltonian.md:32`**
Current: "ZZ끼리는 교환, X끼리도 교환, ZZ와 X는 비교환 — 그래서 양자 상전이가 있고 고전적으로 어렵다."
Wrong: the 1D transverse-field Ising chain is the textbook exactly solvable model (Jordan–Wigner → free fermions,
Pfeuty, *Ann. Phys.* 57, 79 (1970)). It is classically easy for statics and for dynamics. The book never says so:
"Jordan–Wigner / 자유 페르미온" appears once in 17 chapters (ch15:57, as a generic boundary), yet this model is the
running example of ch12, 13, 14, 15, 16 and 17.
Evidence (`mid/R4/ff_check.py`): Heisenberg evolution of 2n Majorana operators with a 2n×2n matrix gives
⟨Z₀⟩(t=6) = 0.178056 for n = 4 — the "정확값 0.178" of ch12 §5 — and agrees with the state vector to 8 digits for
n = 8, 10, 12; n = 1000, t = 6 takes 14 s.
Replacement:
> n = 4에서 항 7개. ZZ끼리는 교환, X끼리도 교환, ZZ와 X는 비교환 — 그래서 양자 상전이가 있다. 단 1차원 TFIM은 Jordan–Wigner 변환으로 자유 페르미온이 되어 고전으로도 다항 시간에 정확히 풀린다(Pfeuty 1970; 15장 3절의 "자유 페르미온" 경계). 이 책이 이 모형을 쓰는 이유는 어려워서가 아니라 정답을 알아 검산할 수 있어서다. 고전으로 어려워지는 것은 세로 자기장 −g Σ Zᵢ를 더하거나 2차원으로 갈 때다.

### B — misleading, incomplete, or contradicts another chapter

**[B1] ch17 §1, §5(A) — `ch17/ch17_discovery.md:11`, `:147`, `:180`** (consequence of A1)
Current: "χ는 t = 1에서 4, t = 6에서 180으로 지수적으로 는다 — 실시간 동역학은 고전이 어려운 쪽이다" and
"후보 표 다섯째 줄(스핀 격자 급랭 동역학, 23a장)이 살아남는 근거다."
The MPS measurement is correct, but the model is free-fermion: the same quench is exact with 2n×2n matrices (see
A1). The experiment shows the cost of one classical method, on a model that trips another boundary listed in
ch15 §3. As written it is the book's only numerical evidence for "dynamics is the hard side" and for keeping ch23a.
Replacement (append to l.147), and re-run (A)-dynamics once with a longitudinal field added (a one-line change in
`ch17_discovery.py`), so the evidence is on a non-integrable model:
> 주의: 1D TFIM은 자유 페르미온 모형이라 이 급랭도 2n×2n 행렬로 정확히 계산된다(n = 1000, t = 6의 ⟨Z₀⟩가 수 초). 위 표가 잰 것은 "MPS라는 한 가지 고전 방법의 비용"이다. "고전으로 어렵다"의 근거로 쓰려면 세로 자기장을 더한 비적분 모형에서 같은 표를 다시 재야 하고, 후보 표의 판정도 그 모형 기준으로 읽는다.

**[B2] ch17 §3, §4, §5 — `ch17/ch17_discovery.md:59`, `:101`, `:159`, `:179`** vs `ch10/ch10_qpe.md:255`, `ch12/ch12_hamiltonian.md:231`
ch17: "QPE가 바닥상태를 뽑을 확률이 γ²이므로 반복 횟수는 1/γ²(진폭 증폭으로 1/γ)".
ch10: "겹침이 γ = |c_0|²인 시험 상태를 넣으면, 한 번 실행당 확률 γ … 1/γ번 반복"; ch12: "γ = |⟨ψ₀|φ⟩|² … 확률 γ", and
ch12 §5 uses 1/γ ≈ 4 for γ = 0.23. Same symbol, same name ("바닥상태 겹침 γ"), amplitude in ch17 and probability in
ch10/12; a reader carrying γ = 0.23 from ch12 into ch17 gets 19 repetitions instead of 4.
Replacement for l.59 (and "1/γ²" → "1/γ" at l.101, 159, 179):
> 그리고 **초기 상태의 바닥상태 겹침 γ = |⟨ψ₀|φ⟩|²**(10·12장과 같은 정의)가 따로 곱해진다 — QPE가 바닥상태를 뽑을 확률이 γ이므로 반복 횟수는 1/γ(진폭 증폭으로 1/√γ)다.

**[B3] Part 3 has no place for variational and kernel methods, on which three of five part-4 cases rest**
- `ch15/ch15_sources.md:5` promises to classify "2부의 알고리즘"; the table (l.35–42) has no row for ch13, and the
  words 변분/VQE/QAOA/커널 do not occur anywhere in ch15.
- `ch16/ch16_primitives.md:30`: the dictionary's premises are "8~12장의 부품"; ch13 is skipped.
- `ch17/ch17_discovery.md:43–51`: no structure→part row; l.227 says the checklist "로는 가능성을 판정할 수 없고".
- Plan §4 candidate table calls "변분" and "커널" 가속 원천 and says the five selected cases cover "가속 원천 네 가지
  (해밀토니안 시뮬레이션·변분·커널·블록 인코딩)"; `ch13/ch13_variational.md:13` says VQE/QAOA are "양자 가속 알고리즘이 아니라 … 휴리스틱".
- Plan ch19: "17장 절차로 후보 10개를 평가" — for candidates 2, 3, 5, 6, 9 the procedure of ch17 returns "cannot judge".
Fix: one row in ch15's table, one in ch17's, and the heuristic track in ch18 (see Part 2, item list).
Row for the ch15 table:
> | 변분·커널(휴리스틱) | 얕은 파라미터 회로를 고전 최적화기로 조정, 또는 특성 맵의 내적 | 6장 feature map, 13장 VQE·QAOA | 증명된 가속 없음 | 훈련 가능하다고 증명된 회로는 대개 고전 시뮬레이션 가능(Cerezo 2023), 커널 집중 |

**[B4] The same Babbush 2021 result is quoted three incompatible ways — `ch09:232`, `ch16:321`, `ch17:11/:170`**
- ch09: "고전 병렬도 10³에서 100일(하한 모델)~880년(SA 모델), 10⁶에서 280년~88만 년 (Table I)" — matches the paper.
- ch16: "제곱 가속만으로는 수일\~수주의 실행에서도 이기기 어렵고" — not a figure from the paper.
- ch17: own calculator, "100초\~32년", presented as "Babbush et al. 2021의 논리", with logical Toffoli τ = 1–10 µs.
Checked against the paper (ar5iv 2011.04149, Table I and text): Toffoli time ≈ 170 µs; quadratic break-even
2.4 h / 320 d (S = 1), 100 d / 880 yr (S = 10³), 280 yr / 880 millennia (S = 10⁶); cubic 12 s … 40 d; quartic
1.4 s … 4.9 h. ch17's formula reproduces the paper's 2.4 h, so only the inputs differ: τ is 17–170× more optimistic
than the source it cites, and also than the ≤ 93 µs per Toffoli implied by the Gidney 2025 figures quoted in the
same chapter (6.5×10⁹ Toffoli in under a week). With τ = 170 µs the four p = 2 scenarios of ch17 §5(C) become
33 days … 9,200 years (T* ∝ τ²).
Replacement for ch16:321:
> Babbush et al.(2021)은 논리 Toffoli 하나에 약 170 µs가 드는 표면 코드 가정에서, 제곱 가속의 손익분기 실행 시간이 단일 코어 상대 2.4시간(하한 모형)\~320일(SA 모형), 고전 병렬도 10³ 상대 100일\~880년임을 보였다(Table I, 9장 7절과 같은 표). 3차는 12초\~40일, 4차는 1.4초\~4.9시간이다.

Sentence to add after the ch17 §5(C) table (l.170):
> 이 표의 τ = 1\~10 µs는 낙관적인 가정이다. Babbush et al.(2021)이 쓴 논리 Toffoli 170 µs를 넣으면 같은 네 시나리오의 p = 2는 33일\~9,200년이 된다 — 9장 7절이 인용한 100일\~88만 년과 같은 자릿수다.

**[B5] ch14 does not deliver the plan's goal, and ten sentences elsewhere still say it runs on hardware**
Plan §4 ch14: "8~13장을 실제 QPU에서 돌리고 노이즈를 해석한다"; principle 4: every algorithm shows simulator and real-QPU
results. `ch14/ch14_hardware.md:5` rewrites the goal to prediction, and l.9, 110, 184, 218 say honestly that nothing
was submitted — ch14 itself is consistent and nothing in it reads as a real result. But §8 has no decision
recording the change, and other chapters still describe ch14 as a hardware run:

| where | current | replacement |
|---|---|---|
| ch04:11 | "14장에서 실제 QPU 결과를 해석할 때 이 숫자들이 노이즈의 예측치가 된다." | 14장에서 노이즈 모델로 결과를 예측할 때(그리고 나중에 실제 QPU 결과를 해석할 때) 이 숫자들이 입력이 된다. |
| ch07:60 | "(14장에서 실측치를 읽는다)" | (14장에서 공개 사양서의 값을 읽는다) |
| ch07:221 | "14장의 하드웨어 실행 전 리허설은 항상 Qiskit Aer의 `NoiseModel.from_backend(backend)`로 한다" | 14장의 리허설은 공개 중앙값으로 만든 numpy 노이즈 모델로 하고, 계정이 열리면 `NoiseModel.from_backend(backend)`(제출 스크립트의 fake 모드)로 한 번 더 한다 |
| ch08:181 | "14장에서 BV를 실제 장비에 올릴 때 이 표가 기대치다." | 14장에서 BV를 Heron r2 노이즈 모델로 예측할 때 이 표가 비교 기준이다(실기 제출은 아직 못 했다). |
| ch09:90 | "14장에서 실제 올리는 회로는 n = 3, k = 2 버전이다." | 14장에서 노이즈 모델로 예측하고 제출 스크립트에 넣어 둔 회로는 n = 3, k = 2 버전이다. |
| ch09:236 | "14장은 k = 2만 올린다" | 14장은 k = 2만 예측한다 |
| ch10:100 | "…SWAP을 없앤 판(…)을 올린다." | …SWAP을 없앤 판(…)을 노이즈 모델로 예측한다(실기 제출은 미실행). |
| ch10:241 | "(14장에서 시도)" | (14장에서 노이즈 모델로 예측, 실기 미실행) |
| ch13:270 / :278 | "그래서 14장의 실기 실행 후보다." / "14장은 HVA p = 1만 올리므로" | 그래서 14장의 예측·제출 목록에 들어간다(실기 실행은 아직 못 했다). / 14장은 HVA p = 1만 예측하므로 |

Also not covered by ch14's five circuits and not stated anywhere: ch12 (Trotter dynamics, energy QPE, amplitude
estimation — no noise simulation in ch12 either), ch13 QAOA, ch08 DJ/Simon. ch11 states its exclusion explicitly.
Proposed §8 line (needs the author's approval, because it touches invariant principle 4):
> 2026-11-xx | 14장 목표를 "실제 QPU 실행"에서 "사양서 기반 예측 + 제출 스크립트"로 낮춘다. 원칙 4는 유지하고, 2부 각 알고리즘의 실기 결과는 계정·네트워크가 열리면 14장 5절 표의 빈 열로 채운다. 12장 Trotter 회로와 13장 QAOA를 14장 예측 목록에 추가한다 | QPU 접근 불가(계정·네트워크 없음). 8~13장의 "14장에서 올린다" 문장을 "예측한다"로 통일

**[B6] Decision 2026-10-27 ("2큐비트 게이트 수 / 깊이 / T-count 세 숫자를 함께") is met by one circuit in part 2**
- Complies: ch09:90, Grover n = 3 (2q 12, T 14, depth ≈ 30 per iteration; total 2q 24, T 28).
- Two of three: ch10:100 QPE t = 4 (no T-count); ch12:88 Trotter (CNOT, depth; no T-count); ch13:84 HEA (2q, depth),
  HVA and QAOA (2q only); ch14:112–118 all five circuits (CZ, SWAP, 1q; no depth, no T-count).
- None: ch08 (DJ, BV, Simon); ch10 IPE and φ = 1/3 circuits; ch11 (N = 15, 21: "수천 게이트", asymptotics); ch12 energy
  QPE ("수만 게이트") and amplitude estimation; ch09 n = 2 exercise and general amplitude amplification.
The rule also has a gap: for circuits with arbitrary-angle rotations (ch10, 12, 13) a T-count exists only relative
to a synthesis precision. Proposed amendment to the decision, and a paragraph for ch08 (numbers are hand-counted
for the chapter's own n = 6 example and should be confirmed by a counter script):
> 2026-11-xx | 임의 각도 회전이 있는 회로의 T-count는 "회전 수 R"과 "T ≈ 3 log₂(1/ε)·R(Ross–Selinger)"로 적는다 | 회전의 T-count는 합성 정밀도 없이는 정해지지 않음

> **게이트 비용 기록(이 책의 규칙).** BV(n = 6, s = 101101): 2큐비트 게이트 4(s의 1의 개수), 깊이 ≈ 7(보조 준비 2 + 보조를 공유하는 CNOT 4 + 마지막 H 1), T-count 0(Clifford 회로). DJ의 균형 오라클도 CNOT 사슬이라 같은 꼴이다. Simon은 이 장에서 오라클을 치환 행렬로 두었으므로 세 숫자가 없다 — 선형 오라클로 짜면 CNOT O(n²), T-count 0이다.

**[B7] Code rule "복사해 그대로 실행되는 완결 스니펫, 생략 표시(...) 금지" fails from ch06 onward**
First (numpy, "실행 완료") cell: SyntaxError in ch06:83 (`inv = [ ...QFT(b) 와 같은 리스트... ]`), ch10:94, ch12:66
(four body-less `def`), ch13:58 (`psi = |0...0>`), ch14:56, ch15:74, ch16:108 (`W = diag(2Π − 1) @ U`). Literal `...`
in code also at ch04:70, ch07:84, ch08:98, ch09:85, ch13:72/76/79, ch14:64, ch16:105. ch01–05 comply.
The label "numpy 구현(실행 완료)" above these cells is true of the `.py` file, not of the cell the reader sees.
Fix, cheapest form — one sentence before each such cell (or fence them as ```text):
> 아래는 `chNN_<slug>.py`에서 핵심만 간추린 발췌다(그대로 실행되지 않는다). 실행 가능한 전체 코드는 스크립트 파일에 있다.

**[B8] Principle 3 / checklist "다른 SDK 병기는 하나뿐이고, 왜 그 SDK인지 한 문장" holds only in ch01–06**
ch01–06 each open section 6 with "이 장에서 X를 고른 이유는 …". From ch07 on, section 6 is a 3–5 column table plus
several snippets and no stated choice: ch07 (QuTiP + Aer, Cirq, PennyLane, CUDA-Q), ch08 (Cirq, PennyLane), ch09
(PennyLane, Cirq, CUDA-Q), ch10 (PennyLane, Cirq), ch11 (Cirq, PennyLane), ch12 (QuTiP, PennyLane, Cirq, CUDA-Q),
ch13 (Cirq, CUDA-Q + QuTiP column), ch14 (Cirq, PennyLane + CUDA-Q, QuTiP columns), ch15 (Cirq, PennyLane, CUDA-Q),
ch16 (Qualtran, pyqsp + CUDA-Q, QuTiP columns), ch17 (Qualtran + CUDA-Q, Azure QRE columns). Qualtran, pyqsp, Azure
QRE, Stim and quimb are not in plan §3 at all. 34 of the 51 unexecuted snippets are non-Qiskit; one per chapter would be 17.
Fix per chapter: keep the table, keep one snippet, and state the choice, e.g. for ch12:
> 이 장에서 QuTiP을 고른 이유는 Trotter 없이 Schrödinger 방정식을 적분하는 독립 검산기이기 때문이다. 다른 SDK는 아래 표의 대응만 적고 코드는 싣지 않는다.

**[B9] Promises made to part 3 that ch15–17 do not keep (each must be kept in ch18 or withdrawn)**
- ch06:211 "3부의 발견 절차에서 '변분 회로를 쓰려면 불모 고원을 먼저 검사한다'가 체크리스트 항목이 된다." — "불모 고원" does not occur in ch15–17.
- ch15:110 (three Aer methods) "17장에서 자동 검사 도구로 쓴다" — ch17 uses a numpy entropy check only.
- ch08:226 "Cirq의 `_decompose_` 패턴이 이 책 3부에서 가장 유용하다" — `_decompose_` does not occur in ch15–17.
- ch11:227 "3부의 '알려진 알고리즘도 다시 볼 수 있다'의 예" — no such theme in ch15–17 (Regev is not mentioned).
- ch15:59 "17장 체크리스트의 두 번째 문항 '이 후보는 어느 고전 경계에 걸리는가'" — in ch17 step 2 is 구조 찾기; the classical boundary is step 5. Replacement:
> 17장 체크리스트의 5단계(고전 베이스라인)가 묻는 "이 후보는 어느 고전 경계에 걸리는가"가 여기서 나온다.
- Plan §4 candidate 10: "옵션 가격 결정 → 16장 진폭 추정 예제" — ch16's AE example is a random 3-qubit A; option pricing appears only as a literature row in ch17.
- ch07:254 promises "(회로의 총 시간) × (1/T1 + 1/Tφ)/3" as ch14's go/no-go formula; ch14's budget is gate-count based and says T1/T2 are negligible for its circuits (l.108).

**[B10] Symbols that change meaning between chapters that build on each other**
- Counting-register size: t (ch10, ch11), t and t_b (ch12, where t is also time: l.48 "ΔE = 2π/(τ·2^t) … t ≈ 20비트" next to "시간 t"), m with M = 2ᵐ (ch16), b (ch17).
- M: number of solutions (ch09:11 "해가 M개면"), number of Q calls 2^t (ch12:50, ch16), confusion matrix (ch14). ch09:234 sends the reader to ch12 for "M을 먼저 재는 양자 계수 세기", where M already means something else.
- γ: see B2. ε: precision, block-encoding error and marked fraction, all inside ch16 §3. κ: sign-threshold steepness, undefined, at ch16:228; condition number at ch16:315 and ch17:183. Σ|cₖ|: w in ch12:48, α in ch16/17.
- 충실도: average gate fidelity (ch07, ch14), state fidelity (ch05, ch13, ch15). In the ch12 Trotter table (l.146–153) the columns "1차/2차 충실도" are state fidelities from |0000⟩ (script l.60) but the text never says so, beside a column that is an operator-norm error.
Fix: one notation page for the book (weakness 5) and, for ch12:144:
> **Trotter 오차 vs 스텝 수(t = 2).** ‖ΔU‖는 연산자 노름 오차, 충실도는 초기 상태 |0000⟩에서의 상태 충실도 |⟨ψ_정확|ψ_Trotter⟩|²이다.

**[B11] ch15 — Forrelation is "1 query" in the table and notes, "two queries" in the experiment — `:42`, `:254` vs `:11`, `:181`**
The circuit of §5(5) queries f once and g once. Aaronson–Ambainis's one-query algorithm uses a control qubit that
ch15 never shows (recalled from the paper, not re-fetched today). Replacement for the table cell (l.42) and the note (l.254):
> 쿼리 2(제어 큐비트를 쓰면 1) vs Ω̃(√N), 최대 분리

**[B12] ch14 §2, §5, §7 — empirical statements about real devices in a chapter with no run and no citation — `:19`, `:133`, `:176`, `:184`**
"최악 연결은 중앙값의 10배", "특정 CZ 연결이 중앙값의 4배(1.2%)인 것은 흔하므로", "일주일 간격으로 내면 성공 확률이 ±5%p 움직이는 것이 정상", "실측이 두 예측보다 나쁜 것이 보통". Unverified by me; they may be true, but principle 6 asks for a primary source and the
chapter has no data of its own. Either cite a backend calibration table or mark them:
> (공개 보정 데이터에서 흔히 보이는 범위이며, 이 책이 직접 잰 값은 아니다.)

### C — minor

**[C1] Seeds not fixed in SDK snippets (code rule).** `StatevectorSampler()` without seed: ch08:133, ch09:110, ch10:117/126, ch11:116, ch16:170 (ch01, 03, 05 use `seed=7`). Unseeded `cirq.Simulator()`/`DensityMatrixSimulator()`: ch07:208, ch08:211, ch09:210, ch10:228, ch11:200, ch14:155. Unseeded Aer runs: ch07:200, ch15:96/104. `cudaq.sample` without `cudaq.set_random_seed`: ch07:218, ch09:221.

**[C2] Version pins.** No chapter pins anything in its first code cell (numpy/scipy versions are only in section 10). `cuda-q==0.16` (ch07, 09, 12, 13, 15) is not an installable spec — the PyPI package is `cudaq`; write `cudaq==0.16.x`. Plan §3 still lists CUDA-Q 0.15.0 (0.16.0 exists per NVIDIA's NGC catalogue; no §8 line records the bump). `qiskit-ibm-runtime>=0.40` (ch14:78) is a floor, not a pin. `qualtran`, `pyqsp` unpinned (ch16:297/306, ch17:203). `pennylane==0.45` (ch17:108) vs `0.45.x` elsewhere.

**[C3] Source grading stops at ch11.** ch01–11 label every source "등급 A/B"; ch12–17 have no labels. Section 10 of ch12, 13, 15, 16 contains no hyperlink at all (checklist: "정의·수치마다 A등급 출처 링크"). ch15:271 lists the Shtetl-Optimized blog without the B label that principle 6 requires.

**[C4] Exercise mix.** Plan §3: "두 책에서 1개씩 변형 + 자작 1개". Only ch02 does that. ch01, 08, 12: one textbook + two 자작; ch14, ch17: three 자작; ch13, 15, 16: papers instead of textbooks (ch13 explains why). Either relax the rule in §8 for chapters the textbooks do not cover, or label the exception in each chapter.

**[C5] Term rule "첫 등장에서 한국어(영어)".** No English gloss anywhere in the book for: 중첩, 얽힘, 간섭, 결어긋남, 텐서망, 탈분극, 판독 오류, 안정자, 충실도, 밀도행렬, 결합 차원, 급랭, 연분수, 손익분기, 진폭 증폭, 진폭 추정, 블록 인코딩 (grep for "(superposition", "(entanglement", "(decoherence", "(tensor network", "(depolariz", "(readout", "(stabilizer", "(density matrix", "(bond dimension", "(quench", "(amplitude a…", "(block encoding" returns nothing in prose). The rule's second half ("이후 영어 용어를 쓴다") is also not what the book does — it keeps the Korean term, per the 6차 review's 용어 통일. Related: "화학 정밀도" (ch12:48, ch13:192/276/302) vs "화학적 정확도" (ch17:241/243) vs plan "화학 정확도"; "orthogonality catastrophe" (ch10, ch12) vs "직교 파국" (ch17:225).

**[C6] Analogies without a "깨지는 곳" line (sample).** ch09:262 ("2배 빠른 차가 시동에 한 시간"), ch13:299 (자동차 평행주차), ch15:253 ("긴 사슬 … 실"), ch16 §2 ("블록 인코딩 = 큰 유니터리의 한 구석", "양자 워크 = 확산 대신 직진"), ch17:249 ("비싼 장비는 많이 써야 본전").

**[C7] Sentence and paragraph length, and a change of voice.** Rough counts on prose only: sentence-final "다." per 1,000 characters is 9–12 in ch01–05, 4–7 in ch08–16, 9.0 in ch17; " — " per 1,000 characters is 1.4–2.3 in ch01–05 and 3.3–4.9 in ch07–16. Units over 25 words between full stops: 2–11 % in ch01–05, 30–43 % in ch09–16. The "numpy로 확인한 것: (1)…(6)" paragraph of section 1 in ch08–17 is one paragraph of 6–11 sentences (rule: 3). Part 1 and ch17 read as prose; ch07–16 read as dense notes.

**[C8] Process narration inside the book text.** "작업계획서 §7 메모 참조" in section 10 of ch01–11; "다음 세션에서 실행해 … 이 자리에 넣는다" (ch01:86); "이 세션" ×5 and "이전 판" ×4 in ch14 (l.33 ×2, 73, 104); "처음 시도에서 약 69 GB 할당 실패 후 재작성" (ch11:257); "이 부호를 처음에 틀려 …" (ch12:86). A reader has no "이전 판". Move these to the plan's §8 or to one "집필 노트" box per chapter.

**[C9] Header drift.** ch17:5 "3부 알고리즘을 찾는 방법론" vs ch15/16 "3부 알고리즘 발견 방법론" (plan wording). Goal sentences of ch01–06 copy plan §4 verbatim; ch07–17 restate and extend them (ch12 adds energy QPE and amplitude estimation, ch13 adds diagnostics) without the plan table being updated.

**[C10] Template divergence not recorded.** Plan §5 has a section "하드웨어 실행(2·4·5부)"; ch08–13 have none (replaced by ch14 and a section "10. 출처와 환경"). Plan item 5 expects "StatevectorSampler/Estimator 출력"; every section 5 is numpy output. Needs one §8 line.

**[C11] Small cross-reference slips.** ch14:33 "위 기여 표의 유휴 몫" — the table is below, in §5. ch12:227 "큐비트화, Low–Chuang 2019 … 게이트 수가 O(αt + …)" vs ch16 "호출 수 … (Low–Chuang 2017)": it is a query count, and the two chapters cite different papers for it. ch08:236 "7장의 T2가 이 책에서 가장 중요한 하드웨어 숫자" vs ch14:108 "게이트 오류와 판독 오류가 지배한다". ch13:9 "10–12장의 알고리즘은 … 오류 정정 전에는 돌릴 수 없다" vs ch14 predicting a 3-bit QPE at 0.88. ch17 §3: "구조 없음" is both a drop signal (l.32) and a mapped row with a part (l.50), while the equally quadratic "확률·기댓값" row proceeds to step 5.

**[C12] ch12 and ch16 summaries on amplitude estimation.** ch12:11 says the error "1/M(Heisenberg)으로 줄고 고전은 1/√M" by comparing the BHMT bound (holds with probability ≥ 8/π²) with a classical standard error; ch16 §5/§7 shows the single-run RMSE does not scale as 1/N. Add a pointer in ch12:
> (이 비교는 확률 0.81의 상한과 표준오차를 견준 것이다. 한 번 실행의 RMSE가 다르게 움직이는 것은 16장 5절.)

---

## Part 2 — Coherence read-through

**Consistent across chapters (checked, no finding).** TFIM convention H = −JΣZZ − hΣX, open chain, J = h = 1 in ch12, 13, 14, 16, 17; E₀(n = 4) = −4.7588 in ch12, 13, 16; α = Σ|cᵢ| = 7 and ‖H‖ = 4.76 in ch16; HVA p = 1 energy −4.646 in ch14 = E₀ + 0.113 from ch13's table; n = 14 ground-state entropy 0.593 bits in ch15 and ch17. Heron r2 figures (CZ 68 ns, 1q 32 ns, readout 1.0–1.6 %, T1/T2) in ch05, ch07, ch14. Grover n = 3: 0.945 at k = 2 and 0.330 at k = 3 in ch09 and ch14; CZ 32 for the 8-CX Toffoli in ch09 and ch14. RSA-2048: 20 M qubits / 8 h (2019) and < 1 M / < 1 week, 1,400 logical qubits, 6.5×10⁹ Toffoli (2025) in ch11 and ch17. Error-budget limits "CZ ≈ 100 → 0.3, ≈ 230 → 0.5" in ch14 and ch16. Trotter-2 O(t^{1.5}/√ε) and qubitization O(αt + log(1/ε)/log log(1/ε)) in ch12 and ch16; both chapters and ch17 agree that Trotter is what runs today and that small problems favour Trotter. 8/π² vs 4/π² in ch10, 12, 16, 17. State-vector memory (17 GB at 30 qubits) in ch03, ch09, ch11. ch17's resource table and exercises 2 and 3 reproduce from its own formulas (recomputed by hand: 603 Toffoli/W, 4.9×10⁶, b = 13/15/16, N* = 10²², the 12 break-even entries).

**Contradictions and double definitions.** A1, B1, B2, B4, B10, B11, C11, C12 above.

**Concepts introduced twice.** Amplitude estimation is fully introduced in ch12 §3 (BHMT bound, experiment, exercise) and again in ch16 §3 with different symbols (t → m, "Q 호출" → "A 호출 N"); ch16:86 acknowledges ch12. Barren plateaus are introduced in ch06 §5/§7 and again in ch13 §3/§5, compatibly. Quantum walk appears in ch15 (Hadamard walk) and ch16 (Szegedy), linked. None of these conflict; ch12's AE material is simply in the wrong chapter for a reader looking for "해밀토니안 시뮬레이션".

**Is the thread 15 → 16 → 17 explicit?** Yes, and it is the best-joined part of the book. ch15:9 ("3부의 출발점 … 17장 체크리스트"), ch16:9 ("15장이 '가속은 어디서 오는가'였다면, 16장은 …"), ch17:9 ("15장은 …, 16장은 … 이 장은 둘을 묶어 체크리스트로") each state their place; ch17:41 rebuilds ch15's table as structure → part; ch16's cost units (calls, α) are the inputs of ch17's calculator. Three breaks in the thread: the wrong step number at ch15:59 (B9), the missing heuristic branch (B3), and the fact that parts 1–2 promised things to "3부" that none of the three chapters picked up (B9).

**What ch17 leaves for ch18 ("탐색 도구와 반례"), concretely.**
Must cover, with the promise it settles:
1. *Dequantization check* — promised at ch15:229, ch16:315, ch17:13 and by the "dequantization 여부" cell of ch17's record form (l.102). A runnable test: given a matrix or data set, measure its rank profile and run a sample-and-query (ℓ²-sampling) estimate against the block-encoding route of ch16 on the same low-rank input.
2. *Automatic classical-boundary triage* — promised at ch15:110 for ch17 and not done there. One function that takes a circuit and reports Clifford-only / T-count / MPS bond growth, plus the free-fermion (matchgate) test that the book's own TFIM example fails (A1, B1).
3. *A heuristic track for the checklist* — ch17:227 hands variational algorithms to ch18; ch06:211 promised a barren-plateau check as a checklist item; ch13 supplies the diagnostics (DLA/orbit dimension, gradient-variance slope, shot budget) and the counter-argument (Cerezo 2023). ch18 has to turn these into steps 3′–5′ so that ch19 can judge plan candidates 2, 3, 5, 6, 9; without it the ch17 procedure cannot evaluate three of the five chosen cases (B3).
4. *Kernel pre-screen* — ch06:45 and :213 defer "is the feature map hard and useful" to ch22; ch18 should give the test ch22 will apply (kernel concentration, comparison with an RBF baseline) and say where kernels sit in ch15's table.
5. *Automated circuit synthesis* — the plan's first tool. Open items it can close: ch04's unexecuted Solovay–Kitaev numbers, ch11's "Beauregard 회로의 게이트 수" (listed as unfinished in plan §7), ch08:226's `_decompose_` promise, and the missing T-counts of B6. A numpy gate counter that emits "2q / depth / T (or rotations)" for the part-2 circuits would do all four.
6. *Variational search as a discovery tool* (plan) — searching for a circuit, as opposed to ch13's optimising a fixed ansatz; with the plan's tools PennyLane and CUDA-Q.
7. *One "known algorithm revisited" counter-example* — ch11:227 names this as a part-3 theme (Regev 2023); ch17 §7 mentions DQI in one sentence.
8. *The sign problem / QMC boundary in one worked check* — ch17's table drops 2D TFIM as "stoquastic → QMC" on the strength of one phrase (ch15:57, ch17:223); if ch19 is to reuse that verdict, ch18 should show how to test it.

Should not repeat: the Tang and random-circuit-sampling stories (ch15 §7); the Babbush break-even argument, already told in ch09 §7, ch16 §7 and ch17 §3/§5 (fix B4 instead); barren-plateau theory and the n = 2–8 variance experiment (ch06 §5, ch13 §5); the stabilizer-tableau, T-branch and MPS-truncation experiments (ch15 §5); the qubitization resource calculator and the QSVT-vs-Trotter table (ch17 §5, ch16 §5); the Kim et al. 2023 episode (ch12 §7).

---

## Five structural weaknesses, each with a one-session fix

**1. The book's main SDK has never run.** 51 SDK snippets (Qiskit family 17, Cirq 11, PennyLane 11, CUDA-Q 5, QuTiP 5, Qualtran 2, pyqsp 1), 0 executed; every result in every section 5 is numpy. Principle 3 ("주 SDK는 Qiskit"), template item 5 and checklist item 1 are unmet in all 17 chapters, and the readable code of ch06–17 is pseudo-code (B7). The chapters are honest about it, but the book currently teaches numpy and annotates SDKs.
*Fix (offline, one session):* write `tools/extract_snippets.py` that pulls every ```python block into `snippets/chNN_k_<sdk>.py`, checks syntax, and turns each trailing `# → value` comment into an assertion. It makes "PyPI opens" a single command, and it would have caught the 7 non-parsing cells today. In the same session cut section 6 to one SDK per chapter (B8), which removes about 17 snippets from the backlog.

**2. Part 2's hardware promise is carried by one chapter that could not keep it, and other chapters have not been told.** ch14 is honest; ch04, 07, 08, 09, 10, 13 still say it "올린다"; ch12 and QAOA have neither noise simulation nor a place in ch14; the plan was not amended.
*Fix:* apply the B5 table, add the §8 line, and add two rows (Trotter-2 one step from ch12, QAOA p = 1 on the 6-ring from ch13) to `ch14_hardware.py` and the §5 table.

**3. The method of part 3 cannot evaluate most of part 4.** Variational and kernel approaches are absent from ch15's taxonomy, ch16's dictionary and ch17's mapping, while ch21, 22, 23 are built on them and the plan's candidate table calls them 가속 원천 (B3).
*Fix:* the two table rows of B3 now; make the heuristic track (item 3 of the ch18 list) the spine of ch18; correct the "가속 원천" column of the plan's candidate table to "근거" with the values 증명된 가속 / 휴리스틱.

**4. The running example is classically trivial and the book does not say so.** 1D TFIM appears in six chapters as the test bed, is called classically hard once (A1) and is the sole evidence for "dynamics is hard" (B1). More broadly, section 7's classical-simulability check is missing or nominal in ch04, 05, 07, 09, 11.
*Fix:* the A1 and B1 texts; one exercise in ch15 (below); re-run ch17 (A)-dynamics with a longitudinal field.
> **문제 (자작).** 12장의 TFIM 사슬에서 Majorana 연산자 γ₂ⱼ = (∏ᵢ<ⱼ Xᵢ)Zⱼ, γ₂ⱼ₊₁ = (∏ᵢ<ⱼ Xᵢ)Yⱼ를 정의하면 H가 γ의 2차식이 됨을 보이고, |0…0⟩에서 ⟨Z₀⟩(t)를 2n×2n 행렬 지수 하나로 계산해 12장 5절의 0.178(n = 4, t = 6)을 재현하라. 이 모형이 15장의 어느 고전 경계에 걸리는가?

**5. Shared conventions were set early and then drifted.** Cost reporting (B6), notation (B10), one-SDK rule (B8), source grades (C3), runnable cells (B7), prose rhythm (C7) and the session narration (C8) all hold in ch01–06 and loosen from ch07; ch17 partly returns to the early style. The result reads as three hands: part 1, ch07–16, ch17.
*Fix:* a front-matter page "표기와 비용 규약" (counting bits t, call count M_Q, number of solutions M, overlap γ = |⟨ψ₀|φ⟩|², normalisation α, precision ε, bond dimension χ, the three cost numbers and the rotation convention of B6), then a mechanical pass over ch10, 12, 16, 17 for symbols and over ch07–16 for the "numpy로 확인한 것" paragraph, broken into a numbered list.

---

## Summary

The 17 chapters are structurally uniform (10 sections, 교과서 비교, three answered exercises, 학습 노트 everywhere) and numerically consistent with one another on the quantities I traced. One statement is factually wrong and propagates: 1D TFIM is called classically hard (ch12) and is then used as the evidence that dynamics is the hard side (ch17), although it is free-fermion solvable — verified numerically. The plan's rules hold in ch01–06 and erode afterwards: one other SDK (fails ch07–17), runnable code cells (7 do not parse), the three cost numbers (1 circuit complies), source grades (stop at ch11). ch14 is honest about having no QPU, but ten sentences in six other chapters still describe it as a hardware run and the plan has no decision recording the change. Part 3 is well joined internally, but it has no slot for variational and kernel methods, which three of the five part-4 cases use; ch18 has to supply that track.

Counts: A 1 · B 12 · C 12.

**Checked and found correct:** section inventory and order in all chapters; 교과서 비교 boxes; Quirk links in section 2 (ch11 the only chapter without one, reason stated); three exercises with in-place answers; 학습 노트; ch14's statements about the missing QPU run (consistent in l.9, 110, 118, 184, 206, 218); the cross-chapter quantities listed under "Consistent across chapters"; ch17 resource table, break-even table and exercises 2–3 (recomputed); ch12's AE bound, crossover ≈ 50 and exercise 3 numbers (recomputed); ch09's Babbush Table I numbers (against the paper); ch16's internal Trotter-vs-QSVT ratios (§5 vs §7); CUDA-Q 0.16 exists.

**Could not check:** any SDK snippet by execution (no PyPI); the hand-counted BV depth in the B6 replacement text; ch14's empirical device statements (B12); 2026 items in ch11 §7 and ch13 §3 (Iceberg Quantum, Google ECC-256 estimate, Marwaha 2026); the html builds and the Docs copies; sentence-length statistics are heuristic (regex on "다."), good for the trend only. I read ch02–11 selectively, so a contradiction confined to their sections 3–5 could have been missed.
