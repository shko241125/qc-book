# R1 findings — exercises and worked answers, ch01–09 (plus §7 and §9 claims)

Reviewer id: R1. Scope: all 27 exercises in "## 8. 연습문제와 해답" of ch01–ch09, solved independently
(numpy/scipy scripts or algebra) before comparing with the book's answer; plus the factual claims in §7 and §9 of
each chapter; plus textbook attributions on exercises.
Scratch code: `/tmp/claude-0/-home-claude/fd63146d-511d-5ae9-8dab-3100fb85f0db/scratchpad/mid/R1/ex01_04.py`, `ex05_06.py`, `ex07_09.py` (all run; outputs quoted below).
Nothing under `/home/claude/qc-book/` was modified.

Headline: **all 27 exercise answers are numerically correct** (every number, sign, phase, gate count I recomputed
agrees). The findings below are one wrong figure in ch09 §7, three logic/example errors inside otherwise-correct
answers, three wrong statements about sources, and a set of minor items.

---

## Severity A

### A1. ch09 §7 — `ch09/ch09_grover.md:232`
Current: "Campbell–Khurana–Montanaro 2019도 SAT에 Grover를 적용할 때 물리 큐비트 수백만 개가 필요하다고 추산했다."

Wrong by about six orders of magnitude. The paper's tables give ~10¹²–10¹⁴ physical (magic-state-factory) qubits:
14-SAT with Grover 1.38×10¹² (optimistic) to 3.14×10¹³ (realistic); 12-SAT backtracking 2.17×10¹³ to 3.85×10¹⁴; graph
colouring 3.59×10¹² to 6.29×10¹³. The text says "very large (e.g. over 10¹²)", and the abstract says "the number of
physical qubits used is extremely large" and that the advantage disappears once classical decoding cost is included.
Evidence: https://ar5iv.labs.arxiv.org/html/1810.05582 , https://quantum-journal.org/papers/q-2019-07-18-167/

Proposed replacement:
> Campbell–Khurana–Montanaro 2019도 SAT·그래프 색칠에 Grover와 양자 백트래킹을 적용하면 물리 큐비트가 10¹²개 이상(마법 상태 공장 기준 10¹²~10¹⁴개) 필요하고, 표면 부호 복호에 드는 고전 연산 비용까지 넣으면 양자 이득이 사라진다고 추산했다.

(The Babbush et al. numbers in the same paragraph — 10³ cores: 100일~880년, 10⁶ cores: 280년~88만 년 — were checked
against Table I of arXiv:2011.04149 and are correct.)

---

## Severity B

### B1. ch08 §8 문제 1 해답 — `ch08/ch08_oracles.md:246`
Current: "일반 균형 f의 위상 오라클은 얽힘 게이트이며(예: f = x₁x₂ → CZ)"

f = x₁x₂ is not a balanced function (2 bits: one 1 out of 4; 3 bits: two 1s out of 8 — computed in `ex07_09.py`), so
it is not a valid DJ input and cannot illustrate "balanced f". A balanced example with an entangling phase oracle is
f = x₁x₂ ⊕ x₃ (four 1s out of 8).

Proposed replacement:
> 일반 균형 f의 위상 오라클은 얽힘 게이트이며(예: 3비트 균형 함수 f = x₁x₂ ⊕ x₃ → CZ₁₂·Z₃)

### B2. ch09 §8 문제 2 해답 and §5 — `ch09/ch09_grover.md:250` and `:152`
Current (250): "M > N/2: θ > π/4, k = 0이 최적 — 그냥 재라."  Current (152): "M이 N/2를 넘으면 θ > π/4라 k = 0이 최적"

As a statement about the success probability this is false. P_k = sin²((2k+1)θ): N = 8, M = 5 gives
P₀ = 0.625, P₁ = 0.156, P₂ = 0.977; N = 16, M = 9 gives P₀ = 0.562, P₄ = 0.952, P₆ = 0.999 (`ex07_09.py`).
What is true: for θ ∈ (π/4, π/2), sin²3θ − sin²θ = sin4θ·sin2θ < 0, so one iteration always lowers the probability;
and plain sampling already succeeds in fewer than 2 expected tries, so amplification is pointless. The conclusion
("그냥 재라") stands, the reason given does not.

Proposed replacement (250):
> M > N/2: θ > π/4이면 sin²3θ − sin²θ = sin4θ·sin2θ < 0이라 반복 한 번은 확률을 오히려 낮춘다. k ≥ 2에서 다시 높아질 수는 있지만(N = 8, M = 5: 0.625 → 0.156 → 0.977) 그냥 재고 검증하면 기대 2회 미만에 끝나므로 증폭할 이유가 없다 — 그냥 재라.

Proposed replacement (152):
> M이 N/2를 넘으면 θ > π/4라 ⌊π/(4θ)⌋ = 0이고 반복 한 번은 확률을 오히려 낮춘다 — 그냥 재는 게 낫다(기대 2회 미만).

### B3. ch09 §8 문제 3 해답 — `ch09/ch09_grover.md:254`
Current: "닿은 뒤는 매 라운드 확률 ≥ 1/4로 성공하므로 기대 라운드 ≤ 4, 각 라운드 ≤ √N 호출. 합쳐 O(√(N/M))"

Non sequitur: "≤ 4 rounds × ≤ √N calls" gives O(√N), not O(√(N/M)), and m keeps growing after the critical stage, so
the round cost is not bounded by m₀ either. The BBHT argument is a geometric series in which the growth λ per round is
beaten by the failure probability 3/4 (this is why λ < 4/3 is required; the paper says any 1 < λ < 4/3 works and gives
3m₀ for the pre-critical stage — https://ar5iv.labs.arxiv.org/html/quant-ph/9605034). The identity
(1/m)Σ sin²((2k+1)θ) = 1/2 − sin(4mθ)/(4m sin2θ) and the m ≥ 1/sin2θ ⇒ P ≥ 1/4 step in the answer are correct
(checked numerically for m = 1, 3, 7, 20). My 200 000-trial rerun of the book's protocol gives 26.1 / 12.7 / 5.76 calls
for M = 1 / 4 / 16; the book's 25.6 / 12.7 / 5.7 (200 trials, s.e. ≈ 1.0 / 0.5 / 0.25) is consistent.

Proposed replacement for the two sentences from "전략:" to "합쳐 O(√(N/M))":
> 전략: m = 1에서 시작해 실패할 때마다 m을 λ = 6/5배(√N 상한)로 키운다. 한 라운드의 기대 호출 수는 약 m/2다. m이 임계값 m₀ = 1/sin 2θ ≈ √(N/M)/2에 닿기 전의 비용은 등비급수 합으로 (1/2)·λ/(λ−1)·m₀ = 3m₀보다 작다. 닿은 뒤 u번째 라운드는 비용이 λ^{u+1}·m₀/2로 커지지만 거기까지 실패가 이어질 확률이 (3/4)^u 이하이므로, 기대 비용은 (λm₀/2)·Σ(3λ/4)^u = 6m₀로 수렴한다(3λ/4 = 0.9 < 1 — λ < 4/3이어야 하는 이유). 합쳐 9m₀ ≈ (9/2)√(N/M), 즉 O(√(N/M)).

(The "(9/2)√(N/M)" constant itself is from my memory of BBHT plus the 3m₀ + 6m₀ arithmetic; the fetched page was
truncated before the final total — medium-high confidence.)

### B4. ch04 §3 and §9 — `ch04/ch04_universality.md:39` and `:226`
Current (39): "교환자가 오차를 제곱으로 줄이는 성질이 핵심이다."  Current (226): "군 교환자 VWV⁻¹W⁻¹가 오차를 제곱으로 줄이는 성질이 재귀의 엔진. 비유: Newton 법이 오차를 제곱으로 줄이는 것과 같은 구조."

The Solovay–Kitaev recursion is ε_n = c·ε_{n−1}^{3/2}, not quadratic. The book's own exponent proves it: sequence
length ×5 per level and c ≈ 3.97 = ln 5 / ln(3/2); quadratic convergence would give ln 5 / ln 2 = 2.32. (What is
quadratic is that the commutator of two elements δ-close to I is δ²-close to I; the net error gain per level is the
3/2 power.) Source: Dawson & Nielsen, arXiv:quant-ph/0505030 (from memory; the arXiv mirror was not fetchable here,
the ln 5/ln 1.5 arithmetic is the verifiable evidence).

Proposed replacement (39):
> V, W의 근사 오차 ε이 교환자에서는 O(ε^{3/2})로 줄어드는 성질이 핵심이다(한 단계마다 길이 5배, 오차 지수 3/2배 → c = ln 5/ln(3/2) ≈ 3.97).

Proposed replacement (226):
> 군 교환자 VWV⁻¹W⁻¹가 오차를 ε → O(ε^{3/2})로 줄이는 성질이 재귀의 엔진. 비유: Newton 법처럼 단계마다 오차의 지수가 곱으로 커지는 초선형 수렴(다만 차수는 2가 아니라 3/2).

### B5. ch08 §9 — `ch08/ch08_oracles.md:262`
Current: "Simon·Shor 후처리와 안정자 표현(Gottesman–Knill, 3장)에 같은 코드가 쓰인다."

Shor's post-processing is continued fractions, not GF(2) elimination; ch11 contains no GF(2) code (grep for
`gf2|GF(2)|가우스` in ch11 md/py: no hits), and ch08 §7 itself says "Shor의 연분수 후처리(11장)".

Proposed replacement:
> Simon 후처리와 안정자 표현(Gottesman–Knill, 3장)에 같은 코드가 쓰인다(Shor의 후처리는 GF(2) 소거가 아니라 연분수다 — 11장).

### B6. ch04 §3, §5 (attribution) — `ch04/ch04_universality.md:41`, `:149`; script `ch04_universality.py:68,85-86`
Current (41): "CNOT 6개, T·T† 7개, H 2개로 정확히 분해된다(N&C 그림 4.9)"  Current (149): "N&C 그림 4.9의 15개 게이트 순서를"

N&C Figure 4.9 is a 16-gate circuit that also contains one S gate (its caption reads "…using Hadamard, phase,
controlled-NOT and π/8 gates"): 6 CNOT, 4 T†, 3 T, 2 H, 1 S. Confirmed from Quipper's `toffoli_NC_at`, documented as
"from Nielsen and Chuang (Figure 4.9)": https://mathstat.dal.ca/~selinger/quipper/doc/src/Quipper/Libraries/GateDecompositions.html
The 15-gate sequence in §4 (`…, T(1), T(2), H(2), CX(0,1), T(0), T†(1), CX(0,1)`) is the S-free variant, identical to
Qiskit's `ccx` definition. The matrix identity is correct (re-verified); only the attribution and count are off.

Proposed replacement (41):
> CNOT 6개, T·T† 7개, H 2개로 정확히 분해된다(N&C 그림 4.9와 동치인 판. N&C 그림은 제어 쪽에 S 하나를 더 쓴 16게이트이고, 여기서는 그 S를 T 배치로 흡수한 15게이트 판 — Qiskit `ccx`의 정의 — 을 쓴다).

Proposed replacement (149):
> N&C 그림 4.9와 동치인 15개 게이트 순서(S 없는 판, Qiskit `ccx` 정의)를 행렬로 곱한 결과가 CCX와 `allclose`로 일치했다.

Also adjust the two print labels in the script ("N&C Fig.4.9", "N&C 4.3절의 표준 분해").

### B7. ch06 §3 (attribution) — `ch06/ch06_circuit_library.md:49`
Current: "산술 회로는 N&C에 없다(연습 5.6은 근사 QFT; 본문은 고전 가역 회로 중심) — Draper 2000 원논문이 1차 출처다."

N&C does contain the QFT adder idea: chapter-end **Problem 5.6 "Addition by Fourier transforms"** (|x⟩ → |x + y mod 2ⁿ⟩
by QFT, single-qubit phase shifts, inverse QFT). The parenthetical is right about *Exercise* 5.6 (approximate QFT) but
overlooks *Problem* 5.6. Evidence: the qsimulatR vignette "Addition by Fourier transform" states "This corresponds to
problem 5.6 in Nielsen & Chuang" (https://archive.linux.duke.edu/cran/web/packages/qsimulatR/vignettes/addbyqft.html).

Proposed replacement:
> 산술 회로는 N&C 본문에 없다(본문은 고전 가역 회로 중심). 다만 5장 장말 문제 5.6 "Addition by Fourier transforms"가 상수를 QFT로 더하는 같은 아이디어를 묻는다 — 연습 5.6(근사 QFT)과는 다른 항목이다. 두 레지스터 덧셈기의 1차 출처는 Draper 2000 원논문이다.

---

## Severity C

### C1. ch09 §9 — `ch09/ch09_grover.md:260`
Current: "더 깊이: Brassard–Høyer–Mosca–Tapp 2002, KLM 8.1." Amplitude amplification is KLM 8.2 by the chapter's own §3
(line 58) and §10 (line 271); 8.1 is plain Grover.
Proposed: `더 깊이: Brassard–Høyer–Mosca–Tapp 2002, KLM 8.2.`

### C2. ch06 §8 문제 2 (attribution) — `ch06/ch06_circuit_library.md:223` (and `:245`)
Current: "문제 2 (N&C 5.2절 준비)". N&C 5.2 is phase estimation; a period comb → peaks at multiples of N/r is the
content of 5.4.1 (period-finding) / 5.3.1 (order-finding). Confidence: high on section titles (from memory).
Proposed: `**문제 2 (N&C 5.4.1절 주기 찾기 준비).**`; in §10: `5.2절(위상 추정), 5.4.1절(주기 찾기 — 연습문제 2)`.

### C3. ch05 §8 문제 1 해답 — `ch05/ch05_measurement.md:186`
Current: "|Ψ⁻⟩ = (X Z ⊗ I)|Φ⁺⟩(전역 위상 제외)이므로 모든 경우에 XZ가 추가로 붙는다."
The table (0,0)→XZ, (0,1)→Z, (1,0)→X, (1,1)→I is correct (verified by brute force, with the extra Pauli on either
half). The reasoning skips a step: XZ⊗I acts on Alice's half (q₁), while the correction is on Bob's qubit; it carries
over only via (A⊗I)|Φ⁺⟩ = (I⊗Aᵀ)|Φ⁺⟩, which the chapter never states.
Proposed:
> |Ψ⁻⟩ = (I ⊗ XZ)|Φ⁺⟩(전역 위상 제외) — Bell 쌍에서는 (A⊗I)|Φ⁺⟩ = (I⊗Aᵀ)|Φ⁺⟩이고 (XZ)ᵀ = −XZ이므로 Alice 쪽 q₁에 걸든 Bob 쪽 q₂에 걸든 같은 상태다. 따라서 Bob의 큐비트에는 모든 경우에 XZ가 추가로 붙는다.

### C4. ch01 §8 문제 2 해답 — `ch01/ch01_qubit_statevector.md:185`
(a) "5절의 2000샷에서 ±0.027이 나온 이유가 여기 있다": the same formula gives a 95 % half-width of
1.96·√(0.25/2000) = 0.0219 at 2000 shots; the observed +0.027 is outside it (2.4σ, as §5 says).
(b) The problem says "`sample()`로 확인하라" but the answer shows no check. Mine: 20 000 repetitions at N = 38 416 →
95.3 % within ±0.005; at N = 2000 → 32.8 %.
Proposed (replace that sentence):
> 같은 식으로 2000샷의 95 % 폭은 ±0.022이므로, 5절의 +0.027은 그 바깥에 놓인 2.4σ짜리 드문 요동이다. 확인: N = 38 416으로 `sample()`을 여러 번 돌리면 약 95 %가 0.5 ± 0.005 안에 든다(2000샷에서는 약 33 %).

### C5. ch07 §8 문제 3 해답 — `ch07/ch07_physical_layer.md:247`
Current: "순도와 Bloch 길이 r의 관계 Tr ρ² = (1 + r²)/2가 맞는다(3장)". The numbers (0.75, 0.707; then 0.625, 0.5)
are correct, but ch03 never states this relation (grep `r²` in ch01–07: only this line).
Proposed: `순도와 Bloch 길이 r의 관계 Tr ρ² = (1 + r²)/2가 맞는다(ρ = (I + r·σ)/2를 제곱해 대각합을 취하면 나온다).`

### C6. ch03 §7 — `ch03/ch03_multi_qubit.md:208`
Current: "이 장의 게이트(H, S, CNOT, CZ, SWAP)는 전부 Clifford 게이트이고". The chapter's own gate table (line 57)
and subtitle include Toffoli, which the same paragraph then lists as non-Clifford; S is not used in the chapter.
Proposed: `이 장에서 Bell·GHZ를 만드는 데 쓴 게이트(H, X, Z, CNOT, CZ, SWAP)는 전부 Clifford 게이트이고(표의 Toffoli만 예외다)`

### C7. ch07 §7 vs §9 — `ch07/ch07_physical_layer.md:231` and `:255`
§7: coherent errors grow "n배(확률은 n²배)", stochastic noise "√n으로 쌓이는"; §9: "보정 오차는 n번에 n²배, 채널 노이즈는
n배". Both are right, but §7 quotes the stochastic case as an amplitude and §9 as a probability, which reads as a
contradiction.
Proposed (231): `독립적인 채널 노이즈가 진폭으로는 √n배(확률로는 n배)로 쌓이는 것과 다르다.`

### C8. ch09 §7 and §9 — `ch09/ch09_grover.md:232`, `:262`
Current: "오류 정정은 게이트당 시간을 고전 연산보다 수백만 배 느리게 만든다" / "게이트당 시간 수백만 배".
Babbush et al. use a Toffoli time of ≈170 μs against a 330 ps clock cycle: ≈5×10⁵, i.e. "수십만 배"
(https://ar5iv.labs.arxiv.org/html/2011.04149). Order-of-magnitude wording only.
Proposed: `수십만 배(Babbush et al.: 논리 Toffoli ≈ 170 μs 대 CPU 클럭 ≈ 0.33 ns)`

### C9. ch06 §8 문제 1 해답 and §9 — `ch06/ch06_circuit_library.md:221`, `:234`
Current: "11장에서 모듈러 지수의 보조 레지스터를 되돌릴 때(uncompute) 이 성질을 쓴다." ch11 implements U_a as a permutation
matrix and mentions the Draper adder only in its cost note (`ch11_shor.md:89`); no uncompute is carried out there.
Proposed: `11장이 비용만 인용하는 Beauregard 회로가 모듈러 덧셈의 비교·보정 단계에서 보조 큐비트를 되돌릴 때(uncompute) 이 성질을 쓴다.`

### C10. ch04 §8 문제 3 해답 — `ch04/ch04_universality.md:221`
Current: "이 패턴이 12장 Trotter 회로의 exp(−iθ ZZ) 게이트이고". Only the CNOT·Rz·CNOT core is exp(−iθZZ); the full
CRz pattern has an extra outer Rz. Verified: CNOT·(I⊗Rz(2θ))·CNOT = exp(−iθ Z⊗Z); the CP = CRz·(P(λ/2)⊗I) remark is
correct.
Proposed: `이 패턴의 가운데 토막 CNOT·Rz·CNOT이 12장 Trotter 회로의 exp(−iθ ZZ) 게이트이고(CNOT·(I⊗Rz(2θ))·CNOT = exp(−iθ Z⊗Z))`

### C11. ch05 §8 문제 3 — `ch05/ch05_measurement.md:192-194`
The problem asks for a channel on a density matrix and the answer uses "Kraus 연산자" and "추적 보존", neither of which
the chapter body introduces (they appear only in this answer and in §9). The answer is correct and self-contained;
a one-line definition of ρ → Σ KₖρKₖ†, Σ Kₖ†Kₖ = I in the problem statement would make it solvable as posed.

### C12. Exercise checks quoted as output but not produced by the chapter scripts (ch02–ch06)
All confirmed by my recomputation, so nothing is wrong, but these have no script behind them (grep of scripts and
`chNN.out`; only ch07–09 scripts carry exercise checks):
- ch02 문제 1 `(True, 0.3927)`; 문제 2 `(True, −0.7854)` — and `SX` is not defined in the ch02 script or §4 snippet; 문제 3 `True`.
- ch03 문제 1 orthonormality; 문제 2 "≈ 0.601 bit" (I get 0.6009); 문제 3 `kron(minus, minus)`.
- ch04 문제 1 BFS count 24 (I get 24; {H,T} exceeds a 5000 cap); 문제 2 and 문제 3 `allclose` checks.
- ch05 문제 1 correction table; 문제 3 Kraus check.
- ch06 문제 1 (b = 011, b = 101); 문제 2; 문제 3 and §7 line 213: K(x,y) = 0.955/0.962/0.981, K(x,z) = 0.071/0.725/0.804
  — the script prints reps = 2 only (I get 0.9552/0.9620/0.9814 and 0.0712/0.7247/0.8040).
Suggest adding a short "연습문제 검산" block to each of these scripts, as ch07–09 already have.

---

## Summary

- 27 of 27 exercise answers reproduce numerically; no wrong number, sign, phase, endianness or gate count in any answer.
- One A: the Campbell–Khurana–Montanaro qubit count in ch09 §7 ("수백만" vs the paper's 10¹²–10¹⁴).
- Seven B: a non-balanced function offered as a balanced example (ch08 문제 1); "k = 0이 최적" for M > N/2 (ch09 문제 2, §5);
  a non sequitur in the BBHT cost argument (ch09 문제 3); Solovay–Kitaev described as quadratic (ch04 §3, §9);
  GF(2) elimination attributed to Shor post-processing (ch08 §9); N&C Fig. 4.9 gate count (ch04); "N&C has no
  arithmetic circuit" (ch06, Problem 5.6).
- Twelve C: section-number slips, a skipped step, cross-references to things other chapters do not contain, wording.
- B6, B7, B4 (§3 part) and C8 lie outside §7–§9 but came up while checking attributions.

## Checked and found correct (own computation or algebra agrees)

- ch01: 문제 1 (⟨+|−⟩ = 0, P(+) = 1); 문제 2 N = 38 416; 문제 3 (θ, φ) = (π/2, −π/2), Bloch (0, −1, 0). §7 H·H interference, global phase.
- ch02: 문제 1 HTH = e^{iπ/8}Rx(π/4); 문제 2 H = e^{iπ/4}Rz(π/2)√X Rz(π/2); 문제 3 SHS† = (Y+Z)/√2, |0⟩ → |+i⟩. §7 Bloch coordinates under Rz(0.3/1.0/2.5), P(1) = sin²(t/2). §9 Rz(2π) = −I.
- ch03: 문제 1 four Bell states, CNOT(H⊗I)|ab⟩ mapping, orthonormality; 문제 2 H₂(cos²θ), 0.601 bit, Ry(2θ)+CNOT; 문제 3 CNOT|+⟩|−⟩ = |−⟩|−⟩. §7 memory figure (2³⁰ × 16 B ≈ 17 GB), CHSH table values, §9 classical-mixture correlation cosθ₁cosθ₂.
- ch04: 문제 1 |⟨H,S⟩| = 24 modulo phase; 문제 2 CNOT = (I⊗H)CZ(I⊗H); 문제 3 CRz from 2 CNOT, CP = CRz·(P(λ/2)⊗I). HTHT angle 0.348886π. §7 Ω(2ⁿ log(1/ε)/log n).
- ch05: line-42 Bell-basis expansion; 문제 1 table; 문제 2 (1/√2 vs 1/2); 문제 3 K₀ = |0⟩⟨0|, K₁ = |0⟩⟨1|, trace-preserving, output |0⟩⟨0|. §5/§7/§9 F ≈ 1 − (4/3)p (exact 1 − (2/3)(2p − p²): 0.9867, 0.9606, 0.9350, 0.8733 — book's 300-sample values consistent); uncorrected mean fidelity 1/2.
- ch06: QFT = DFT for n = 2, 3, 4; Draper listing = adder on all 64 inputs, 18 gates; 문제 1 subtractor (b = 3, b = 5); 문제 2 peaks (0.25 at k = 0, 2, 4, 6; 0.5 at k = 0, 4); 문제 3 six kernel values; kernel is identical under either sign convention of the feature-map exponent. §9 parameter-shift formula.
- ch07: 문제 1 Rx(2π) = −I and controlled version sends |+⟩ → |−⟩; 문제 2 derivation, T2 = 85.714, |ρ₀₁|(70) = 0.22095 = 0.5e^{−70/T2}, Kraus ≡ Lindblad to 1e−15 with γ = 0.5034; 문제 3 (0.3536, 0.75, 0.7071) then (0.25, 0.625, 0.5). §5 1−F table reproduced (0.0055, 0.0270, 0.0526, 0.0996, 0.1794); 2.7×10⁻⁴ estimate.
- ch08: 문제 1 U_φ = ⊗Z^{sᵢ}, output |s⟩; 문제 2 (1 − k/2^{n−1})² for k = 0…8, 1/16 on the other seven outcomes; 문제 3 (a) 0.3076, limit 0.2888, (b) 5.543, constant 1.6067, (c) exact 4.58×10⁻⁴ < 2⁻¹¹ (Monte Carlo 200 000 runs: 5.542, 0.306, 3.8×10⁻⁴). BV noise formula values. §7 2^{−19}, QSD count.
- ch09: 문제 1 P = 1 for all four marked items with the stated circuit, CZ 2, single-qubit gates 10–14, classical average 2.25; 문제 2 M = N/4 → 1, M = N/2 → 1/2; 문제 3 averaging identity and simulation. §3 Toffoli/ancilla counts for Cⁿ⁻¹Z, §4 gate tallies (12/14, 24/28, 32, 16), §7 Babbush Table I numbers.
- Attributions confirmed (from memory, high confidence unless noted): N&C Ex. 4.13/4.14 (HTH = Rx(π/4)); N&C 1.3.6 (Bell states); N&C Ex. 4.20 (CNOT in ± basis); Corollary 4.2; N&C Ex. 1.2 and Box 12.1 (cloning); N&C 8.3.5–8.3.6 (amplitude/phase damping); N&C 6.1.4 (doubling the search space); N&C Ex. 5.6 = approximate QFT; KLM 6.5 (Simon); KLM 8.4 (BBHT). Journal references in §9/§10 of ch02–ch09 (volume/page/year) match what I know.

## Could not check

- "KLM 3장 연습 변형" (ch01 문제 3), "KLM 4장 연습 변형" (ch02 문제 3), "KLM 6장 변형" (ch03 문제 3): no exercise number is given and I have no access to KLM's exercise list.
- ch05 §9 "N&C 12.1.1절" for no-cloning: Box 12.1 sits in section 12.1 near the start of 12.1.1 (Holevo bound); unverified suspicion that "12.1절, Box 12.1" is the better pointer.
- ch01 §9 "KLM 3.4절 'The measurement postulate'": I recall the section title as "Measurement"; unverified.
- BBHT's final constant (9/2)√(N/t) and Dawson–Nielsen's ε^{3/2} wording: from memory plus arithmetic; primary PDFs were not fetchable (arXiv not on the allowlist; ar5iv page truncated).
- SDK snippets and Quirk links (not run; outside this angle).
