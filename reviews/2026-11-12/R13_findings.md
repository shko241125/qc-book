# R13 — independent technical review of chapter 21 (`ch21/ch21_hp_folding.*`)

Reviewer: R13 · 2026-11-12 (UTC) · read-only on the repository. All scratch code and outputs are in
`/tmp/claude-0/-home-claude/fd63146d-511d-5ae9-8dab-3100fb85f0db/scratchpad/r13/` (my own code; the chapter script was read but never imported or re-run).
Line numbers: `md:N` = `ch21_hp_folding.md`, `py:N` = `ch21_hp_folding.py`, `out:N` = `ch21_stdout.txt`.

**Counts: A 6 · B 13 · C 15.**

## Verdict

The classical half of the chapter is sound: cost tables, λ bounds, Walsh–Hadamard counts, random-sampling TTS, the simulated-annealing baseline, the small SDK case and the exercise arithmetic all reproduce with my own code. Every number in the md that the script prints matches `ch21_stdout.txt`.

The quantum half is not. **One sign in `ramp()` (py:250–253) sends every QAOA optimisation in the chapter to the wrong basin.** With the chapter's convention e^{−iγC/σ}·e^{−iβΣX} and a cost that is *minimised*, an annealing-like start needs γ and β of opposite sign. `ramp` returns both positive, which is the start for *maximising* C (it is the convention ch13 used for MaxCut maximisation). At that start ⟨C⟩ is above the uniform mean and p\* is 0.24–0.31× random. Nelder–Mead then stops in poor local minima, and those are what the chapter reports as "이상적인 경우". Changing `0.6 * (1 - f)` to `-0.6 * (1 - f)` in my replica of the chapter's own procedure gives ⟨C⟩ = 1.2965 / 0.4740 / 0.1031 at N = 8 (chapter: 2.958 / 0.651 / 0.327).

Consequences, all recomputed:

- "p = 1 lowers ⟨C⟩ but makes p\* worse than random" is an artefact. At the real p = 1 minimum p\* is 4.8× / 6.4× / 10.6× random (N = 8 / 10 / 12).
- "Even with zero training cost SA wins at N = 12 (24,724 vs 16,045)" is wrong. The ⟨C⟩-optimal p = 3 circuit has p\* = 78.8× random, TTS₉₉ = 15,321 — a statistical tie with SA (16,045; the script's own ±1 s.e. range is 15,351–16,795).
- The answer-aware value at N = 12 is not 81.2× but at least 262.8× (TTS₉₉ 4,592), shorter than the tuned SA (9,865).
- "Parameter transfer does not beat SA either" (md:343) is wrong. p = 3 angles fixed on the smaller instances and used unchanged at N = 12 (zero evaluations at N = 12) give 16,215–19,452 (⟨C⟩-tuned) or 5,963–9,011 (p\*-tuned on N = 8/10, where the answer is known by enumeration). Read literally, the pre-registered criterion (TTS including training shots) is **met** by the second protocol.
- The pre-registered verdict "불통과" still holds for the protocol the chapter registered as (다) — SPSA at N = 12 — also after the sign is fixed (my 4 seeds: 125k–275k). But the chapter's claim that the verdict does not depend on the training choices is false; it depends on exactly that.

What survives: the unit caveat (one shot = 3 layers × 10⁶–10⁷ two-qubit gates vs one table lookup) and the "방법 사례" role. The chapter needs a re-run and a rewritten §1, §5 (5)–(7), §7, §9, exercise 3 and figure; the author has to decide how the registered criterion treats classically pre-computed angles.

---

## A — wrong number/claim, or misleads about the verdict

### A1. `ramp()` has the wrong sign; all "(가) 이상적인 경우" values are poor local minima

- Location: py:250–253 (`ramp`), py:276–299 (N = 8 starts are `ramp(p)` and the stretched previous solution; N = 10, 12 start from the smaller size, so the bad basin propagates), py:327 (SPSA start), py:436 (case file start). md:14–15, md:104, md:186–198 (table and text), md:216–218, md:226–230, md:332 (exercise 3), md:342–343, figure panels (2)(4)(5), asserts py:511–513.
- What is wrong: under e^{−iγC}·e^{−iβΣX}, |+⟩ⁿ is the *top* state of +ΣX. With γ, β > 0 the layers follow the top of βΣX + γC, i.e. toward maximum cost. The chapter's own printed angles show the optimiser fighting this: N = 8, p = 3 is (+0.317 −0.546 −1.321 | +1.099 +0.189 +0.118); β₁ = 1.099 ≡ −0.472 (β → β − π/2 is the global bit flip, a symmetry of C — I checked numerically that it leaves ⟨C⟩, P(valid), p\* unchanged), so layer 1 has (γ, β) = (+, −) and layers 2–3 have (−, +): a mixed-orientation local minimum.
- Evidence (`b_qaoa8.py`, `c_multi.py`, `o_fix.py`, `g_big.py`, `p_more12.py`):
  - Start point `ramp(2)`: ⟨C⟩ = 5.742 / 8.182 / 9.744 at N = 8 / 10 / 12 against uniform means 3.562 / 5.252 / 6.348; p\*/random = 0.27 / 0.31 / 0.24. The script prints this (out:66, "출발점에서의 p\*: 3.910e-04; 4.395e-04", random is 1.465e-03 / 1.434e-03) but the md never mentions it.
  - My replica of the chapter's N = 8 procedure (same starts, Nelder–Mead, maxfev 400) with the chapter's sign reproduces the chapter exactly: 2.9582 / 0.6511 / 0.3272, ratios 0.55 / 7.52 / 10.36, angles `[0.317, -0.546, -1.321, 1.099, 0.189, 0.118]`. With the β sign flipped the same procedure gives 1.2965 / 0.4740 / 0.1031, ratios 4.78 / 9.71 / 14.48, angles `[0.449, 0.893, 1.615, -0.19, -0.153, -0.109]`.
  - p = 1 at N = 8 on a 241 × 121 grid over one full period (γ ∈ [−πσ, πσ], β ∈ [−π/4, π/4]): global minimum 1.2965 at (γ, β) = (0.935, −0.204). The chapter's 2.9582 is the third-best local minimum.
  - Multi-start BFGS at N = 8 (58–62 starts per p) finds nothing below 1.2965 / 0.4740 / 0.1031.

  | N | p | chapter ⟨C⟩ | R13 ⟨C⟩ | chapter P(valid) | R13 | chapter p\*/random | R13 | chapter TTS₉₉(훈련 0) | R13 |
  | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
  | 8 | 1 | 2.958 | 1.2965 | 0.190 | 0.478 | 0.55 | 4.78 | 5,678 | 655 |
  | 8 | 2 | 0.651 | 0.4740 | 0.641 | 0.687 | 7.52 | 9.71 | 416 | 322 |
  | 8 | 3 | 0.327 | 0.1031 | 0.730 | 0.802 | 10.36 | 14.48 | 301 | 215 |
  | 10 | 1 | 3.958 | 2.1218 | 0.121 | 0.315 | 1.07 | 6.36 | 2,999 | 502 |
  | 10 | 2 | 1.176 | 0.9369 | 0.475 | 0.532 | 10.75 | 13.56 | 296 | 235 |
  | 10 | 3 | 0.670 | 0.3578 | 0.586 | 0.663 | 15.17 | 20.5 | 209 | 154 |
  | 12 | 1 | 4.731 | 2.6604 | 0.066 | 0.183 | 0.62 | 10.6 | 1,931,997 | 114,382 |
  | 12 | 2 | 1.511 | 1.1746 | 0.305 | 0.360 | 25.56 | 35.5 | 47,235 | 33,973 |
  | 12 | 3 | 0.794 | 0.3863 | 0.432 | 0.498 | 48.82 | 78.8 | 24,724 | 15,321 |

  N = 12 values are Nelder–Mead from my N = 8 angles (50 / 110 / 420 evaluations). p = 3 is converged (a further 220 evaluations moved ⟨C⟩ by < 10⁻⁴); p = 2 may still improve slightly. All are local searches in the annealing-like basin; I do not claim global minima for N ≥ 10.
- Proposed fix:
  - py:250–253:
    ```python
    def ramp(p):
        """출발점: 선형 경사. e^{−iγC}·e^{−iβΣX} 에서 C 를 줄이려면 γ 와 β 의 부호가 반대여야 한다(|+⟩ 는 +ΣX 의 최고 상태)."""
        f = (np.arange(p) + 0.5) / p
        return np.concatenate([0.6 * f, -0.6 * (1 - f)])
    ```
  - Add guards so this cannot recur: `assert summary(N, ramp(p))["mean"] < TAB[N]["mean"]` for every (N, p); for N = 8, p = 1 compare the Nelder–Mead result with a coarse grid minimum over one period (γ period 2πσ because C is integer-valued, β period π/2); use several starts at N = 8 (ch13 used 20).
  - Re-run and replace every QAOA number in §1, §5 (5)–(7), exercise 3, §9, and the figure.
  - md:65 add one sentence, e.g. "비용을 줄이는 쪽으로 가려면 γ와 β의 부호가 반대여야 한다(13장의 MaxCut은 비용을 키우는 문제라 둘 다 양수였다)."

### A2. "p = 1은 ⟨C⟩를 줄이는데도 p\*가 무작위보다 낮다 … 최저 접힘의 진폭은 오히려 깎인다" is an artefact

- Location: md:198 ("첫째, …"), md:308 ("p = 1에서 ⟨C⟩는 내려갔는데 p\*는 무작위보다 낮아졌다"), md:342 ("이 장에서 p = 1은 평균을 줄이면서 p\*를 깎았고").
- Evidence: at the true p = 1 minimum p\*/random = 4.78 (N = 8), 6.36 (N = 10), 10.6 (N = 12), and P(valid) = 0.478 / 0.315 / 0.183 (not 0.190 / 0.121 / 0.066). The mechanism sentence ("한 층은 겹침을 조금 줄여 … 진폭은 오히려 깎인다") describes a local minimum, and was stated as fact without a supporting number.
- The underlying lesson (⟨C⟩-optimal ≠ p\*-optimal) is real and survives, with different numbers: p = 1, N = 8: 4.78× at the ⟨C⟩ minimum vs 5.65× at the p\* maximum; p = 3, N = 12: 78.8× vs ≥ 262.8×.
- Proposed replacement for md:198 "첫째": "첫째, ⟨C⟩-최적 각도는 p\*-최적 각도가 아니다. p = 1에서도 ⟨C⟩의 최솟점은 무작위의 4.8배·6.4배·10.6배이고 p\*를 직접 키운 각도는 5.7배·8.2배·…배다((나)). 틈은 p와 크기가 클수록 벌어진다." (numbers from the re-run). Delete the "진폭은 오히려 깎인다" sentence. md:308 and md:342 accordingly.

### A3. "훈련 비용을 0으로 쳐도 큰 크기에서는 담금질이 앞선다" and the bound "QAOA의 TTS₉₉는 24,724 이상" are wrong

- Location: md:15 (item 3), md:216–218 (column "QAOA 이상"), md:228 ("잔기 12개에서 뒤집힌다(24,724 대 16,045)"), md:230 ("훈련을 아무리 잘해도 (가)의 이상적인 값보다 좋을 수 없다고 보면 QAOA의 TTS₉₉는 24,724 이상이고 담금질(16,045)보다 길다"), out:81, out:83, assert py:512.
- Evidence: ⟨C⟩-optimal p = 3 at N = 12: ⟨C⟩ = 0.3863, P(valid) = 0.4975, p\* = 3.005×10⁻⁴ (78.8×), TTS₉₉ = 15,321. Re-checked with a second, structurally different mixer implementation (qubit-by-qubit) to 4 digits. The zero-training break-even against SA is p\* > ln(100)/16,045 = 2.87×10⁻⁴ = 75.2× random. 78.8× is above it; 15,321 vs 16,045 is inside the SA sampling error (15,351–16,795 at ±1 s.e., out:35; my own SA: 15,806 with ±2 s.e. 14,859–16,861).
- So with zero training cost the ⟨C⟩-optimal QAOA and the pre-fixed SA are tied at N = 12; QAOA is behind only the tuned SA (9,865). The "24,724 이상" argument fails twice: (가) was not the ⟨C⟩ minimum, and — as the md itself concedes — ⟨C⟩-optimal is not p\*-optimal, so (가) is not a bound in any case.
- Proposed replacement for md:15: "3. **훈련 비용을 0으로 치면 담금질과 비긴다.** ⟨C⟩를 정확히 최소화한 p = 3의 QAOA는 잔기 8·10개에서 담금질보다 5.3배·4.2배 짧고(215 대 1,144, 154 대 644), 12개에서는 15,321 대 16,045로 표본 오차 안이다. 일정을 맞춘 담금질(9,865)보다는 길다." (numbers from the re-run.) Delete the "24,724 이상" argument at md:230; see A5 for what should replace it.

### A4. "(나) 정답을 아는 각도 탐색" at N = 12 is 81.2× → at least 262.8×; three sentences built on it are false

- Location: md:200, md:216–218 (last column), md:226 ("어느 것도 90배에 닿지 않는다"), md:228 ("정답을 아는 각도 탐색(14,866)만이 … 비슷하고, 일정을 맞춘 담금질(9,865)보다는 길다"), md:308, md:342 ("절반쯤"), out:58–60, assert py:513.
- Evidence (`g_big.py`, `p_more12.py`, `k_transfer.py`): Nelder–Mead on −p\* from my ⟨C⟩-optimal point, 360 evaluations: p\* = 1.0024×10⁻³ = 262.8× random, TTS₉₉ = 4,592, angles (1.514, 2.840, 2.922 | −0.234, −0.205, −0.186); re-evaluated with the second mixer: 259.5× at the intermediate point, identical to 4 digits. N = 10, p = 3: 29.9× (TTS 105; chapter 26.0×, 121). N = 8, p = 3: 22.86× (135) — here the chapter's 22.9× is right. N = 8/10, p = 1 and 2: 5.65× / 12.74× and 8.2× / 17.6× (chapter out:58–59: 1.1× / 8.9× and 1.5× / 13.6×).
- Consequences: 262.8× > 90×, so "어느 것도 90배에 닿지 않는다" is false; 4,592 < 9,865, so "일정을 맞춘 담금질보다는 길다" is false; the ⟨C⟩-optimal p\* is 78.8/262.8 = 30% of the answer-aware value at N = 12 (63% at N = 8, 69% at N = 10), not "절반쯤".
- Note also that md:226 compares the two training-free values with the 90× threshold, which assumes 2,560 training shots. For training-free rows the threshold is 75×.
- Proposed fix: re-run (나) from the corrected (가) with a few hundred evaluations at N = 12 (the present 34 are far too few), and rewrite md:200, 226, 228 from the new numbers.

### A5. The verdict is not robust to the training protocol; parameter transfer meets the registered criterion as written

- Location: md:9 (verdict), md:14 (item 2), md:19, md:222–226, md:230 ("판정이 이 선택에 달려 있는지는 이렇게 확인했다"), md:343 ("다만 이 장에서는 그렇게 해도 담금질을 넘지 못했다"), md:198 ("셋째 … 21.8배 … 41.0배").
- Evidence (`k_transfer.py`, `l_more.py`; "0 evaluations at N = 12" means the angles were fixed using N = 8 and N = 10 only, where the optimum is known from the chapter's own enumeration):

  | p = 3 angles used unchanged at N = 12 | p\*/random | TTS₉₉ (0 training shots at N = 12) |
  | --- | --- | --- |
  | ⟨C⟩-optimal at N = 8 | 62.1 | 19,452 |
  | ⟨C⟩-optimal at N = 10 | 70.5 | 17,116 |
  | linear extrapolation 8 → 10 → 12 of the ⟨C⟩-optimal angles | 74.4 | 16,215 |
  | p\*-optimal at N = 8 | 161.8 | 7,458 |
  | p\*-optimal at N = 10 | 202.4 | 5,963 |
  | linear extrapolation of the p\*-optimal angles | 133.9 | 9,011 |

  p = 2: 29.3–36.2× (⟨C⟩-tuned, 33,328–41,131) and 55.4–69.0× (p\*-tuned, 17,502–21,801). SA: 16,045 (pre-fixed schedule), 9,865 (tuned; my own scan of ten schedules bottoms out at about 10,200).
- So: (i) with p\*-tuned p = 3 angles from the small instances, QAOA is shorter than SA and than the tuned SA with zero training shots at N = 12. The registered criterion is "QAOA의 TTS(훈련 샷 포함)가 담금질과 무작위 추출 둘 다보다 짧다"; nothing in ch19 says the angles must be trained on the N = 12 instance, and the chapter itself recommends transfer as "실용적인 방법" (md:343). Under that reading the criterion is met. (ii) With ⟨C⟩-tuned transferred angles it is a near-tie (16,215–19,452 vs 16,045). (iii) With SPSA at N = 12 it fails, sign fix or not (B1).
- Caveats on my side: single sequence; N = 8/10 are prefixes of the N = 12 sequence; I tried six transfer variants per p and report all of them; a stronger classical baseline (valid-fold moves, PERM) is not included; the unit still counts ≈ 3 × 10⁶–2 × 10⁷ two-qubit gates as 1.
- What this requires: not a specific verdict from me, but (a) removing the false robustness claim, (b) reporting the verdict per protocol, (c) an explicit author decision on whether classically pre-computed angles count as "훈련 샷 0". The ch19 prediction ("불통과(담금질에 진다)") is then not simply confirmed, and §1 cannot say "이 단위는 양자 쪽에 크게 유리한데도 졌다" without qualification — the unit argument becomes the main reason the result is not an advantage, not a side remark.
- Proposed replacement for md:230 from "판정이 이 선택에 달려 있는지는…" on: "판정은 훈련 방식에 달려 있다. 잔기 12개에서 샷으로 훈련하면(SPSA, 2,560샷) 시드에 따라 10만\~30만으로 담금질의 7\~19배다. 잔기 8·10개에서 정한 p = 3의 각도를 12개에 그대로 쓰면 12개에서의 훈련은 0샷이고, ⟨C⟩로 정한 각도는 16,215\~19,452로 담금질(16,045)과 비슷하며, 작은 크기의 정답으로 p\*를 키운 각도는 5,963\~9,011로 담금질보다 짧다. 등록한 기준의 글자대로면 마지막 방식은 통과다. 이 장은 등록할 때 훈련 방식을 정해 두지 않았다. 그것이 이 사전 등록의 구멍이다." (all numbers from the re-run.) md:343: delete "다만 이 장에서는 그렇게 해도 담금질을 넘지 못했다".

### A6. Exercise 4: "Perdomo-Ortiz et al. 2012는 아미노산 6개의 문제에서 이 방식으로 81큐비트를 썼다" misattributes the 81 qubits

- Location: md:336; §10 md:352 ("81큐비트 실험에서 바닥상태 13/10,000회. 본문을 확인했다. 등급 A").
- Source (PMC3417777, fetched): the Figure 1 caption says the 81 qubits are coloured so that "nodes with the same color represent the same logical qubit from the original 19-qubit Ising-like Hamiltonian". The 81 are *physical* qubits after embedding a 19-logical-variable problem in the hardware graph; the degree-reduction ancillas the paper reports are 5 (Experiment 1) and 1 (Experiment 3). The full six-residue problem did not fit ("Scheme 1 requires more resources than the number of qubits available"); the experiments are pieces of divide-and-conquer schemes. The model is Miyazawa–Jernigan, not HP. "13 out of 10,000" and the turn encoding are confirmed.
- Proposed replacement for the parenthesis: "(Perdomo-Ortiz et al. 2012는 아미노산 6개의 문제를 나눠 풀었고, 가장 큰 실험은 논리 변수 19개의 2차식을 장치의 연결에 맞춰 물리 큐비트 81개에 심은 것이다. 3차·5차 항을 2차로 낮추는 데는 보조 변수 1개·5개를 썼다)". In §10 add "(Miyazawa–Jernigan 모형; 81은 임베딩 뒤의 물리 큐비트 수)".

---

## B — should fix

### B1. (다) is one seed of a noise-dominated procedure; the start and gains are undisclosed

- Location: md:106, md:108, md:202–210, md:222, md:226; py:324–347.
- My own SPSA with the chapter's gains (a_k = 0.15/(k+3)^0.602, c_k = 0.15/(k+1)^0.101, gradient/σ, 64 shots × 2 × 20), different random streams (`h_spsa.py`, `h12.py`):

  | start | N = 8 (30 seeds) p\*/random median (min–max) | N = 10 (20 seeds) | N = 12 (3–4 seeds): p\*/random → TTS₉₉ |
  | --- | --- | --- | --- |
  | chapter `ramp(2)` | 4.09 (1.72–8.81) | 4.09 (2.17–9.42) | 5.2, 9.5, 4.0, 11.5 → 235k, 129k, 302k, 108k |
  | sign-corrected ramp, p = 2 | 5.63 (2.06–8.96) | 5.33 (3.03–9.25) | 6.7, 4.4, 9.9, 8.9 → 182k, 275k, 125k, 139k |
  | sign-corrected ramp, p = 3 | 5.56 (2.94–9.89) | 5.92 (2.15–12.44) | 9.2, 17.8, 7.8 → 133k, 70k, 157k |
  | N = 8 ⟨C⟩-optimal p = 3 angles, then SPSA | — | 17.6 (12.8–21.5); start was 19.3 | 53, 45, 30 → 25k, 29k, 42k; start was 62.1 (19,452) |

- Reading: the chapter's 7.38× / 5.65× / 6.23× are inside these spreads, so its SPSA code is consistent with mine. But "196,189", "12배", "6.2배" are single draws from a range of roughly 108k–302k (7–19× SA). With 2,560 shots SPSA mostly performs a random walk; started at good angles it makes them worse. The useful statement is "훈련하지 않고 옮겨 온 각도를 쓰는 편이 낫다", which the chapter cannot make because of A1.
- Not disclosed in the md: the SPSA gains, that the start has p\* below random, that only one seed was run, and that only p = 2 was trained at N = 12 although the hypothesis says p ≤ 3 (at N = 10 the chapter's own sensitivity row has p = 3 better than p = 2: 2,980 vs 3,126).
- Fix: report median and range over ≥ 10 seeds at N = 8/10 and ≥ 5 at N = 12; print the start-point ⟨C⟩ and p\* in the table; state the gains; run p = 1 and p = 3 at N = 12 as well or say why not.

### B2. md:14: "같은 회로에서 기댓값을 정확히 계산해 최적화하면 p = 3에서 48.8배" — not the same circuit

The trained circuit is p = 2; the comparable ideal value is the p = 2 one (chapter 25.6×, R13 35.5×). Proposed: "같은 p = 2 회로에서 기댓값을 정확히 최소화하면 …배, p = 3이면 …배다."

### B3. SA baseline: disclosures and small accounting points

- Best L is chosen after the fact from six lengths using the known answer (py:209–212; md:81 says so, md:230 does not list it among the chapter's choices). It is the Rønnow-style "optimal annealing time" and is standard, but it is uncharged tuning for SA while QAOA's training is charged. Say so in md:230.
- The md quotes 16,045 without its error; out:35 prints 15,351–16,795. Add "(성공률의 표준오차로 15,351\~16,795)". After the re-run the QAOA–SA gap is of this size, so the error matters.
- "성공" is "visited an optimum at any step ≤ L" (py:199–200), not "final state is optimal". Consistent with the QAOA and random definitions, but the md never says it. Add to md:174.
- A run of L steps looks at L + 1 bitstrings; the initial one is free (py:191). Negligible at L = 300/1000 (0.1–0.3%), 3% at L = 30. One clause or count L + 1.
- The move (py:194–196) always proposes a *different* direction (`+ integers(1, 4)`), so there are no null proposals; every proposal is one evaluation whether or not it is accepted. md:31 and md:230 describe it adequately; "한 걸음" there means a chain step while "담금질 1걸음" means an SA step (see C2).
- Seed collisions: `1000 * N + Lr` gives 11000 for (N = 8, L = 3000) and (N = 10, L = 1000), and 13000 for (N = 10, L = 3000) and (N = 12, L = 1000) — the verdict's SA number shares its seed with another run. Harmless here (different table sizes give different draws) but use `np.random.default_rng([N, Lr])`.
- md:238 gives SA = 600 for "본 실험" (N = 10, λ = 2) while md:179 gives 644 for the same setting (1,000 vs 2,000 runs, other seeds, other L set). Add "(1,000회 실행; (4)의 644와는 표본 오차 안)".

### B4. Pre-registration: undisclosed or unverifiable choices

md:230 lists λ, the schedule and move, SPSA (p = 2, 64 × 20) and σ normalisation. Missing: post hoc L (B3), the SPSA start and gains, single seed, p = 2 only at N = 12 (B1), and that (가)/(나) at N = 12 were budget-limited local searches from one start (partly said). "p, 샷 수, 반복 수는 실행 전에 정했고" (md:108) and "온도 일정은 실행 전에 정한 값이다" (md:182) cannot be checked from the repository (ch19 does not contain them; there is no dated record). Either point to a record or write "이 장을 쓰면서 정했고 결과를 보고 바꾸지 않았다".

### B5. §5 (7) λ sweep uses the poor local minimum as its start

- Location: md:232–242; py:404 starts from `IDEAL[(10, 2)]["x"]`.
- R13 (`i_lam.py`, best of four starts): λ = 0.5 / 1 / 2 / 5 / 20 → P(valid) 0.114 / 0.378 / 0.532 / 0.595 / 0.615, p\*/random 17.88 / 18.28 / 13.56 / 9.77 / 8.00 (chapter: 0.117 / 0.351 / 0.475 / 0.523 / 0.539 and 15.92 / 15.46 / 10.75 / 7.88 / 6.63). The trend the text draws survives; the numbers change.
- "비용의 차이가 거의 다 벌점에서 나오게 되어, 위상이 '겹쳤는가'만 구별하고 …" is a plausible reading, not something the script shows. Write "…때문으로 보인다" or print P(opt | valid) per λ.

### B6. "QAOA가 주로 한 일은 벌점을 피하는 것이다" holds only for ⟨C⟩-trained angles and is not shown

- Location: md:198 ("둘째"), md:310, md:340.
- Decomposition p\*/random = [P(valid) gain] × [P(opt | valid) gain] at N = 12, p = 3 (`l_more.py`): chapter angles 15.0 × 3.2; ⟨C⟩-optimal 17.3 × 4.5; p\*-optimal 11.5 × 22.6. So the sentence is supported for ⟨C⟩-trained QAOA and reversed for p\*-tuned angles. Print the two factors in the script and write "⟨C⟩로 훈련한 QAOA가 주로 한 일은…". The arrow chain "0.029 → 0.305 → 0.432" skips p = 1 (uniform → p = 2 → p = 3); label it.

### B7. md:172: "항이 2²⁰ − 1개로, 대칭이 허용하는 항이 거의 다 살아 있다"

The only symmetry the text names (flip all bits = reflection in the anti-diagonal; 00 ↔ 11, 01 ↔ 10 — correct) allows all even-degree terms: 2²¹ − 1 = 2,097,151, so half are alive. There is a second one: flipping the high bit of every step (00 ↔ 10, 01 ↔ 11) is the 180° rotation. Together they require an even number of Z on high bits and on low bits: 2¹⁰ · 2¹⁰ − 1 = 2²⁰ − 1. I checked that every non-zero term at N = 6, 8, 10, 12 obeys both parities and that C is invariant under all three flips (`j_sym.py`). At N = 12 (and N = 6) *all* allowed terms are alive; at N = 8 and 10 it is 43.5% and 43.7%.
Proposed: "홀수 차수는 없다. 모든 비트를 뒤집으면 접힘이 대각선에 대해 뒤집히고, 걸음마다 상위 비트만 뒤집으면 180° 돌 뿐이라 비용이 같기 때문이다. 두 대칭을 합치면 상위 비트와 하위 비트에 걸린 Z가 각각 짝수 개인 항만 남는다. 잔기 12개에서는 그런 항 2²⁰ − 1개가 전부 살아 있다(8·10개에서는 44%)." (print the 44% from the script.)

### B8. QUBO and gate-count wording

- ch19 (line 262) and ch13 (line 13) announce this chapter as the one that "QUBO로 옮긴다". It does not: the cost is a degree-20 polynomial applied as a table. "‘QUBO로 옮긴다’는 말의 비용을 숫자(항 100만 개)로 보여 줬다" (md:344) conflates two things — 1,048,575 terms is the cost of *not* reducing; the cost of reduction (ancilla count, extra penalties) was not computed. Say "이 장은 QUBO를 만들지 않았다. 차수를 낮추지 않으면 항이 100만 개이고, 낮추는 데 드는 보조 큐비트 수는 세지 않았다."
- md:17 and md:304 state "QAOA 한 층이 2큐비트 게이트 20,971,522개다" flatly; md:110 correctly calls it the simplest construction. Any 22-qubit diagonal unitary can be done with at most 2²² − 2 = 4,194,302 CX (Bullock & Markov 2004; from memory, not re-fetched), and Gray-code ordering of the Walsh terms brings it to about one CX per term. Add "가장 단순한 구성으로" in §1 and §7. The conclusion does not change.
- CLAUDE.md asks for gate costs via `tools/gate_cost.py` and a numerical check that the gate list equals the circuit. The script counts by formula only (it adds `tools` to `sys.path` and never uses it). My gate-by-gate emulation of the N = 6 case (CX ladder + rz + ladder for 255 terms, then rx) reproduces the diagonal-unitary distribution to 1.4×10⁻¹⁷ with 2,050 CX = Σ 2(k − 1); the script could do the same check on the case file.

### B9. Numbers in the md that the script does not print

All derived by hand from printed values; I recomputed each and none is wrong, but the house rule is that the script prints them: "75배" (md:182, 228; 1,207,215/16,045 = 75.2), "2.7배, 5.0배" (md:228; 2.746, 4.98), "2.9배" (md:210; 1,874/644 = 2.91), "rG = 55" (md:306; 18,388 × 3×10⁻³ = 55.2), "100배 넘게" (md:118; 18,388/170 = 108), "2,097만 개" (md:304), "97%" (md:310, 340), "6.2배 짧지만" for the TTS ratio (md:9, 226; out:79 prints 6.15), "절반쯤" (md:200, 308, 342), "전수 열거가 1초 안에" (md:312; the script prints only the table-building time, 1.07 s at N = 12; ch19 says 0.3 s for the enumeration), "ㄷ자 모양 8개 = 네 방향 × 좌우 두 가지" (count 8 is printed). Double-rounding hazards: 0.4315 → 0.432, 0.0655 → 0.066, 0.1165 → 0.117 are rounded from already-rounded 4-decimal prints; print 3 decimals or quote 4.

### B10. md:306 noise remark

"rG = 55다. 출력은 균등 분포와 구별되지 않는다" is a prediction of a depolarising-type model, not a result of this chapter, and for p layers it is 55p. By the book's own convention write it as "예측": "14장의 오류 모형으로는 한 층에서 이미 rG ≈ 55이고 무오류 확률이 e^{−55}이다. 출력이 균등 분포에 가까워진다고 예측할 수 있고, 그러면 …". (Amplitude damping would push toward 00…0, which here is the straight chain — valid, energy 0 — so "균등" is model-dependent.)

### B11. md:312 PERM / replica exchange

No source in §10. "잔기 100개 안팎의 2차원 문제를 푼다" should be "알려진 최저 에너지에 닿는다" (those results are best-known energies, not proofs of optimality). PERM is a chain-growth Monte Carlo with pruning and enrichment, not "가지치기 탐색". Candidates to cite (titles confirmed by search, details not re-fetched): Hsu, Mehlhorn, Nadler, Grassberger, "Growth algorithms for lattice heteropolymers at low temperatures", J. Chem. Phys. 118, 444 (2003), arXiv:cond-mat/0208042; Thachuk, Shmygelska, Hoos, "A replica exchange Monte Carlo algorithm for protein folding in the HP model", BMC Bioinformatics 8, 342 (2007). Likewise "유효한 접힘만 도는 고전 방법(pivot 이동)은 더 빠를 것이다" (md:310) is a guess; keep the "것이다" and add "이 장은 돌리지 않았다".

### B12. SDK snippets (not executed; no SDK here)

- CUDA-Q logic is right: my numpy emulation of the kernel exactly as written (H; per term CX ladder, `rz(2γw)` on the last qubit, ladder back; `rx(2β)`; standard rz/rx definitions; `(mask >> (n − 1 − q)) & 1`; `int(bits, 2)` with qubit 0 leftmost) gives ⟨C⟩ = 0.290386, p\* = 0.031291, identical to the diagonal version. 20,000-shot spread: ±0.0012 on p\*, ±0.010 on ⟨C⟩, so the comments "0.29 근처", "±0.001" are right. A wrong rz sign would give ⟨C⟩ = 3.69, p\* = 0.0019, so "첫 출력으로 바로 드러난다" is right.
- md:260 "Z 항 256개를 게이트로 풀어 넣는다": 255. The identity is skipped; 256 is "항등 포함" (out:110). The snippet's layer is 2,050 CX + 255 rz + 10 rx.
- md:298 fallback: "항마다 `cudaq.SpinOperator`를 만들어 `exp_pauli`로" is not how it works. The 0.16 docs pass `list[cudaq.pauli_word]` to the kernel and call `exp_pauli(coefficients[i], q, words[i])`, and the documented example reads it as e^{+iθP}. Proposed: "안 되면 항마다 길이 n의 문자열(`"ZIZI…"`)을 만들어 `list[cudaq.pauli_word]` 인자로 넘기고 `exp_pauli(-gamma * weights[t], q, words[t])`로 건다(`exp_pauli(θ, …)`는 e^{+iθP}라 부호가 반대다)."
- Kernel syntax: `list[int]`/`list[float]` arguments, `for`/`range`, local integer assignment and arithmetic are in the 0.16 kernel specification. I could not confirm from the docs `range(a, b - 1)` with run-time bounds, indexing an argument list with a computed index, or re-using the name `k` (a loop variable in the first inner loop, an assigned local in the second). I expect them to work; renaming the second `k` to `m` costs nothing. `x.ctrl`, `rz(angle, qubit)`, `cudaq.sample(..., shots_count=...)`, `.items()`, q0-leftmost bitstrings (the repo's `env/smoke_test.py` tests exactly this) and the default target `qpp-cpu` without a GPU are as documented.
- The stated reason for choosing CUDA-Q (GPU state vector beyond 18 qubits, md:250) cannot be exercised by this kernel: at N = 10 it needs 28,671 terms and 409,602 CX per layer, at N = 12 a 10⁷-entry argument list. Say that the snippet is for the 10-qubit case only.
- PennyLane: `qp.DiagonalQubitUnitary(D, wires)` (0.45.1 docs), RX(φ) = e^{−iφX/2}, `qp.probs` with wire 0 as the most significant bit — all as the md says. The table row "큰 크기: `lightning.qubit`·`lightning.gpu`" (md:257) is unverified for a 22-wire `DiagonalQubitUnitary`; if a device falls back to a dense matrix it cannot work. Mark it "확인하지 못했다".
- Status lines `# pennylane==0.45.x — 작성됨, 미실행` and `# cudaq==0.16.x — 작성됨, 미실행 (…)` follow CLAUDE.md and the repo's pins (`requirements/sdk.txt`: pennylane 0.45.1; `gpu.txt`: cudaq 0.16.0.post1).
- ch19 gave "기울기 계산이 편하고" as the reason for PennyLane; the snippet uses no gradient. One sentence, or drop the reference to the ch19 reason.

### B13. Regression asserts

They pin the table, λ, Z expansion, random, SA, the (가)/(나) ratios, (다) TTS and the verdict — i.e. they currently lock in the A1 artefacts. Not pinned: (다) ⟨C⟩ and p\*, the ten sensitivity rows, the λ sweep, the transfer values (21.8×, 41.0×), SA p_L. The N = 12 asserts depend on 15–37-evaluation Nelder–Mead paths through BLAS-backed `einsum`; a different BLAS may take another path and trip `atol = 0.02` (not verifiable here). After the fix pin converged values and add the two guards from A1. `rng = np.random.default_rng(21)` (py:33) is never used; "seed = 21" in the docstring and md:360 describes nothing (see C3).

---

## C — minor / wording

- **C1.** md:320 "문제 4는 N&C 연습문제 4.x가 아니라 12장의 Z 항 전개의 변형" reads like a leftover note. The rule is "교과서 변형 최소 1"; this chapter has a paper variant (Rønnow) and a variant of the book's own ch12 exercise, no N&C/KLM variant. ch19 handled the same situation by saying so. Proposed: "문제 2는 Rønnow et al.(2014)의 TTS 정의의 변형, 문제 4는 12장 연습문제 1(Z 항의 지수)의 변형이고 나머지 둘은 자작이다. QAOA와 TTS를 다루는 연습문제는 두 교과서에 없다." Or add a real N&C chapter 6 variant (Grover iterations for 16 marked items among 2²²: ≈ 402 oracle calls vs the TTS values) — the number must come from the script.
- **C2.** "걸음" means a chain step (md:27), an SA step (md:17, 96, 304) and both in one sentence (md:31, 230). Use "걸음" for the chain and "제안" or "스텝" for SA. "같은 기간에" (md:228) should be "같은 크기에서". N is the residue count here but "√N을 보장" (md:87) uses N for the search-space size (ch00's default meaning); write "탐색 공간 크기의 제곱근".
- **C3.** md:360 "seed = 21(담금질과 훈련은 따로 정한 seed)": nothing is drawn from seed 21. Write the actual seeds.
- **C4.** md:61: the product ∏_{l=1}^{p} does not say which layer acts first; add "(l = 1이 먼저 걸린다)". Side note outside this chapter: ch13 md:34 says "오른쪽 B가 먼저 적용" for e^{−iγC}e^{−iβB}, but ch13's code (and ch21) apply the cost layer first, so "13장의 정의 그대로다" is true of the code, not of that sentence.
- **C5.** md:230 "(가)의 p = 3은 37회": 35 optimiser evaluations plus two reporting evaluations (py:289, 296).
- **C6.** md:328: "227.95 … 227.9" are the same number printed with different precision (out:104).
- **C7.** Figure: panel (2) is cut at C = 25 (the maximum is 53) without a note; panel (4) has no bar for the answer-aware value or for "훈련 샷도 기회로". The six panel descriptions at md:244 match the png. Panel (1) checks: H at 0, 2, 5, 6, 8, 11 and five dotted contacts.
- **C8.** md:224 "세 방법 모두 예산 안에서 최저 접힘을 본다": for QAOA this is computed from p\*, not observed; training saw the optimum 0 times at N = 12. "TTS₉₉가 예산 안이다" is accurate.
- **C9.** md:302 Boulebnane: the abstract says deep circuits are *required* for accurate results and that QAOA's performance "can be matched by random sampling up to a small overhead". "깊은 회로에서도 … 수준" shifts the meaning slightly. Proposed: "정확한 결과에는 깊은 회로가 필요하고, QAOA의 성능은 작은 덧비용의 무작위 추출로 따라잡을 수 있다고 보고하고". Note for the author: my p = 3 numbers (hundreds of times random) are for this one small HP instance and are not in conflict with that paper.
- **C10.** md:310 "유효한 접힘만 도는 섞는 층 … (접힘에는 Fingerhuth et al. 2018)": the abstract confirms a split into hard and soft constraints, with the mixer keeping the hard ones. Whether the hard set includes full self-avoidance I could not confirm (my recollection is that it does not). Safer: "제약의 일부를 섞는 층이 지키게 하는 방법".
- **C11.** md:81 and `tts_from_p`: the rule "P ≥ 0.99이면 TTS = R" (which makes the 1.000 entries 3,000 and 10,000) is in the code only.
- **C12.** md:13 "무작위 추출과 같아진다": 2,984 vs 3,141 and 3,126 vs 3,208; "거의 같아진다".
- **C13.** md:23 defines contacts in words as "사슬에서는 떨어져 있는데 격자에서 이웃한", md:38 as j − i ≥ 3. Equivalent on the square lattice (j − i = 2 cannot be adjacent, by the parity statement at md:42, which is correct), but a newcomer will ask why 3; one clause.
- **C14.** md:104 calls (가) "훈련이 완벽했을 때의 값". Even after the fix it is the best *found* ⟨C⟩ minimum; write "찾은 최솟값".
- **C15.** md:162 "하한은 λ > 0, λ > 1, λ > 1": fine, but for N = 8 "λ > 0" means any positive λ works; λ = 0 does not. Say so in half a sentence, since the table below starts at 0.5.

---

## Independent recomputations

All with my own code in `scratchpad/r13/`. "✓" = agrees to the printed precision.

| Quantity | Chapter | R13 | Match |
| --- | --- | --- | --- |
| N = 8: valid folds / E_min / optimal bitstrings | 2,172 / −2 / 24 | same (table + my own DFS) | ✓ |
| N = 10 | 16,268 / −2 / 376 | same (table + DFS) | ✓ |
| N = 12 | 120,292 / −5 / 16 | same (table; no parity shortcut) | ✓ |
| max contacts among overlapping z (8/10/12) | 2 / 3 / 7 | 2 / 3 / 7 | ✓ |
| λ lower bound | 0 / 1 / 1 | 0.0 / 1.0 / 1.0 | ✓ |
| λ = 0.5, 1, 2 minima (N = 10; N = 12) | −2.5, 0/208; −2, 376/584; −2, 376/376; −6, 0/8; −5, 16/56; −5, 16/16 | same | ✓ |
| cost mean / σ at λ = 2 | 3.562/3.073; 5.252/4.220; 6.348/4.990 | 3.5625/3.0730; 5.2522/4.2197; 6.3477/4.9895 | ✓ |
| random TTS₉₉ | 3,141 / 3,208 / 1,207,215 | 3,141.5 / 3,208.4 / 1,207,215.4 | ✓ |
| Z terms / max degree / degree ≤ 2 / CX per layer / rotations, N = 8 | 1,782 / 12 / 38 / 18,388 / 1,796 | same | ✓ |
| N = 10 | 28,671 / 16 / 70 / 409,602 / 28,689 | same | ✓ |
| N = 12 | 1,048,575 / 20 / 110 / 20,971,522 / 1,048,597 | same | ✓ |
| degree histogram N = 8 | 38, 393, 804, 485, 60, 2 | same | ✓ |
| QAOA at the printed N = 8, p = 3 angles: ⟨C⟩ / P(valid) / p\* | 0.3272 / 0.7304 / 1.518e-2 | 0.32726 / 0.7301 / 1.5157e-2 (two mixers) | ✓ (3-decimal angles) |
| QAOA at the printed N = 12, p = 3 angles | 0.7939 / 0.4315 / 1.862e-4 | 0.7940 / 0.4311 / 1.8616e-4 | ✓ |
| chapter's N = 8 procedure, chapter's sign | 2.9582 / 0.6511 / 0.3272 | 2.9582 / 0.6511 / 0.3272, same angles | ✓ (reproduces the defect) |
| ⟨C⟩ minimum, N = 8, p = 1 / 2 / 3 | 2.958 / 0.651 / 0.327 | 1.2965 / 0.4740 / 0.1031 | ✗ (A1) |
| ⟨C⟩ minimum, N = 10 | 3.958 / 1.176 / 0.670 | 2.1218 / 0.9369 / 0.3578 | ✗ (A1) |
| ⟨C⟩ minimum, N = 12 | 4.731 / 1.511 / 0.794 | 2.6604 / 1.1746 / 0.3863 | ✗ (A1) |
| ideal p = 3 TTS₉₉ (훈련 0) | 301 / 209 / 24,724 | 215 / 154 / 15,321 | ✗ (A3) |
| answer-aware p = 3, p\*/random | 22.9 / 26.0 / 81.2 | 22.86 / 29.9 / ≥ 262.8 | N = 8 ✓, others ✗ (A4) |
| N = 10 angles used at N = 12, p = 2 / 3 | 21.8× / 41.0× | 33.4× / 70.5× (my ⟨C⟩-optimal angles) | ✗ (A1, A5) |
| SA N = 12, L = 1000: p_L → TTS₉₉ | 0.249 → 16,045 | 0.2527 ± 0.0069 → 15,806 (4,000 runs); 0.257 → 15,479 (3,000 runs) | ✓ within error |
| SA N = 12, L = 300 / 3000 | 0.079 → 16,899 / 0.528 → 18,376 | 0.085 → 15,553 / 0.530 → 18,298 | ✓ within error |
| SA N = 8, L = 300 | 0.701 → 1,144 | 0.686 ± 0.007 → 1,192 | ✓ (1.2 σ) |
| SA N = 10, L = 300 | 0.883 → 644 | 0.884 → 641 | ✓ |
| SA tuned, N = 12 | 9,865 (T = 0.5 fixed) | ≈ 10,200–10,600 (T = 0.4–0.5 fixed) | ✓ consistent |
| SPSA (다) p\*/random | 7.38 / 5.65 / 6.23 (one seed each) | medians 4.1 / 4.1, range 1.7–9.4; N = 12: 4.0–11.5 | consistent; not exactly reproducible (B1) |
| λ sweep p\*/random (N = 10, p = 2) | 15.92 / 15.46 / 10.75 / 7.88 / 6.63 | 17.88 / 18.28 / 13.56 / 9.77 / 8.00 | trend ✓, numbers ✗ (B5) |
| case HPHPPH: table, 8 optima, σ, 256 terms incl. identity | json | identical; Σ w_S(−1)^{S·z} = C exactly | ✓ |
| case: ⟨C⟩ / p\* at γ = −0.8199, β = 0.2145 | 0.2904 / 0.0313 | 0.290386 / 0.031291 (diagonal and gate-by-gate) | ✓ |
| case: is that the p = 1 minimum? | — | yes; grid minimum 0.2907 at the conjugate point (+0.818, −0.218) | ✓ |
| exercise 1 | 64 / 36 / 8 / 34.5 | 64 / 36 / 8 / 34.5 (28 overlapping) | ✓ |
| exercise 2 | 664.4 / 2,000.0 / 227.9 | 664.4 / 2,000.0 / 227.95 | ✓ |
| exercise 3 | 2,840 / 2,999 / 1,182,491 | arithmetic ✓; inputs are the A1 values | ✗ inputs |
| exercise 4 | 22 CX; 1,744 of 1,782 | 22; 1,744 | ✓ |

Not recomputed: "훈련 샷도 기회로 센 TTS₉₉" (1,103 / 1,365 / 195,176), "훈련 중 최저 접힘을 본 샷" (19 / 9 / 0) and the ten sensitivity rows — they depend on the chapter's random streams. I read the code for them (py:330–347) and found no accounting error: every training shot is counted (2 × 64 × 20 = 2,560), the final exact evaluation is reporting only, and the conditional expression at py:346 parses as intended.

## Verified OK

- Every table entry and in-text number that the script prints matches `ch21_stdout.txt` (§1, §2, §4 gate table, §5 (1)–(7), §8); the only numbers not from the script are those in B9.
- HP model definition; parity statement (even chain distance for overlaps, odd for contacts); "flip all bits = reflection" (00 ↔ 11, 01 ↔ 10); only even-degree Z terms; HPPH hand example.
- TTS₉₉ formula and its use agree with ch00 and ch19; the three methods use the same success event (optimum seen at least once) and the same unit. The extra column "훈련 샷도 기회로" is the more favourable accounting for QAOA, as the md says.
- Mixer and phase conventions: md formula, `qaoa_probs`, PennyLane `RX(2β)`, CUDA-Q `rx(2β)` / `rz(2γw)` are mutually consistent (checked by gate-level emulation).
- SA implementation: Metropolis rule, geometric schedule, vectorisation, first-hit bookkeeping are correct; no overflow in `exp`.
- Pre-registered items reported as registered: sizes 8/10/12, unit, baselines, exhaustive answer, 10⁷ budget, auxiliary metric, N = 8/10 as reference only; tools PennyLane (main) and CUDA-Q (other) with a reason sentence.
- Format: ten sections in order; `\~ \_ \*` escaping; no paragraph starting with `|`; table column counts consistent; "실측" not used; snippet status lines; versions agree with `requirements/`.
- Citations confirmed against a fetched page: Boulebnane et al. (arXiv 2204.01821 abstract; DOAJ: npj Quantum Information 9, 2023, DOI 10.1038/s41534-023-00733-5 — article number 70 not shown there); Perdomo-Ortiz et al., Sci. Rep. 2, 571 (2012) (but see A6); Fingerhuth, Babej, Ing, arXiv:1810.13411 (2018); Robert et al., npj Quantum Information 7 (2021), 10 residues on 22 qubits and 7 on 9 (article number 38 not shown on DOAJ); Babbush et al., Adv. Chem. Phys. 155, 201–243 (2014); Barkoutsos et al., Quantum 4, 256 (2020) (journal URL). Title-level search hits only, bibliographic details from memory and consistent with the md: Hadfield et al., Algorithms 12, 34 (2019); Spall, IEEE TAC 37, 332 (1992); Rønnow et al., Science 345, 420 (2014); Kirkpatrick et al., Science 220, 671 (1983); Lau & Dill, Macromolecules 22, 3986 (1989); Crescenzi et al., J. Comput. Biol. 5, 423 (1998); Berger & Leighton, J. Comput. Biol. 5, 27 (1998); Farhi et al., arXiv:1411.4028. N&C/KLM having no QAOA: correct.

## Could not verify

- Execution of either SDK snippet (no SDK, no installs). CUDA-Q points open: run-time `range(a, b − 1)`, computed list indices, name re-use; PennyLane `lightning.*` with a 22-wire diagonal.
- Whether the SA schedule and the SPSA settings were fixed before running (B4).
- Global optimality of my own QAOA minima for N ≥ 10 (local searches in the annealing-like basin; N = 8 is backed by a grid for p = 1 and 58–62 starts for p = 1–3).
- Article numbers 70 and 38 for the two npj Quantum Information papers, and the hard-constraint set in Fingerhuth et al. (C10).

## Scratch files (`scratchpad/r13/`)

`a_tables.py` (tables, λ, FWHT), `qaoa_lib.py` / `fast.py` / `fast2.py` (three mixer implementations), `b_qaoa8.py` (chapter angles, p = 1 grid), `c_multi.py` (multi-start N = 8), `o_fix.py` (chapter procedure with both signs), `g_big.py` / `p_more12.py` (N = 10, 12 refinement and answer-aware search), `k_transfer.py` / `l_more.py` (transfer, decomposition), `h_spsa.py` / `h12.py` (SPSA), `f_sa.py` / `m_sa_tune.py` (SA), `i_lam.py` (λ sweep), `j_sym.py` (symmetries), `e_case.py` (case file, CUDA-Q emulation), `n_md.py` (md format); outputs `*.out`, `multi_N8.json`, `big_N10.json`, `big_N12.json`.
