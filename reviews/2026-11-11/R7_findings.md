# R7 — ch19 "후보 문제 선정" independent review

Reviewed: `/home/claude/qc-book/ch19/ch19_candidates.md` (line numbers below refer to this file), `/home/claude/qc-book/ch19/ch19_candidates.py`, against ch17/ch18, the plan's candidate table, and primary sources.
Scratch code and outputs: `/tmp/claude-0/-home-claude/fd63146d-511d-5ae9-8dab-3100fb85f0db/scratchpad/mid/R7/` (`rerun.txt`, `kern_indep.py`, `kern_exp.py`, `bw_check.py`, `heis_indep.py`, `hp_indep.py`, `md_cell.py`, `schwinger.py`, `mod/ch19_mod.py`).

Rerun result: every seeded number in the script stdout reproduces bit-for-bit. Only timings moved (FFT 9.0 → 9.1 µs, HP N = 12 0.25 → 0.21 s, N = 30 extrapolation 112 → 91 days).

Counts: **A 2, B 11, C 14.**

---

## A — factually or numerically wrong

### A1. ch19 §5(4) — md:239 (propagates to md:16, 250, 290, 308; script:206–208)
Quote: "χ는 2^S 꼴로 자란다. MPS가 버거워지는 기준을 χ ≈ 10⁴(S ≈ 13.3비트)로 잡고 같은 포화 비율을 쓰면 n ≳ 42, t ≈ 4.6이 경계다."

What is wrong: the text equates χ = 10⁴ with S = log₂10⁴ = 13.3 bits, i.e. it uses χ = 2^S. But χ is defined in §4 (md:138) as "discarded weight ≤ 10⁻⁶", and the table directly above (md:234–237) shows that this χ is 4–7 times 2^S: t = 1.0: 2^2.98 = 7.9 vs χ = 43; t = 1.5: 2^4.30 = 19.7 vs χ = 133; t = 2.0: 2^5.03 = 32.7 vs χ = 214. ch17 itself warns that 2^S is only a lower bound (ch17 md:260: "2^S ≈ 14, χ = 180").
Evidence (`heis_indep.py`, independent Hamiltonian construction, same numbers as the script): χ/2^S = 3.6, 4.2, 4.6, 5.5, 6.1, 6.8, 6.7, 6.5 for t = 0.25 … 2.0. Fit of the chapter's own χ column on t = 0.25–1.5: log₂χ = 1.76 + 3.58 t (R² = 0.999). With that, χ = 10⁴ is reached at t ≈ 3.2 where S ≈ 9.4 bits, giving n ≳ 2·⌈9.4/0.646⌉ = 30, not 42. The boundary circuit is then n = 30, 33 steps, 3·(15·34 + 14·33) = 2,916 CNOT, (1 − 0.003)^2916 = 1.6×10⁻⁴ instead of 5,721 CNOT and 3.4×10⁻⁸.
The qualitative conclusion (dead without mitigation) survives, but the headline numbers in the summary, the verdict table and the pre-registration (c) are derived from an identity the chapter's own table contradicts by a factor of ~6.

Proposed replacement (md:239, from "긴 사슬에서는" to the end of the paragraph):

> 긴 사슬에서는 선형 증가가 계속되고 χ는 시간에 지수로 자란다. 다만 χ = 2^S는 하한일 뿐이다. 위 표에서 버린 가중치 10⁻⁶ 기준의 χ는 2^S의 4\~7배다(t = 2에서 2^5.03 ≈ 33인데 χ = 214). MPS가 버거워지는 기준을 χ ≈ 10⁴으로 잡으면 경계를 두 가지로 어림할 수 있다. 하한 χ = 2^S를 쓰면 S ≈ 13.3비트가 필요하고, 같은 포화 비율(65%)에서 n ≳ 42, t ≈ 4.6이다. 표의 χ를 그대로 외삽하면(t = 0.25\~1.5에서 log₂χ ≈ 1.76 + 3.58t) t ≈ 3.2, S ≈ 9.4비트, n ≳ 30에서 이미 χ ≈ 10⁴에 닿는다. 아래에서는 고전 쪽에 유리한 앞의 값(n ≈ 42, t ≈ 4.6)을 경계로 쓴다. 두 값 모두 n = 16 하나에서 얻은 거친 외삽이고, 23a장이 여러 크기로 다시 잰다.

And append to md:250 after "…= 3.4×10⁻⁸이다.":

> 표의 χ로 외삽한 경계(n = 30, t = 3.2, 33단계)로 계산하면 CNOT 2,916개, 무오류 확률 1.6×10⁻⁴이다. 어느 쪽이든 결론은 같다.

(If the second estimate is added, also add the fit and the n = 30 line to the script so the numbers have a source. In md:308 (c) change "n ≈ 42, t ≈ 4.6" to "n ≈ 30\~42, t ≈ 3.2\~4.6".)

### A2. ch19 §7 — md:342
Quote: "1번과 2번은 뒤의 것이다. … 판정표의 "탈락"은 뒤의 경우에만 썼다."

What is wrong: the verdict table (md:287–296) also gives "탈락" to candidates 6, 9 and 10. Candidate 10 (option pricing) does not lose in scaling — it has a proven quadratic speed-up and is dropped because of the constant (break-even 33 days–9,200 years, md:281, ch17 md:182). Candidate 9 (Max-Cut QAOA) is dropped because a p = 1 guarantee is below a classical guarantee, which is not a scaling statement either. So the sentence is false for two of the five "탈락" rows.

Proposed replacement for the last sentence of md:342:

> 판정표의 1번·2번·6번 "탈락"은 뒤의 경우다. 9번과 10번의 "탈락"은 이유가 다르다. 9번은 증명된 가속이 없고 고전 보장(0.878)에 못 미치며, 10번은 제곱 가속이 있지만 호출 하나의 비용 때문에 손익분기가 33일\~9,200년으로 밀린다.

---

## B — misleading / incomplete / contradicts another chapter

### B1. ch19 §1, §5(2) — md:14, md:206 (script:106–118)
Quote (md:206): "둘째, 그 결과 같은 화자와 다른 화자의 구분이 사라진다. 1-NN 정확도 0.31은 화자 4명에서 찍어 맞히는 수준(약 0.24)에 가깝다. 특성을 더 주었는데 성능이 떨어졌다."

What is wrong: with 4 speakers the centred centroids span 3 dimensions, so PCA components 4–10 carry no speaker information at all, and the script stretches each of them to the full [0, π] range before using it as a rotation angle. Evidence (`kern_indep.py`): between-class variance fraction per PC = 0.984, 0.973, 0.960, 0.002, 0.001, 0.000, 0.000, 0.001, 0.001, 0.001; explained variance per PC = 0.146, 0.114, 0.109, then ≈ 0.017 each. This is why the quantum 1-NN peaks at n = 3 (0.98) and falls from n = 4 on — the summary's "0.91에서 0.31로 떨어진다" hides that non-monotonic step. The collapse to chance is therefore mostly "pure-noise features fed at full bandwidth", not concentration as such.
Control experiment (`kern_exp.py`, same generator, noise 0.085, 11 speakers so that all 10 PCs are informative): at n = 10 the kernel is exactly as concentrated (mean 0.0021, std 0.0022, alignment +0.200 vs RBF +0.967) but quantum 1-NN is 0.77 against chance 0.088. With 8 speakers (7 informative PCs) n = 10 gives 0.45 against chance 0.12. What robustly degrades with n is the kernel magnitude (shot cost) and the alignment; "1-NN falls to chance" is specific to the 4-speaker design. The RBF comparison itself is fair (same features, RBF stays at 1.00), so the verdict stands, but the causal sentence does not.

Proposed replacement for the three quoted sentences (md:206):

> 둘째, 같은 화자와 다른 화자의 구분이 흐려진다. 정렬이 +0.646에서 +0.174로 줄고, 1-NN 정확도 0.31은 화자 4명에서 찍어 맞히는 수준(약 0.24)에 가깝다. 다만 1-NN의 붕괴에는 이 데이터의 사정이 섞여 있다. 화자가 4명이면 화자 중심의 차이는 3차원에 다 들어가므로 넷째 주성분부터는 화자 정보가 없는 잡음이다(주성분별 화자 간 분산 비율이 1\~3번째는 0.96\~0.98, 4\~10번째는 0.002 이하). 그 잡음 성분도 \[0, π\] 전체로 늘려 회전 각도로 썼다. 1-NN이 n = 3에서 0.98로 가장 높고 n ≥ 4에서 떨어지는 이유다. 같은 잡음 특성을 받은 RBF 커널은 1.00을 유지했으므로 "잡음 특성에 약하다"는 것은 양자 커널 쪽의 약점이 맞지만, n과 함께 반드시 나빠지는 것은 커널 값의 크기와 정렬이지 1-NN이 아니다.

Summary (md:14), replace "화자를 가르는 1-NN 정확도가 0.91에서 0.31로 떨어진다." with:

> 화자를 가르는 1-NN 정확도가 0.91에서 0.31로 떨어진다(n = 3의 0.98이 가장 높고, 화자 정보가 없는 넷째 주성분부터 넣으면 떨어진다).

(Optional: add the 8- or 11-speaker control to the script; ch22's pre-registration already says "화자 4\~8명".)

### B2. ch19 §5(4) — md:232, md:250, md:290 vs ch17 md:192
Quote (md:232): "18장의 자유 페르미온 검사를 통과하지 못한다는 뜻이고, 2n×2n 행렬 풀이가 없다." and (md:250) "판정은 17장과 같은 "가능"이되 조건이 붙는다."

What is wrong: ch17's "가능" for spin-chain quenches is explicitly conditioned on non-integrability: "비적분 모형에서 χ가 t에 지수(A). … 가능 — 23a장(비적분·2D 모형, 얽힘이 자라는 초기 상태)". The nearest-neighbour spin-1/2 Heisenberg chain is not free-fermion, but it is Bethe-ansatz integrable, and ch19 never says so (no occurrence of "적분" or "Bethe" in the chapter). Integrability does not stop entanglement growth (the table shows it), but it gives classical handles that a "이득 후보" must disclose: the Néel-quench steady state is solved exactly (Wouters et al., PRL 113, 117202 (2014), https://arxiv.org/abs/1405.0172; Pozsgay et al., PRL 113, 117203 (2014)), and the entanglement growth after the Néel quench in XXZ is predicted by the quasiparticle picture (Alba & Calabrese, PNAS 114, 7947 (2017), https://arxiv.org/abs/1608.00614). So "판정은 17장과 같은" is not accurate: ch19 applies ch17's verdict to a model ch17's condition excludes.

Proposed insertion after "2n×2n 행렬 풀이가 없다." (md:232):

> 다만 이 사슬은 자유 페르미온은 아니어도 Bethe 가설풀이(Bethe ansatz)로 풀리는 적분 가능 모형이다. 17장의 "가능" 판정은 비적분 모형을 조건으로 달았다. 적분 가능성은 얽힘 증가를 막지 못하지만(아래 표) 고전 쪽에 다른 길을 준다. Néel 급랭의 정상 상태는 정확히 풀렸고(Wouters et al. 2014), 얽힘 엔트로피의 증가도 준입자 그림으로 예측된다(Alba & Calabrese 2017). 유한 시간의 일반 상관 함수에는 여전히 수치 계산이 필요하다. 그래서 23a장은 이 사슬을 기준 사례로 쓰되, 다음 이웃 결합 J₂를 더해 적분 가능성을 깬 사슬을 함께 돌린다.

md:250, replace "판정은 17장과 같은 "가능"이되 조건이 붙는다." with:

> 판정은 "가능"이되 조건이 둘 붙는다. 하나는 위의 적분 가능성이고(17장의 "가능"은 비적분 모형이 조건이다), 다른 하나는 잡음이다.

and add to §10: Wouters, De Nardis, Brockmann, Fioretto, Rigol, Caux, *Phys. Rev. Lett.* 113, 117202 (2014); Alba & Calabrese, *PNAS* 114, 7947 (2017).

### B3. ch19 §5(7) — md:263, md:293
Quote: "가장 유리한 축에서도 4차이므로 상호작용 모형이고 4번과 같은 판정을 받는다."

What is wrong: §3 (md:52) says the deciding numbers for a dynamics candidate are two — Majorana degree and entanglement growth rate. For candidate 7 only the first was computed, yet it receives "가능(4번과 같은 조건)". The gauge-eliminated Schwinger model is the textbook case where the second number can fail: confinement suppresses entanglement growth, and MPS real-time simulation works well there (Buyens et al., PRL 113, 091601 (2014), https://arxiv.org/abs/1312.6654; ch17 already cites Kormos et al. 2017 on confinement).
Evidence (`schwinger.py`, 12 sites, H = Σ(XX+YY) + m Σ(−1)ⁿZ + J ΣL², quench from the bare vacuum): with J = 0 (free hopping) the half-chain entropy reaches 4.58 bits at t = 2; with J = 2, m = 1 it never exceeds 0.91 bits up to t = 3 (χ ≤ 32); with J = 4 it stays below 0.65 bits. So in a wide parameter range candidate 7 is caught by ch18's "낮은 얽힘" row and would not get candidate 4's verdict.
Also (wording) "15쌍 중 10쌍" reads as if 10 non-neighbour pairs carry ZZ; the 10 ZZ pairs are all pairs among sites 0–4, of which 4 are nearest neighbours and 6 are not.

Proposed replacement for md:263 from "가장 유리한 축에서도" to the end:

> 가장 유리한 축에서도 4차이므로 자유 페르미온은 아니다. 그러나 3절의 표대로 동역학 후보의 결정적 숫자는 둘이고, 둘째(얽힘 증가 속도)는 이 후보에서 재지 않았다. Schwinger 모형은 가둠(confinement) 때문에 매개변수에 따라 얽힘이 느리게 자라고, 그런 영역은 MPS로 긴 시간까지 계산된다(Buyens et al. 2014). 판정은 "가능(얽힘이 빠르게 자라는 매개변수 영역에 한함)"으로 적는다. 4번과 다른 점은 ZZ가 이웃이 아닌 쌍에도 걸려(ZZ 10쌍 가운데 6쌍) 회로가 깊어진다는 것이다. 따로 장을 두지 않고, 얽힘 증가를 재는 일까지 묶어 23a장의 연습문제로 돌린다.

Verdict table row 7 (md:293), cells 4–5: "Majorana 4차, 장거리 ZZ 10개(얽힘 증가는 미측정)" / "가능(얽힘이 빨리 자라는 영역에 한함)".

### B4. ch19 §3, §5(2), verdict table — md:50, md:218, md:288
Quote (md:218): "그때 양자 커널은 RBF 커널과 상관 0.74로 닮아 있고 … 18장 휴리스틱 트랙의 4′(측정 예산)와 5′(고전 대리 계산)에 차례로 걸린 것이다."

What is wrong, three points.
(i) The decisive number cited in the verdict table ("n = 10에서 1-NN 0.31(RBF 1.00)") is an exact-kernel performance comparison with a classical kernel on the same features. In ch18's track (ch18 md:139–142) that is 6′ (실측 비교: "같은 예산에서 고전이 같거나 나음"), not 4′ or 5′. §3's table (md:50) makes the same mis-assignment ("성능 차" → 4′, 5′).
(ii) ch18's 4′ fail signal is "예산이 17장의 손익분기를 넘음". 5,058 shots per entry at n = 10 (2.5×10⁷ shots for a 100-utterance matrix) is not shown to exceed anything; what fails is the scaling s ∝ 2ⁿ. Say so.
(iii) The 5′ claim rests on a Pearson correlation of 0.74 (r² ≈ 0.55), and the same chapter's pre-registration (md:306) uses 0.9 as the line between "resembles RBF" and "does not". By the chapter's own threshold the c = 0.05 kernel does not count as RBF-like. In addition the similarity is not monotonic in bandwidth: `bw_check.py` gives 0.742 (c = 0.10), 0.739 (c = 0.05), 0.616 (c = 0.02), so "대역폭을 줄이면 … 그만큼 RBF 커널과 닮아 간다" (md:14) holds only down to c ≈ 0.1. (Minor: at c = 1.0 the best-matching RBF width sits on the edge of the search grid, g = 100, script:126.)

Proposed replacement for the last three sentences of md:218:

> 성능이 돌아오는 구간에서 양자 커널은 고전 커널 쪽으로 다가간다. 다만 상관 0.74는 양자 커널 값의 분산 가운데 절반쯤(0.74² ≈ 0.55)을 RBF 하나가 설명한다는 뜻이고, 22장의 사전 등록이 "닮았다"의 문턱으로 둔 0.9에는 못 미친다. c를 0.02까지 더 줄이면 상관은 0.62로 도로 내려간다. 이 데이터에서 시험한 설정 가운데 양자 커널이 RBF보다 나은 것은 없었다(c = 0.05의 1-NN만 1.00으로 같다). 18장 휴리스틱 트랙으로 말하면 4′(측정 예산: 샷 수가 2ⁿ에 비례해 는다)와 6′(실측 비교: 같은 특성의 고전 커널이 같거나 낫다)에 걸렸고, 5′(고전 대리 계산)는 방향만 보였다.

md:50 third cell: "18장 4′(측정 예산), 5′(고전 대리 계산), 6′(실측 비교)". md:288 third cell: "18장 4′·6′(5′는 부분적)". md:14 last sentence: "대역폭을 줄이면 정확도가 돌아오지만 그 구간에서 RBF 커널과의 상관이 0.08에서 0.74까지 오른다."

### B5. ch19 pre-registration, ch22 row — md:306
Quote: "NMI가 RBF보다 0.02 이상 높고, 그 설정에서 가장 닮은 RBF 커널과의 상관이 0.9 미만"

What is wrong:
(i) Ceiling: on the chapter's data the RBF baseline already has NMI = 1.000 at every n (`kern_exp.py`, spectral clustering on the RBF kernel, n = 2, 3, 4, 8, 10; cosine also 1.000). "0.02 higher than RBF" is then impossible by construction, so the test fails trivially. §7 (md:330) promises "더 어려운 설정" but the pre-registration does not fix the difficulty.
(ii) Resolution: one misassigned utterance changes NMI by 0.030 at 4 × 25 utterances (1 error → 0.970) and by 0.010 at 8 × 25. A 0.02 threshold on a single seed is below (or at) one-utterance resolution; no seeds or error bars are specified.
(iii) "그 설정": the set of settings (n, bandwidth c, speakers, shots) is not enumerated, so a pass can be found by searching settings after the fact — the thing pre-registration is meant to prevent.
(iv) The correlation clause never binds: the largest correlation the chapter measured is 0.742.

Proposed replacement cells for the ch22 row (크기와 예산 / 통과 기준):

> 화자 8명 × 발화 25개(4명은 참고용), n ∈ {2, 4, 6, 8, 10}, 대역폭 c ∈ {1, 0.5, 0.25, 0.1, 0.05}, 무잡음 커널과 1,000샷 커널, 데이터 시드 10개. 임베딩 잡음은 RBF 커널의 NMI가 0.8\~0.9가 되도록 먼저 맞추고 고정한다(이 장의 데이터는 RBF의 NMI가 1.000이라 넘을 여지가 없다)

> (n, c) 설정 하나 이상에서 시드 평균 NMI가 RBF보다 0.02 이상 높고, 그 차이가 시드 간 표준오차의 2배를 넘으며, 그 설정에서 가장 닮은 RBF 커널과의 상관이 0.9 미만. 설정은 시드 절반으로 고르고 나머지 절반으로 판정한다

### B6. ch19 pre-registration, ch21 row — md:305
Quote: "잔기 8개 이상에서 QAOA의 TTS가 담금질보다 짧다" (예산 10⁶ 평가, 단위 "회로 실행 1샷 = 담금질 1걸음")

What is wrong: at the registered sizes the instances are so easy that blind guessing decides the outcome. From the chapter's own table (md:256–258) the fraction of bitstrings that are optimal folds is 8/4⁵, 24/4⁷, 376/4⁹; uniform random bitstrings therefore have TTS₉₉ = 587 (N = 6), 3,141 (N = 8), 1,794 (N = 9), 3,208 (N = 10) evaluations (`bw_check.py`; N = 9 from the script stdout, 168 optimal folds). Any QAOA whose training uses more than ~3×10³ shots (a few cost evaluations at 1,000 shots each) loses to random sampling regardless of how good its final distribution is. The test measures budget accounting, not the algorithm; it is a trivially-failing criterion and it does not list the random-sampling baseline that makes this visible. Two smaller gaps: "잔기 8개 이상" does not say whether all of N = 8, 9, 10 must pass, and nothing says what TTS is when the optimum is never sampled within the budget.

Proposed replacement cells (지표와 고전 베이스라인 / 통과 기준 / 예상):

> TTS₉₉, 단위는 목적함수 평가 1회(회로 실행 1샷 = 담금질 1걸음 = 무작위 비트열 1개). 베이스라인은 담금질과 균등 무작위 비트열 추출, 정답은 전수 열거. 보조 지표: 훈련 비용을 뺀 한 샷 성공 확률 P와, 그 P가 무작위 추출의 몇 배인가

> 잔기 8·9·10개 모두에서 QAOA의 TTS(훈련 샷 포함)가 담금질과 무작위 추출 둘 다보다 짧다. 예산 안에 최적해를 한 번도 못 보면 TTS = ∞로 적는다

> 불통과. 이 크기에서는 무작위 추출의 TTS₉₉가 587(N = 6)\~3,208(N = 10)이라, 훈련에 수천 샷 이상을 쓰는 QAOA는 최종 분포의 품질과 무관하게 진다. 그래서 보조 지표(배율)를 함께 보고한다

### B7. ch19 pre-registration, ch20 row — md:304 vs ch17 md:195
Quote: "잘 조건화된 희소 행렬에서 QSVT 역행렬의 블록 인코딩 호출 수는 κ에 거의 선형이고 크기 N과 무관하다 … 차수의 κ 지수가 1 ± 0.2이고, κ를 고정하고 N을 16배 키울 때 차수 변화가 10% 미만"

What is wrong:
(i) The second criterion is true by construction: the degree of the 1/x polynomial on [1/κ, 1] depends only on κ and ε, never on N. It cannot fail, so it is a sanity check, not a hypothesis.
(ii) "호출 수" is identified with the polynomial degree, which leaves out the success probability. ch17's own candidate table says "QSVT 차수 κ·log(κ/ε), 진폭 증폭 포함 총 κ²·log(κ/ε)". As written, the hypothesis ("호출 수는 κ에 거의 선형") contradicts ch17 unless it says "degree only".
(iii) No κ grid is given. The only well-conditioned family in this chapter has κ = 2.62–2.79, which cannot determine an exponent; the Poisson families at N = 16–256 have κ = 117–2.7×10⁴ (1D), i.e. degrees of order 10⁴–10⁶.
(iv) κ must be the effective condition number of the block encoding (subnormalisation α/σ_min), not κ(A).

Proposed replacement cells (가설 / 크기와 예산 / 통과 기준):

> 잘 조건화된 희소 행렬에서 QSVT 역행렬 다항식의 차수(진폭 증폭을 뺀 블록 인코딩 호출 수)는 블록 인코딩의 유효 조건수 κ에 거의 선형이고, 성공 확률을 포함한 총 호출 수는 κ²에 가깝다(17장)

> N = 16\~256(4\~8큐비트와 보조 큐비트), κ는 I + γS 족의 γ를 바꿔 {2, 4, 8, 16, 32}로 맞춘다. 상태 벡터 시뮬레이션

> 차수의 κ 지수가 1 ± 0.2이고, 성공 확률을 포함한 총 호출 수의 κ 지수가 2 ± 0.3. 차수가 N과 무관한지는 검산 항목으로만 보고한다(다항식이 κ와 ε만으로 정해지므로 가설이 아니다)

### B8. ch19 pre-registration, ch23a row — md:308
Quote: "(b) 결정계수 0.95 이상 (c) 추정치가 이 장의 n ≈ 42, t ≈ 4.6과 2배 안에서 일치"

What is wrong:
(b) "포화 전까지" is not defined, so the fit window is a free choice; and R² ≥ 0.95 of log χ against t does not separate exponential from polynomial growth on a short window. A pure power law χ ∝ t^p gives R² = 0.974 on t ∈ [0.5, 1.5] and 0.936–0.941 on [0.25, 1.5] and [0.5, 3.0], for any p (`heis_indep.py`). The chapter's own data gives R² = 0.999 on [0.25, 1.5] but 0.975 on [0.25, 2.0] and 0.84 on [0.25, 3.0], so pass or fail is decided by the window.
(c) is not a well-posed test. The hypothesis is "the circuit at the MPS limit loses the signal", the metric is "MPS 한계와 잡음 한계의 교차점 추정", and the criterion compares against (n ≈ 42, t ≈ 4.6), which is the MPS limit alone, not a crossing; "잡음 한계" has no definition; "2배 안" of a pair (n, t) is undefined; and with n = 10–20, t ≤ 3 ch23a can only repeat the same extrapolation, so agreement is self-referential. See also A1: the reference value itself moves to n ≈ 30, t ≈ 3.2 under the chapter's χ definition.
(a) passes on data already in this chapter (ratios 16.2 and 16.2 for 10→20→40 steps; fitted exponent 4.02 on r = 10, 20, 40 and 3.85 if r = 5 is included). Fine as a check, but it carries no information.

Proposed replacement cells (지표 (b)(c) / 통과 기준 (b)(c)):

> (b) 포화 전 구간(S가 반쪽 사슬 큐비트 수의 1/3에 닿기 전)에서 log χ 대 t 직선 맞춤과 log χ 대 log t 직선 맞춤의 잔차 비교 (c) χ = 10⁴에 닿는 (n, t)의 외삽값과, 그 크기 회로의 rG(r = 0.003)

> (b) 지수 맞춤의 잔차 제곱합이 거듭제곱 맞춤의 절반 이하이고, 기울기가 n = 12, 16, 20에서 10% 안에서 일치 (c) n = 12, 16, 20에서 따로 구한 외삽값이 서로, 그리고 이 장의 값(n ≈ 30\~42, t ≈ 3.2\~4.6)과 n·t 모두 2배 안에서 일치하고, 그 크기에서 rG > 5

### B9. ch19 §3 eq.(3) — md:82 vs md:360
Quote: "변분 알고리즘은 훈련에 쓴 회로 실행(평가 수 × 평가당 샷 수)을 R에 넣어야 한다(연습문제 3)."

What is wrong: in eq. (3) R multiplies the number of repetitions ln(0.01)/ln(1 − P). Exercise 3's solution does not put training into R; it adds it once: 2,000 × 1,000 + 228 = 2,000,228. Following §3 literally gives (2,000,000 + 1) × 228 ≈ 4.6×10⁸, a factor 228 off from the chapter's own answer.

Proposed replacement (md:82):

> R은 양쪽이 같은 단위여야 한다. 변분 알고리즘은 훈련에 쓴 회로 실행(평가 수 × 평가당 샷 수)도 같은 단위로 세어 TTS에 한 번 더한다. 훈련은 표본을 뽑을 때마다 다시 하지 않으므로 R에 곱해지지 않는다(연습문제 3).

### B10. ch19 §2, §5(4), verdict table — md:35, md:250, md:290, md:308
Quote (md:250): "그 간격을 메우는 것이 오류 완화(25장)와 오류 정정이다." / (md:290) "가능(오류 완화·정정 필요)" / (md:35) "이득 후보 | 고전 경계 검사를 빠져나가고 자원 추정이 닫히는 문제 | 17장 5단계를 모두 통과"

What is wrong: the boundary circuit has rG = 0.003 × 5,721 ≈ 17. The chapter's own eq. (5) remark (md:100) says "rG가 1을 크게 넘으면 오류 완화의 표본 비용이 지수로 든다", and ch18's classical-boundary table lists "게이트 수 × 오류율 … 1을 크게 넘음" as a caught signal. So by the book's own rules error mitigation is not what closes this gap; a lower error rate or error correction is. The verdict cell and the expected result in md:308 ("오류 완화 없이는 닿지 못함") suggest mitigation would reach it, which no number in the chapter supports.
Related: the role definition (md:35) requires that "자원 추정이 닫히는" and that all five ch17 steps pass. For candidate 4 the only resource number shown does not close (3.4×10⁻⁸) and no fault-tolerant estimate is given; for candidate 8 no step-4 estimate is given at all. Neither "이득 후보" meets the stated definition.

Proposed replacements. md:250 last two sentences:

> 고전이 버거워지는 크기와 지금의 잡음 장치가 닿는 크기 사이에 간격이 있다. 경계 회로는 rG ≈ 17이라 3절 식 (5)의 기준에 걸리고, 오류 완화만으로는 표본 비용이 지수로 든다. 이 간격을 메우는 것은 오류율을 한 자릿수 낮추는 일(r = 3×10⁻⁴이면 rG ≈ 1.7, 무오류 확률 0.18)이나 오류 정정이다(25장).

md:290 판정 cell: "가능(오류율 10배 개선 또는 오류 정정 필요)". md:35 second and third cells:

> 고전 경계 검사를 빠져나가고, 자원 추정에서 남은 간격을 숫자로 적을 수 있는 문제 | 17장 1\~3단계와 5단계 통과, 4단계의 남은 조건을 명시

md:308 예상 cell, last clause: "(c)는 일치하되 결론은 "오류율 0.3%로는 닿지 못함"".

### B11. ch19 §5(3) — md:230
Quote: "이득이 걸린 크기는 표의 아래 두 줄인데, 그쪽은 VQE가 아니라 QPE의 영역이고"

What is wrong: N₂/cc-pVDZ (1.4×10¹² determinants, all-electron, no symmetry) is not beyond exact classical diagonalisation. The largest reported FCI is 1.3×10¹² determinants (C₃H₈/STO-3G; Gao et al., J. Chem. Theory Comput. 20, 1185 (2024), as described in https://arxiv.org/abs/2402.13111: "by far, the largest FCI calculation reported to date … 1.3×10¹² determinants"). With the 1s cores frozen N₂/cc-pVDZ is C(26,5)² = 4.3×10⁹, a standard FCI benchmark for two decades. Only the FeMoco row is clearly out of reach.

Proposed replacement for that sentence:

> 이득을 물을 수 있는 크기는 표의 맨 아래 줄 쪽이다. N₂/cc-pVDZ의 1.4×10¹²은 지금까지 보고된 가장 큰 FCI 계산(1.3×10¹² 행렬식, Gao et al. 2024)과 같은 크기이고, 1s 궤도를 얼리면 4.3×10⁹으로 줄어 FCI 기준값이 오래전부터 있는 분자다. 그보다 큰 쪽은 VQE가 아니라 QPE의 영역이고 초기 상태 겹침이 조건이다(17장).

---

## C — minor

### C1. §1, §5(1) — md:13, md:193 (script:43, 52)
Quote: "99%에 32샷이 든다." The script tests only powers of two. Independent Monte Carlo on the same folded distribution (20,000 trials, `kern_indep.py`): 20 shots 0.987, 22 shots 0.990, 24 shots 0.992, 32 shots 0.996. The threshold is about 22–24 shots. Also "기본 주파수를 99% 맞히려면" is at ±1 bin, i.e. ±15.6 Hz (the fundamental at bin 11.52 is split 0.275/0.278 between bins 11 and 12), whereas ch17's exercise set the accuracy target at "수 Hz".
Replacement (md:193 first sentence): "2의 거듭제곱만 시험했을 때 99%를 넘는 첫 값이 32샷이다(촘촘히 재면 문턱은 약 24샷이고, CNOT은 24 × 1,022 ≈ 2.5만 개다). 여기서 "맞힌다"는 ±1구간, 곧 ±15.6 Hz 안이라는 뜻이다. FFT는 봉우리 보간으로 구간 폭보다 훨씬 잘게 읽는다."

### C2. §5(1) — md:193
Quote: "게이트 하나가 수십 ns인 장치에서도 FFT보다 100배 넘게 느리며". 32,704 sequential CNOTs give 0.65 ms at 20 ns (72×), 0.98 ms at 30 ns (109×). The claim holds only above ~28 ns per gate. With the book's own number (ch14: t_CZ = 68 ns) it is 2.2 ms, about 250×.
Replacement: "14장의 t\_CZ = 68 ns로 차례로 돌리면 2.2 ms로 FFT(9 µs)의 약 250배이며(판독·초기화 시간은 뺀 값이다)".

### C3. §5(1) — md:187 (script:34)
Quote: "기본 주파수 180 Hz와 배음 10개". The script sums k = 1…10, i.e. the fundamental plus nine overtones, with amplitude 1/k.
Replacement: "기본 주파수 180 Hz와 그 2\~10배음(진폭 1/k)".

### C4. §5(1) — md:187, md:191
"기본 주파수 구간(±1)에 실린 확률은 0.564" and the table's 1-shot entry 0.578 are the same quantity; the difference is Monte Carlo noise (2,000 trials, σ ≈ 0.011). Add after the table: "(1샷의 0.578은 2,000회 추정값이고 정확한 값은 위의 0.564다.)"

### C5. §4 — md:181 (script:259–260)
Quote: "스크립트의 외삽으로 약 110일". Timing-dependent: my rerun gave 0.21 s → 91 days (saved stdout 0.25 s → 112 days). The extrapolation also multiplies by 2.638¹⁸, which undercounts: OEIS gives a(27)/a(11) = 881,317,491,628/120,292 = 7.33×10⁶ against 2.638¹⁶ = 5.50×10⁶ (ratio 1.33, the n^{γ−1} correction; https://oeis.org/A001411).
Replacement: "이 코드의 속도로는 잔기 30개쯤에서 몇 달이 걸린다(스크립트의 외삽으로 3\~4달, 실행마다 달라진다)."

### C6. §1, §5(8) — md:18, md:271, md:277 (script:287–291)
"대각 우세 희소(행당 비영 원소 5개)". Checked on the exact matrices the script generates (`mod/ch19_mod.py`): 2, 4 and 28 rows (N = 256, 1,024, 4,096) violate strict diagonal dominance (max off-diagonal row sum 1.10–1.50 against a diagonal of 0.8–1.2), and nonzeros per row average 5.0 but range from 1 to 12–16. The matrices are positive definite (eigenvalues 0.52–1.46), which is what CG needs. For block encoding the maximum row sparsity (12–16), not the mean, sets the cost.
Replacement: "대각이 지배적인 희소 행렬(I + 0.2·S, 행당 비영 원소 평균 5개·최대 16개)" and in md:18 "대각이 지배적인 희소 행렬에서는".

### C7. §5(4) — md:250
(i) "단계를 2배로 늘리면 비충실도가 16배 줄지만": true for 10→20 (16.2×) and 20→40 (16.2×), but 5→10 is 11.0×. Replace with "단계를 2배로 늘리면 비충실도가 약 16배 줄지만(5 → 10단계는 아직 점근 영역이 아니어서 11배)".
(ii) "경계 크기(n = 42, dt = 0.1, 46단계)": dt = 0.1 gave infidelity 4.3×10⁻³ at n = 10, t = 2. Error norm grows roughly like n·t·dt², so at n = 42, t = 4.6 the same dt gives a state infidelity of order (4.2 × 2.3)² × 4.3×10⁻³ ≈ 0.4 (unverified estimate; local observables do better). Add: "dt = 0.1은 n = 10에서 고른 값이라 이 크기의 상태 충실도까지 보장하지는 않는다. CNOT 수는 하한에 가깝다."

### C8. Verdict labels for candidate 3 — md:51, md:123, md:289
§4 puts candidate 3 on the heuristic track (consistent with ch18 md:135), but §3's table and the verdict table cite "17장 5단계", which the heuristic track replaces. Replace the cell in md:289 with "18장 6′(17장 5단계에 해당: 고전 정확 대각화)", and md:51 third cell with "17장 5단계(변분 후보는 18장 6′)".

### C9. §5(8) second table — md:273–277
Header "양자 호출 수(κ에 선형)" with entries "N^2.00 · polylog N". The call count of the optimal algorithm is κ·log(1/ε) with no N-dependence beyond κ; polylog N belongs to the gate cost of one block-encoding call. Replace the header with "양자 게이트 수(호출 수 κ·log(1/ε) × 호출당 polylog N)".

### C10. §3 textbook comparison — md:113
Quote: "입력 적재·출력 읽기·고전 베이스라인을 한 표에서 따지는 일은 Aaronson의 "Read the fine print"(2015)가 선형계와 기계학습에 대해 처음 정리했다." The essay has no table and does not claim priority; it says the original papers already convey the caveats and that they get lost in second-hand accounts (https://scottaaronson.com/papers/qml.pdf). The four caveats themselves match (load b, apply e^{−iAt}/sparsity, condition number, output).
Replacement: "입력 적재·출력 읽기·고전 베이스라인을 함께 따지는 일은 Aaronson의 "Read the fine print"(2015)가 선형계와 기계학습에 대해 널리 알려진 형태로 정리했다."

### C11. §9 fourth note — md:367
Quote: "깨지는 곳: RBF는 폭을 넓히면 그만이지만, 양자 커널은 넓힐수록 양자만의 특징이 함께 사라진다." The same note has just said the cure is to reduce the bandwidth; "넓힐수록" now refers to the RBF-style width (the inverse quantity) without saying so.
Replacement: "깨지는 곳: RBF는 폭을 넓히면 그만이지만, 양자 커널에서 같은 일(회전 각도의 범위를 줄이는 것)을 하면 고전 커널과 다른 부분도 함께 줄어든다."

### C12. Verdict table row 9 — md:295
"18장 6′ | 근사비 0.6924 대 0.878": both numbers are worst-case guarantees (p = 1 on 3-regular graphs vs all graphs), not a same-budget measured comparison, which is what 6′ is. Suggested cell: "18장 6′(보장값 비교)".

### C13. Pre-registration, ch23 row — md:307
"결합 길이 곡선 전체" has no range, and ADAPT has no stopping rule (with an unbounded operator count it reaches FCI trivially at 4 and 12 qubits). Add to 크기와 예산: "결합 길이 H₂ 0.4\~2.5 Å, LiH 1.0\~3.5 Å(0.1 Å 간격), ADAPT는 연산자 수 상한을 UCCSD의 파라미터 수로 둔다". (Whether UCCSD stays within 1.6 mHa for stretched LiH: not checked.)

### C14. Small wording points
- md:263 "도약 항(XX + YY) 10개": 5 bonds, 10 Pauli terms. Suggest "도약 항 10개(결합 5개 × XX, YY)".
- md:16 / md:288: after A1 and B4 are applied, keep the summary and table cells in sync.
- ch18 md:146 says the concentration test "22장에서 한다"; ch19 now does it first. Not an error in ch19; a pointer in ch18 could be updated later.

---

## Summary

The script is sound and every seeded number in the md matches its stdout; independent re-implementations (explicit-gate ZZ feature map, bit-level Heisenberg Hamiltonian, brute-force HP over all 4^(N−1) direction strings, the md code cell run verbatim) reproduce the kernel, quench, Trotter and HP numbers. The problems are in what the numbers are said to mean. (1) The MPS boundary n ≈ 42, t ≈ 4.6 comes from χ = 2^S, which the chapter's own table contradicts by ~6×; with the stated χ definition it is n ≈ 30, t ≈ 3.2. (2) The kernel "collapse to chance" is driven by seven pure-noise PCA features stretched to full bandwidth. (3) The Heisenberg chain is Bethe-integrable while ch17's "가능" was conditioned on non-integrability, and the Schwinger model gets the same verdict without its second deciding number. (4) Four of the five pre-registered criteria are trivially decided or ill-defined as written (ch20 N-independence, ch21 vs random guessing, ch22 NMI ceiling, ch23a (b)/(c)). (5) Two internal contradictions: §7's statement about "탈락" and §3's "R에 넣어야 한다".

## Checked and found correct

- All numbers in md §1, §5, §8, §9 against a fresh run of the script (tables at md:189–191, 197–204, 210–216, 222–228, 234–237, 243–248, 254–259, 267–277, verdict table, exercises E1–E3). Markdown tables all have consistent column counts; no garbled Hangul or stray code points.
- Pitch model: folding of negative frequencies, ±1-bin tolerance, p_fund = 0.564 (bins 11–13: 0.275 + 0.278 + 0.011), 2ⁿ − 2 CNOT for real-amplitude preparation by uniformly controlled Ry (Σ_{k=1}^{n−1} 2ᵏ; Möttönen et al.), 45 controlled phases, N·log₂N = 10,240.
- `zz_feature_states` equals the gate-by-gate circuit H, P(2xᵢ), CX–P(2(π−xᵢ)(π−xⱼ))–CX with full entanglement and reps = 2 (ch06's definition): rows n = 2, 3, 4, 6 recomputed with explicit matrices agree to all printed digits. RBF uses the same scaled features with the median heuristic. Shot formula s ≥ 8K(1−K)/Δ² is right for a 2-standard-error separation of two independent estimates. Chance 1-NN = 24/99.
- FCI dimensions and orbital/electron counts: H₂ (2, 2) → 4; LiH (6, 4) → 225; H₂O (7, 10) → 441; N₂/cc-pVDZ (28, 14) → 1.40×10¹²; FeMoco (54e, 54o) → 3.79×10³⁰, and Reiher et al. do use 54 electrons in 54 spatial orbitals (108 spin orbitals).
- Heisenberg: Majorana degree 4 on every string axis; S(t), χ(t), rate 2.91 bits per unit time, saturation 5.17 bits (65%; n = 12 gives the same 65%); Trotter infidelities; CNOT counts 27r + 15 for n = 10 and 5,721 for n = 42 with merged half-steps; 3 CNOT per bond; (1 − 0.003)^5721 = 3.4×10⁻⁸; infidelity ∝ dt⁴ asymptotically.
- HP: (120292, −5, 16) and all table rows; OEIS A001411 values through 11 steps; valid fraction 2.868%; 4^(N−1) bitstrings and walk counts both use absolute directions including the first step, so they are consistent; "약 1/8" after symmetry is right (15,037 classes); growth 2.719; connective constant 2.638.
- Schwinger: 10 hopping Pauli terms, 6 Z terms, 10 ZZ pairs (all pairs among sites 0–4), degrees 11/11/4.
- Linear systems: closed-form κ for 1D and 2D Poisson, CG iteration counts and residuals (≤ 10⁻⁸), fitted exponents; CG bound nnz·√κ·log(1/ε); Costa et al. O(κ log(1/ε)) and optimal; tridiagonal and multigrid O(N) statements.
- Exercises: 4,094; 131,008; 49,152; 2.7× and 3.2×; N = 2³² ≈ 4×10⁹ → 74.6 h; 4,950 pairs; 2.5×10⁷; 228; 2,228; 65,000; 29× and 31×.
- Consistency with the plan (selection 2, 3, 4, 5, 8 → ch22, 23, 23a, 21, 20; tool assignment in §6; fates of 1, 6, 7, 9, 10) and with reused numbers (ch13 0.6924/0.878, ch14 CZ error 3×10⁻³, ch17 33 days–9,200 years, ch17 χ definition, ch18 Majorana function).
- Citations confirmed on the web: Thanasilp et al. (arXiv:2208.11060, Nat. Commun.); Shaydulin & Wild (PRA 106, 042407; abstract supports bandwidth mitigating kernel decay); Costa et al. (PRX Quantum 3, 040303); Möttönen et al. (QIC 5, 467–473); Robert et al. (10 residues / 22 qubits simulated, 7 residues / 9 qubits on hardware, tetrahedral lattice; DOI s41534-021-00368-4); Reiher et al.; Tindall et al. (PRX Quantum 5, 010308; arXiv 26 Jun 2023); Begušić, Gray, Chan (Sci. Adv. 10, eadk4321; arXiv 9 Aug 2023; sparse Pauli dynamics and tensor networks), so "몇 주에서 몇 달 사이에" is accurate; Aaronson's four caveats; OEIS A001411.
- §9: all five notes and the §2 analogy state where the analogy breaks.

## Could not check

- No Qiskit in the sandbox: `ZZFeatureMap` was compared with ch06's numpy definition and with the documented formula from memory, not executed.
- Citations accepted from memory without opening the source: Huang et al. (Nat. Commun. 12, 2631), Liu–Arunachalam–Temme (Nat. Phys. 17, 1013), Crescenzi et al. and Berger & Leighton (J. Comput. Biol. 5, 423 and 5, 27), Lau & Dill, Martinez et al., Kim et al. (publication date and "beyond brute-force classical computation" wording), Havlíček et al., HHL, the journal volume/page of Thanasilp (15, 5200) and Robert (7, 38), N&C section 4.7 title, Pozsgay et al. PRL 113, 117203.
- OEIS a(28)–a(30) (the fetched page lists terms to a(27) only); the N = 30 time estimate was checked through a(27).
- Whether UCCSD stays within 1.6 mHa for stretched LiH (ch23 expectation), the PEC/ZNE sampling exponent at rG ≈ 17 (only the chapter's own qualitative rule was used), the n = 42 Trotter error estimate in C7(ii), and the "GPU·클러스터에서는 χ 수만" statement in §7 (plausible, no source cited in the chapter).
- The control experiments in B1 (8 and 11 speakers) and B3 (Schwinger entropy) are reviewer scratch runs with their own seeds and parameter choices; they show direction, not book-quality numbers.
