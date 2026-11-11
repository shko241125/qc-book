# R3 findings — citations, links and attributed facts (ch01–ch17)

Reviewer id: R3. Date of web checks: 2026-11-09 (UTC).
Scope: every markdown link, arXiv id, "Author (year)" citation, textbook location, Quirk link and cross-chapter pointer in the 17 `chNN_*.md` files.

Companion files (same directory as this report):
- `/tmp/claude-0/-home-claude/fd63146d-511d-5ae9-8dab-3100fb85f0db/scratchpad/mid/R3_inventory.csv` — 854 extracted items (links, arXiv ids, author–year, textbook locations, cross-refs) with a status column.
- `/tmp/claude-0/-home-claude/fd63146d-511d-5ae9-8dab-3100fb85f0db/scratchpad/mid/R3_quirk_decoded.csv` — the 32 Quirk circuit URLs, decoded JSON.
- `/tmp/claude-0/-home-claude/fd63146d-511d-5ae9-8dab-3100fb85f0db/scratchpad/mid/R3_xref.csv` — 493 chapter-to-chapter pointers (source line, target chapter, context).
- Scripts: `/tmp/claude-0/-home-claude/fd63146d-511d-5ae9-8dab-3100fb85f0db/scratchpad/mid/R3/` (`extract.py`, `quirk.py`, `qsim.py`, `xref.py`, `final_csv.py`).

Method notes. The shell cannot reach arxiv.org (proxy 403), so all web checks went through WebSearch + WebFetch. Evidence URLs are given per finding. Bibliographic data that I confirmed only from my own knowledge (journal/volume/page of well-known papers, N&C/KLM tables of contents) is labelled "memory".

Counts: **A = 4, B = 2, C = 13.** No arXiv id points to a non-existent paper; one id carries the wrong title.

---

## A — wrong paper, or a number/claim the source does not make

### A1. ch13 §10 — `ch13/ch13_variational.md:313` (and §8, line 290)
Current text:
> Gokhale et al., "Minimizing state preparations in VQE by partitioning into qubit-wise commuting cliques", arXiv:1908.11857 (2019) — O(n³) 측정 그룹.

What is wrong. arXiv:1908.11857 is **Gokhale & Chong, "O(N³) Measurement Cost for Variational Quantum Eigensolver on Molecular Hamiltonians"** (two authors, submitted 2019-08-30). The title printed in the book belongs (approximately) to a different paper, arXiv:1907.13623, whose real title is "Minimizing State Preparations in Variational Quantum Eigensolver by Partitioning into **Commuting Families**" (Gokhale, Angiuli, Ding, Gui, Tomesh, Suchara, Martonosi, Chong). The phrase "qubit-wise commuting cliques" is in neither title. The content attributed in line 290 (general commutativity → O(n³) groups) is the 1908.11857 paper, so the id is right and the title/author string is wrong; "et al." is also wrong for a two-author paper.
Evidence: https://arxiv.org/abs/1908.11857 , https://arxiv.org/abs/1907.13623 (both fetched).

Proposed replacement (line 313):
> - Gokhale & Chong, "O(N³) measurement cost for variational quantum eigensolver on molecular Hamiltonians", arXiv:1908.11857 (2019) — 일반 교환 그룹화로 O(n³) 측정 그룹. 선행 논문: Gokhale et al., "Minimizing state preparations in variational quantum eigensolver by partitioning into commuting families", arXiv:1907.13623 (2019).

Proposed replacement (line 290, inside the 해답): replace `(Gokhale et al. 2019)` with
> (Gokhale & Chong 2019)

### A2. ch09 §7 — `ch09/ch09_grover.md:232`
Current text:
> Campbell–Khurana–Montanaro 2019도 SAT에 Grover를 적용할 때 물리 큐비트 수백만 개가 필요하다고 추산했다.

What is wrong. The paper does not say "millions". Its text says the number of physical qubits needed "is very large (e.g. over 10¹²)"; Table 4 (Grover on 14-SAT) lists 1.4×10¹² – 3.1×10¹³ factory qubits, Table 6 (backtracking, graph colouring) 3.6×10¹² – 6.3×10¹³. The book understates by about six orders of magnitude. The paper's actual headline is also different: speedups of 10³–10⁵ over a desktop exist on paper, but "the quantum advantage disappears if one includes the cost of the classical processing power required to perform decoding of the surface code". The paper covers Grover and a quantum backtracking algorithm, for k-SAT and graph colouring.
Evidence: https://quantum-journal.org/papers/q-2019-07-18-167/ (abstract), https://ar5iv.labs.arxiv.org/html/1810.05582 (tables, text).

Proposed replacement:
> Campbell–Khurana–Montanaro 2019는 k-SAT·그래프 채색에 Grover와 양자 백트래킹을 적용하는 비용을 표면 부호까지 넣어 추산했다: 물리 큐비트가 10¹²개 이상 필요하고, 표면 부호 복호에 드는 고전 계산 비용까지 넣으면 양자 이득이 사라진다.

### A3. ch10 §7 — `ch10/ch10_qpe.md:239` (echoed in §9, line 264)
Current text:
> 계수 레지스터를 균등 중첩 대신 "창 함수"(예: 사인 창, Kaiser 창)로 초기화하면 꼬리가 지수적으로 줄어든다(Rendon et al. 2022).

What is wrong. Rendon, Izubuchi, Kikuchi (PRD 106, 034503; arXiv:2110.09590) analyse only the **cosine** tapering window against the rectangular one, and claim a **cubic** (polynomial) improvement: "an upper bound in the cost that implies a cubic improvement with respect to the algorithm's error rate". Per-bin probability goes from ~1/k² to ~1/k⁴ (cumulative tail 1/k → 1/k³). The paper does not mention the Kaiser window and does not claim exponentially decaying tails; it leaves "the exploration of other windows" to future work. Exponential-type suppression is a property of the Kaiser window, which is a different literature (e.g. Greenaway, Pol, Sim, arXiv:2404.01396: "the Kaiser window function is currently the most practical choice for realizing QPE with high success probability").
Evidence: https://arxiv.org/pdf/2110.09590 (fetched), https://arxiv.org/abs/2404.01396v2 (fetched).

Proposed replacement (line 239):
> 계수 레지스터를 균등 중첩 대신 "창 함수"로 초기화하면 꼬리가 훨씬 빨리 줄어든다. 코사인 창은 칸별 확률을 1/k²에서 1/k⁴로(누적 꼬리는 1/k에서 1/k³으로) 줄여 오류율에 대한 비용을 세제곱으로 개선하고(Rendon et al. 2022), Kaiser 창은 그보다 더 강하게 누설을 억제해 현재 가장 실용적인 선택으로 평가된다(Greenaway, Pol, Sim 2024, arXiv:2404.01396).

Proposed replacement (line 264, the parenthesis only): `창 함수로 줄인다(Rendon 2022)` →
> 창 함수로 줄인다(코사인 창: Rendon 2022, 누적 꼬리 1/k → 1/k³)

and add to ch10 §10:
> - Greenaway, S., Pol, W., Sim, S., "A case study against QSVT: assessment of quantum phase estimation improved by signal processing techniques", arXiv:2404.01396 (2024) — Kaiser 창. 등급 A.

### A4. ch11 §7 and §9 — `ch11/ch11_shor.md:221` and `:253`
Current text (221):
> Smolin–Smith–Vargo 2013은 이 논리를 끝까지 밀어 "큐비트 한 개로 어떤 N도 인수분해"하는 풍자 논문을 썼다(컴파일 단계에서 답을 쓰면 회로는 항등이 된다).

Current text (253):
> (Smolin 2013의 풍자: 큐비트 1개로 무엇이든)

What is wrong. The paper says **two** qubits: the method "requires only two coherent qubits" / "can be performed using only two qubits"; the modular exponentiation "reduces to a controlled-NOT gate" (not the identity), and the output distribution equals "tossing an unbiased coin". It also is restricted to N = pq with unequal odd primes, not "any N".
Evidence: https://arxiv.org/abs/1301.7007 (abstract), https://ar5iv.labs.arxiv.org/html/1301.7007 (text). Published as Nature 499, 163 (2013).

Proposed replacement (221):
> Smolin–Smith–Vargo 2013은 이 논리를 끝까지 밀어 "큐비트 두 개로 어떤 N = pq도 인수분해"하는 풍자 논문을 썼다(답을 알고 위수 2인 밑을 고르면 모듈러 지수 회로가 CNOT 하나로 줄고, 출력은 동전 던지기와 같아진다 — RSA-768도 이렇게 "인수분해"된다).

Proposed replacement (253):
> (Smolin 2013의 풍자: 큐비트 2개로 무엇이든)

---

## B — misleading attribution / contradiction between chapters

### B1. ch13 §3 — `ch13/ch13_variational.md:48`
Current text:
> Farhi 2014의 추측을 Mbeng 2019·Bravyi 2020이 상한으로 증명했고, QAOA가 실제로 그 값에 도달함은 Marwaha 2026이 증명했다.

What is wrong. Bravyi–Kliesch–Koenig–Tang did not prove the (2p+1)/(2p+2) upper bound; their paper credits it to earlier work ("Ref. [22] proved that the optimal approximation ratio achieved by level-p QAOA for the ring of disagrees is bounded above by (2p+1)/(2p+2) for all p and conjectured that this bound is tight") and proves a different statement (Theorem 2: any Z₂-symmetric circuit of range R < n/4 gives ⟨H⟩/n ≤ (2R+1/2)/(2R+1)). The upper bound is Mbeng–Fazio–Santoro's (valid for 2P < N). Marwaha's 2026 paper also gives the upper bound in a general form (any depth-p local algorithm that cannot see the whole cycle) together with the achievability proof. Farhi et al. state M_p = n(2p+1)/(2p+2) as a conclusion from numerics for p ≤ 6, which later papers call a conjecture — fine.
Evidence: https://arxiv.org/pdf/1910.08980 , https://ar5iv.labs.arxiv.org/html/1906.08948 , https://arxiv.org/abs/2606.29562 , https://ar5iv.labs.arxiv.org/html/1411.4028 (all fetched).

Proposed replacement:
> Farhi 2014가 p ≤ 6의 수치에서 읽어낸 이 값은 Mbeng 2019가 상한으로 증명했고(2p < n; Bravyi 2020은 Z₂ 대칭인 국소 회로 일반에 대해 비슷한 꼴의 상한을 따로 증명), QAOA가 실제로 그 값에 도달함은 Marwaha 2026이 증명했다.

### B2. ch15 §3 → ch17 — `ch15/ch15_sources.md:59`
Current text:
> 17장 체크리스트의 두 번째 문항 "이 후보는 어느 고전 경계에 걸리는가"가 여기서 나온다.

What is wrong. In ch17 the checklist has five steps (`ch17/ch17_discovery.md:31–35`): 1 문제 정의, 2 구조 찾기, 3 원시 요소 대응, 4 자원 추정, 5 고전 베이스라인. The "which classical boundary" question is step 5, not the second item; step 2 is "구조 찾기". The quoted wording does not appear in ch17 either (grep for "고전 경계" in ch17: no hits).

Proposed replacement:
> 17장 체크리스트의 다섯 번째 단계(고전 베이스라인 — "같은 입력·출력·정확도에서 최선의 고전 방법은 무엇인가")가 여기서 나온다.

---

## C — minor (section numbers, link targets, missing bibliography entries)

### C1. ch13 §10 — `ch13/ch13_variational.md:307`
> arXiv:1411.4028 (2014) — §V 고리(ring of disagrees), 3-정규 0.6924.

In the paper the ring of disagrees is Section IV; Section V is "MaxCut on 3-Regular Graphs" (0.6924). Evidence: https://ar5iv.labs.arxiv.org/html/1411.4028 .
Replacement: `— §IV 고리(ring of disagrees), §V 3-정규 0.6924.`

### C2. ch13 §3/§7 — `ch13/ch13_variational.md:48` and `:278`
> (Bravyi et al. 2020은 일부 그래프에서 상수 p의 QAOA가 고전 국소 알고리즘을 못 이김을 증명)

The paper's headline (abstract, Corollary 1) is that **Goemans–Williamson** outperforms QAOA at any constant level on certain MaxCut instances (bipartite Ramanujan graphs). The "1-local classical algorithm" comparison exists only as Appendix Lemma D.1 for Ising-ring instances. Since both sentences sit next to the GW number, name GW.
Replacement (48): `(Bravyi et al. 2020은 일부 MaxCut 사례에서 어떤 상수 p의 QAOA도 Goemans–Williamson을 못 이김을 증명)`
Replacement (278): `Bravyi 2020은 대칭성 보호 때문에 상수 p 전체가 Goemans–Williamson에 지는 그래프족을 증명했다.`

### C3. ch10 §7 → ch15 — `ch10/ch10_qpe.md:241`
> (a)는 Heisenberg 한계(오차 ∝ 1/회로 길이)를 유지하고 (b)는 샷 노이즈 한계(오차 ∝ 1/√샷)로 후퇴한다 — 이 차이가 15장 "가속의 원천"에서 다시 나온다.

ch15 never discusses Heisenberg vs shot-noise scaling (the only "Heisenberg" in ch15 is "Heisenberg 그림" in an exercise answer). The topic actually returns in ch12 (amplitude estimation 1/M vs 1/√M) and ch16 §3/§5.
Replacement: `— 이 차이가 12장(진폭 추정의 1/M 대 1/√N)과 16장에서 다시 나온다.`

### C4. ch09 §9 — `ch09/ch09_grover.md:260`
> 더 깊이: Brassard–Høyer–Mosca–Tapp 2002, KLM 8.1.

The bullet is about amplitude amplification; the same chapter (line 58 and §10 line 271) says KLM 8.1 = Grover, 8.2 = amplitude amplification. Replacement: `KLM 8.2.`

### C5. ch08 §3 vs §10 — `ch08/ch08_oracles.md:73` vs `:271`
Line 73: "KLM 6장은 네 알고리즘을 모두 상세히 다루고(… 6.4 DJ·BV)". Line 271: "6.3–6.4 Deutsch·DJ(BV는 연습)". The two statements contradict each other. I could not open KLM; from memory the section titles are 6.3 "The Deutsch algorithm", 6.4 "The Deutsch–Jozsa algorithm", 6.5 "Simon's algorithm", with no Bernstein–Vazirani section, so line 271 is the safer one (unverified suspicion).
Replacement (73): `KLM 6장은 Deutsch·Deutsch–Jozsa·Simon을 상세히 다루고(6.2 위상 되돌림, 6.3 Deutsch, 6.4 DJ; BV는 연습문제 수준)`

### C6. ch09 §10 — `ch09/ch09_grover.md:270`
> Nielsen & Chuang, 6.1절(기하적 유도, 해 수 모를 때), 6.3(양자 계수 세기)…

Inconsistent with line 58 of the same chapter ("해 수를 모를 때(6.3, 양자 계수 세기로 해결)"), which matches N&C. Replacement: `6.1절(절차·기하적 유도 6.1.3·성능 6.1.4), 6.3(양자 계수 세기 — 해 수 모를 때), 6.6(최적성), 연습 6.2.`

### C7. ch05 §9 — `ch05/ch05_measurement.md:201` (unverified suspicion)
> 더 깊이: N&C 12.1.1절, Wootters & Zurek (1982).

N&C 12.1.1 is "The Holevo bound"; the no-cloning theorem is Box 12.1 inside section 12.1 (as the chapter's own §10 and exercise 2 say). From memory of the TOC only. Replacement: `N&C 12.1절 Box 12.1`.

### C8. ch04 §10 — `ch04/ch04_universality.md:242`
> [IBM Quantum 프로세서 유형과 기저 게이트](https://docs.quantum.ibm.com/guides/processor-types) — Heron의 {rz, sx, x, cz}.

The page exists (now https://quantum.cloud.ibm.com/docs/en/guides/processor-types) but as fetched today it no longer lists native gates per processor; it says to query the backend ("To find the supported instruction set for a backend, use [this code]") and links a table. The cited fact is not on the cited page any more.
Replacement: `- [IBM Quantum 프로세서 유형](https://quantum.cloud.ibm.com/docs/en/guides/processor-types) — Heron 계열 개요. 기저 게이트 {rz, sx, x, cz}는 이 페이지가 아니라 백엔드의 `backend.operation_names`(또는 문서의 지원 명령 표)에서 확인한다. 등급 A.`

### C9. ch02 §9/§10 — `ch02/ch02_single_qubit_gates.md:227` and `:239`
> Cirq 문서 "Gates and operations — EigenGate exponent" / [Cirq Gates and operations](…/cirq/build/gates) — `EigenGate` 거듭제곱 정의.

The page exists and explains `cirq.X**0.5`, `XPowGate(exponent=t)`, but never names `EigenGate`; there is no section called "EigenGate exponent". Replacement (239): `- [Cirq Gates and operations](https://quantumai.google/cirq/build/gates) — 게이트 거듭제곱(`cirq.X**t` = `XPowGate(exponent=t)`); `EigenGate` 자체의 정의는 API 레퍼런스 `cirq.EigenGate`. 등급 A.` Replacement (227): `Cirq 문서 "Gates and operations"(게이트 거듭제곱)`.

### C10. ch07 §10 — `ch07/ch07_physical_layer.md:271`
> [CUDA-Q — Noisy Simulation](https://nvidia.github.io/cuda-quantum/latest/using/backends/simulators.html)

The URL is the "CUDA-Q Circuit Simulation Backends" overview. Noise is on a separate page it links to, `using/backends/sims/noisy.html` ("Noisy Simulators"). Replacement: `[CUDA-Q — Noisy Simulators](https://nvidia.github.io/cuda-quantum/latest/using/backends/sims/noisy.html)` (path taken from the overview page's own link; I did not fetch the target).

### C11. ch17 §10 — `ch17/ch17_discovery.md:264`
> [qml.estimator](https://docs.pennylane.ai/en/stable/code/qml_estimator.html)

Still returns 200, but the page is now titled "qp.estimator — PennyLane 0.45.1" and the canonical URL is `…/code/qp_estimator.html`; the book's own convention is `qp`. All seven names the book lists (estimate, PauliHamiltonian, Qubitization, SelectPauli, UniformStatePrep, QPE, TrotterPauli) exist. The demo's canonical URL is https://pennylane.ai/demos/re_how_to_use_pennylane_for_resource_estimation (title "How to estimate the resource requirements of quantum algorithms").
Replacement: `[qp.estimator](https://docs.pennylane.ai/en/stable/code/qp_estimator.html)`.

### C12. IBM docs domain — ch01:208, ch02:237–238, ch03:244–245, ch04:241–242, ch05:215, ch06:252, ch08:276, ch09:280, ch10:283 vs ch14:210
Fifteen links use `docs.quantum.ibm.com/...`; every one I tested answers 302 → `quantum.cloud.ibm.com/docs/en/...`. ch14 already uses the new domain. Not broken today, but inconsistent and one redirect away from breaking. Mechanical fix: replace `https://docs.quantum.ibm.com/` with `https://quantum.cloud.ibm.com/docs/en/` in all chapters. Note also that the "latest" API pages are now built from `stable/2.6`, while the chapters say "Qiskit(2.5.x)".

### C13. Works cited in running text but absent from the chapter's §10 list
Bibliographic strings below are from memory, not web-checked.
- ch09:46 "Barenco 1995" → `- Barenco, A. et al., "Elementary gates for quantum computation", Phys. Rev. A 52, 3457 (1995) — 다중제어 게이트 분해. 등급 A.`
- ch11:229 "Hallgren 2002" → `- Hallgren, S., "Polynomial-time quantum algorithms for Pell's equation and the principal ideal problem", STOC 2002. 등급 A.`
- ch11:227 "Ragavan–Vaikuntanathan" → `- Ragavan, S. & Vaikuntanathan, V., "Space-efficient and noise-robust quantum factoring", CRYPTO 2024. 등급 A.`
- ch11:252 "Hardy & Wright 10장" → `- Hardy, G. H. & Wright, E. M., *An Introduction to the Theory of Numbers*, 10장(연분수). 등급 A.`
- ch13:290 "Huggins et al. 2021" → `- Huggins et al., "Efficient and noise resilient measurements for quantum chemistry on near-term quantum computers", *npj Quantum Inf.* 7, 23 (2021).`
- ch13:299 "Ragone et al. 2023" → `- Ragone et al., "A Lie algebraic theory of barren plateaus for deep parameterized quantum circuits", *Nat. Commun.* 15, 7172 (2024) (arXiv 2023).`
- ch13:48, :278 Goemans–Williamson 0.878 → `- Goemans & Williamson, "Improved approximation algorithms for maximum cut and satisfiability problems using semidefinite programming", *J. ACM* 42, 1115 (1995).`
- ch15:229 "Google 2019" → `- Arute et al., "Quantum supremacy using a programmable superconducting processor", *Nature* 574, 505 (2019).`
- ch16:64 "Low–Yoder–Chuang 2016" → `- Low, Yoder, Chuang, "Methodology of resonant equiangular composite quantum gates", *Phys. Rev. X* 6, 041067 (2016).`
- ch12:268 Montanaro has no article number → `*Proc. R. Soc. A* 471, 20150301 (2015)`.
- ch13: "Cerezo 2021" is used for two different papers (line 46/213/274 = Nat. Commun. 12, 1791; line 50 = the Nat. Rev. Phys. review). Suggest writing `Cerezo 2021a`(비용함수 불모 고원) / `Cerezo 2021b`(리뷰) in text and in §10.

---

## Summary

Nineteen findings: 4 A, 2 B, 13 C. The reference lists are in much better shape than is typical — every one of the 17 arXiv ids resolves to a real paper, and the two post-2026-06 ids I could not know (2606.29562, 2607.11637) exist and say what the book says. The real failures are attributed facts: one title/id mismatch (ch13 Gokhale), one number off by 10⁶ (ch09 Campbell et al.), one claim the source does not make (ch10 Rendon "exponential"/Kaiser), one wrong qubit count (ch11 Smolin), and one result credited to a paper that itself credits someone else (ch13 Bravyi). All 32 Quirk circuits decode to valid JSON and, simulated in Quirk's little-endian convention, match the surrounding text. Cross-chapter pointers are sound except ch15→ch17 (wrong step number) and ch10→ch15 (topic not there).

## Checked and found correct

Web-verified (fetched primary page):
- arXiv 2505.15917 (Gidney 2025): <1M physical qubits (897,864), <1 week (4.96 days), 0.1 % gate error, 1 µs cycle, Toffoli ≈ 6.5×10⁹, logical qubits 1,399 (Table 5); three techniques named in ch11:89 are the paper's own. Used consistently in ch11 and ch17.
- Babbush et al. 2021 Table I: 10³ cores → 100 days / 880 years, 10⁶ cores → 280 years / 880 millennia (ch09:232 exact).
- Chakrabarti et al.: arXiv v2 "7.5k logical qubits, T-depth 46 million, 10 MHz" vs Quantum 5, 463 "8k, 54 million, order of a second" (ch17:69, :181, :239, :259 exact).
- The Quantum Insider 2026-03-31: Iceberg Quantum <100k (QLDPC, simulation-based), Google ECC-256 <500k qubits, minutes, circuits withheld, zero-knowledge proof (ch11:223 exact).
- arXiv 2606.29562 (Marwaha, 2026-06-28) proves QAOA reaches (2p+1)/(2p+2) on the ring; arXiv 2607.11637 Table 4 lists ibm_fez / ibm_marrakesh (Heron r2) with CZ 3×10⁻³, readout 1.25–1.45×10⁻², T1 100–300 µs, T2 80–200 µs (ch14:211 exact).
- Mbeng–Fazio–Santoro 1906.08948: bound ε ≥ 1/(2P+2) for 2P < N. Farhi 1411.4028: 0.6924 for 3-regular p = 1. Bravyi et al. abstract: GW beats constant-level QAOA.
- IBM docs: Open Plan 10 min per 28-day rolling window, us-east only; 3-hour job cap; `quantum_seconds` (ch14:210). Pages for Statevector, UGate (matrix matches), quantum_info `partial_trace`/`entropy`, bit-ordering, transpile guide, SolovayKitaev pass, classical-feedforward guide, circuit_library (all 23 names the book uses are listed, none marked deprecated).
- Cirq pages: build/gates, build/circuits (Moments), build/custom_gates, transform/transformers (`optimize_for_target_gateset`), transform/routing_transformer (`RouteCQC`), cirq/noise; Shor tutorial confirms `examples/shor.py` with `ArithmeticGate`.
- PennyLane 0.45.1 pages: templates, GroverOperator, QuantumPhaseEstimation (parameter names match), estimator.estimate, resource-estimation demo. Qiskit Aer 0.18.2 `thermal_relaxation_error`. qutip-qip "Pulse-level circuit simulation". QuTiP `#module-qutip.measurement` anchor (on v5.0.0 docs).

Confirmed from memory only (not fetched): journal/volume/page/year of all other papers in the §10 lists of ch02–ch17 (about 95 entries — CHSH, Aaronson–Gottesman, Dawson–Nielsen, Ross–Selinger, BBCJPW, Wootters–Zurek, Draper, Havlíček, Kandala, McClean, Schuld, Krantz, Lindblad, GKS, Li et al., Rauch, Werner, Wallman–Emerson, Deutsch–Jozsa, BV, Simon, Beals, Grover, BBBV, Zalka, BBHT, BHMT incl. Theorem 12, Babbush 2018/2021, Regev–Schiff, Kitaev, Cleve et al., Dobšíček, Wiebe–Granade, Lee et al., Shor, Beauregard (2n+3 qubits, O(n³ log n) gates), Vandersypen, Regev, Lloyd, Suzuki, Childs et al. 2018/2021, Low–Chuang 2017/2019, Montanaro, Kim et al., Peruzzo, Wecker, Mitarai, Bravyi 2016/2019/2020, Cerezo ×3, Bittel–Kliesch, Larocca, Wierichs, Gonthier, Magesan, Nielsen 2002, SABRE, Nation, ADH, Stim (abstract numbers), Vidal, Jozsa–Linden, Page, Fattal, Ambainis, Aaronson–Ambainis, Raz–Tal, Tang, Kerenidis–Prakash, Pan–Chen–Zhang (512 GPUs, 15 h), Konno, Valiant, Terhal–DiVincenzo, Childs–Wiebe, Berry et al., GSLW, Martyn, Dong et al., Suzuki 2020, Grinko, Szegedy, MNRS, Aaronson 2015, Lin Lin, Hoefler et al., Calabrese–Cardy, Jordan et al.). No inconsistency of year/authors/id for the same paper across chapters.

Textbook locations, all from memory of the tables of contents: N&C 1.2, 1.2.1, 1.3.2, 1.3.6, 1.3.7, 1.4.3–1.4.4, 2.1.7, 2.2.3–2.2.6, 2.4.3, 2.6, 4.2 (Thm 4.1, Cor 4.2), 4.3 (Fig 4.9), 4.4, 4.5.2–4.5.5, 4.7 (Thm 4.3, Fig 4.19, Ex 4.49–4.51), 5.1 (Fig 5.1), 5.2/5.2.1 (Figs 5.2–5.3, bound 0.75 for p = 2), 5.3 (Thm 5.1, Thm 5.3 = A4.13), 5.4, 6.1.3, 6.1.4, 6.3, 6.6, 7.7.2, 8.2–8.4, 8.3.3–8.3.6, 10.5.2, 10.5.4, App. 3, App. 4, Box 12.1; exercises 1.2, 4.13–4.14, 4.20, 5.6, 5.8–5.9, 5.13, 6.2; KLM ch 2–3, 3.3–3.5, 4.2–4.4, 5.1–5.2, 6.2–6.5, 7.1–7.5, 8.1–8.4, 9, 10; Hardy & Wright ch. X; Aaronson *Democritus* ch. 10 — confirmed or plausible. Suspicious ones are C4–C7.

Quirk: 32 circuit links, all valid JSON; simulated with a small numpy model of Quirk (top wire = least significant bit; `QFT3` on the little-endian integer). Every stated outcome matches: ch08 kickback/Deutsch/BV(101), ch09 Grover amplitudes 0.354 → 0.884/0.177, P = 0.781 and 0.945, ch10 QPE 001 and Bloch longitudes 45°/90°/180°, ch10 φ = 1/16 split 0.411/0.411, ch06 `QFT3` ≡ H/controlled-phase/SWAP expansion (overlap 1.0) with the stated per-wire phases, ch13 P(ON) = 0 and 0.146, ch16 LCU ancilla-0 probabilities 1/4 + 1/4, ch17 peaks at 0 and 4, ch05 teleportation, ch14 SWAP chain.

Cross-chapter pointers with a section or exercise number (17) all land on the right content: 2장 7절, 4장 연습문제 3, 6장 2절, 6장 연습문제 2 (×2), 7장 연습문제 1 (×2), 8장 7절, 10장 3절, 10장 연습문제 2, 10장 7절, 12장 6절, 12장 7절 (×2), 14장 4절, 15장 5절, 16장 5절. Spot-checked numeric pointers also agree (ch05 → ch14 readout 1.0–1.6×10⁻², ch07 → ch14 32 ns, ch09/ch14 C³Z = 1 ancilla + 3 Toffoli, ch13 → ch12 E₀ = −4.7588, ch11 → ch15 table first row).

## Could not check

- KLM and N&C page-level content (no access to the books); all textbook locations rest on memory.
- `qutip.readthedocs.io/en/stable/...` (two links) and `qutip.readthedocs.io/` root, `github.com/quantumlib/Cirq/blob/main/examples/shor.py`, `algassert.com/quirk` home: not fetchable in this session; checked indirectly (other doc versions / the Cirq tutorial).
- docs.quantum.ibm.com `compiler#transpile` anchor and the old-domain `guides/*` URLs individually (only the redirect pattern and the new-domain pages were confirmed).
- Whether the SDK snippets call any deprecated argument: only noted that PennyLane 0.45.1 deprecates the `id` keyword of `GroverOperator` / `QuantumPhaseEstimation` (removal in v0.46); I did not audit snippet bodies.
- The remaining ~430 un-numbered "N장" pointers were grouped and skimmed, not individually verified.

---

## Promises to future chapters (18–27)

Collected from every "N장" pointer with N ≥ 18 (source line in brackets).

**18장**
- dequantization 검사 [ch15:229]; 입력 모형을 고전에도 주면 가속이 사라지는 경우 [ch16:315].
- 체크리스트 각 단계의 자동화 도구와 반례(dequantization) [ch17:13].
- 탐색 도구(회로 합성·변분 탐색) — 휴리스틱 후보를 "실제로 돌려 보는" 경우 [ch17:227].

**19장**
- 이 절차로 **후보 10개**를 다시 평가 [ch17:13] — note ch17 itself evaluates eight.
- 체크리스트 왼쪽 기둥을 통과한 후보만 19장으로 [ch17:39]; 자원 기록 표를 후보마다 채움 [ch17:94].

**20장 (HHL·선형계)**
- 진폭 증폭을 부품으로 재사용, "후선택 성공률 올리기" [ch09:13, :54, :238, :260].
- QPE가 HHL의 엔진 [ch10:13]; 블록 인코딩으로 행렬 넣기, QSVT 역행렬 [ch16:34, :36]; 16장 부품의 "첫 응용" [ch16:9].
- 출력이 기댓값일 때만 지수 이득이 남는다는 조건 [ch17:183, :229]; 분류표 "선형대수" 행 [ch15:41, ch17:48].

**21장 (QUBO·QAOA)**
- `n_local` ansatz 재사용 [ch06:9]; 13장 진단 도구 위의 QUBO 사례 [ch13:13]; "작은 p의 QAOA는 고전에 진다"를 전제로 출발 [ch13:278].
- HP 격자 QAOA [ch17:227].

**22장 (양자 커널)**
- feature map을 양자 커널로 사용 [ch06:9]; Havlíček의 "고전 시뮬레이션이 어렵다" 추정 검증 [ch06:45].
- 커널 집중 현상 [ch06:147]; 같은 데이터에 RBF 커널 베이스라인 [ch06:213]; 세 하이퍼파라미터를 실제 임베딩 데이터로 탐색 [ch06:229].
- 탈락 후보(입력 적재가 이득을 지우는 주기성 후보)의 반례 박스 [ch17:235].

**23장 (분자 바닥상태·VQE)**
- UCCSD ansatz [ch06:47, :211, ch13:34]; `n_local`/VQE 사례 [ch06:9, ch13:13].
- 진폭 증폭으로 바닥상태 준비 [ch09:13, :54, :260].
- Hartree–Fock 겹침 γ를 "분자 하나로 직접 재 본다" [ch10:237]; γ가 지수로 작아질 때의 검증 [ch10:255]; 상태 준비 문제 재론 [ch12:231]; 화학에는 공짜 입력 상태가 없음 [ch11:250].
- α/ε가 커지는 양자 화학 자원 추정 [ch17:243]; 해밀토니안 시뮬레이션 부품 위의 화학 [ch12:13].

**23a장 (스핀 동역학)** — a chapter numbered "23a" is promised in three places and is the surviving candidate of ch17: [ch12:13], [ch17:147], [ch17:180]. Check that the plan really has a 23a.

**24장**
- 백엔드별 트랜스파일 설정 [ch04:207].
- 다른 클라우드(IonQ, Quantinuum, IQM)의 사양서를 14장과 같은 표로 비교 [ch14:182].

**25장 (노이즈·오류 완화·오류 정정)** — the most heavily loaded target; it is variously called "노이즈", "오류 완화", "오류 정정".
- QuTiP의 밀도행렬 관점 사용 [ch01:139]; 측정 오류율 [ch01:175]; 양자 궤적 관점 [ch05:147].
- POVM·일반 측정·약한 측정 [ch05:47, :180]; Kraus로 쓴 노이즈 채널 [ch05:194, :202].
- 측정 오류 완화(혼동 행렬), "25장의 첫 번째 오류 완화" [ch05:178, ch14:198].
- 오류 정정: 중간 측정이 전제 [ch05:11], 복제 불가와의 관계 [ch05:190], 신드롬 복호 [ch05:200], T가 비싼 이유 [ch15:55, :252], Stim 기반 시뮬레이션 [ch15:219], 논리 대 물리 자원 [ch17:251].
- 동적 디커플링 (`PadDynamicalDecoupling`) [ch07:229, ch14:178]; 무작위 컴파일링 [ch07:231, :255]; 영잡음 외삽(ZNE) [ch14:129]; "완화는 편향을 줄이고 분산을 키운다" [ch14:180].
- 실제 모듈러 곱셈 회로 + 7장 노이즈로 N = 15 Shor 성공률 실험 [ch11:231].
- 14장의 예측 표 위에 서는 장 [ch07:13, ch14:13].

**26장 (GPU·텐서망 시뮬레이션)**
- 상태벡터 한계 다음 단계로서 GPU·텐서망 [ch03:210]; (2^t, 2^n) 배열 기법 [ch11:55].
- CUDA-Q로 n ≥ 30 Grover [ch09:224]; n = 30 스핀 사슬 (`exp_pauli` 커널) [ch12:223]; n = 25 불모 고원 실험 [ch13:266]; 타깃만 바꿔 다른 경계 쓰기 [ch15:213].

**27장 (재현 가능 워크플로)**
- 14장 표 위에 서는 장 [ch14:13]; 백엔드 이름·날짜·사용 큐비트 기록 규칙 [ch14:176]; "예측 → 가짜 백엔드 → 실기" 순서 고정 [ch14:206].

Un-numbered promises worth tracking: "후속 하드웨어 장"에서 텔레포테이션 실기 비교 [ch05:141] 및 "컴파일된" N = 15 Shor 실기 실행 [ch11:221]; 14장 노이즈 모델로 예측해 볼 과제로 남긴 것 — Grover k = 1 대 k = 2 [ch09:236], 텔레포테이션 표 [ch05:141], QAOA/HVA의 p 의존성 [ch13:278].
