"""10장 QFT와 위상 추정 — 제어 U^{2^j} 위상 되돌림 + 역 QFT 로 고유값 위상을 읽는 회로를 numpy 로 검증한다.
실행: python ch10_qpe.py
환경: Python 3.13, numpy 2.5 (SDK 없이 동작). 큐비트 순서: 빅 엔디안(3장 규약). 계수 레지스터 q[0:t] (q[0] 최상위), 고유상태 레지스터 q[t:].
"""
import numpy as np

np.set_printoptions(precision=4, suppress=True, linewidth=120)
rng = np.random.default_rng(seed=7)

I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
P0 = np.diag([1, 0]).astype(complex); P1 = np.diag([0, 1]).astype(complex)
def Pgate(lam): return np.diag([1, np.exp(1j*lam)])

def kron(*ops):
    out = np.array([[1]], dtype=complex) if ops[0].ndim == 2 else np.array([1], dtype=complex)
    for o in ops: out = np.kron(out, o)
    return out
def single(g, q, n): ops = [I]*n; ops[q] = g; return kron(*ops)
def cphase(c, t_, lam, n):
    a = [I]*n; a[c] = P0; b = [I]*n; b[c] = P1; b[t_] = Pgate(lam); return kron(*a) + kron(*b)
def cnot(c, t_, n):
    a = [I]*n; a[c] = P0; b = [I]*n; b[c] = P1; b[t_] = X; return kron(*a) + kron(*b)
def swap(a_, b_, n): return cnot(a_, b_, n) @ cnot(b_, a_, n) @ cnot(a_, b_, n)
def Hn(n): return kron(*([H]*n))
def seq(gates):
    M = np.eye(gates[0].shape[0], dtype=complex)
    for g in gates: M = g @ M
    return M
def probs(psi): return abs(psi)**2
def bits(i, n): return format(i, f'0{n}b')

# ── 6장의 QFT 재사용 ───────────────────────────────────────────────
def qft_gates(n):
    g = []
    for j in range(n):
        g.append(single(H, j, n))
        for k in range(j+1, n): g.append(cphase(k, j, 2*np.pi/2**(k-j+1), n))
    for j in range(n//2): g.append(swap(j, n-1-j, n))
    return g
def qft(n): return seq(qft_gates(n))
def iqft(n): return qft(n).conj().T

# ── 1. 위상 추정 회로 ────────────────────────────────────────────────
def controlled_U_pow(U, c, t, m, power):
    """계수 큐비트 c 가 제어, 고유상태 레지스터(m 큐비트, 위치 t..t+m-1) 에 U^power. 전체 t+m 큐비트."""
    N = t + m; Up = np.linalg.matrix_power(U, power)
    a = [I]*t; a[c] = P0; A = kron(*a, np.eye(2**m))
    b = [I]*t; b[c] = P1; B = kron(*b, Up)
    return A + B
def qpe_state(U, u, t):
    """U: 2^m 유니터리, u: 고유상태(또는 임의 상태), t: 계수 비트. 반환: 최종 상태 (t+m 큐비트)"""
    m = int(np.log2(len(u))); N = t + m
    psi = kron(np.zeros(2**t, dtype=complex), np.zeros(2**m, dtype=complex)); psi = np.zeros(2**N, dtype=complex)
    psi[:] = 0; psi = np.kron(np.eye(2**t)[0], u)                   # |0...0>|u>
    psi = kron(Hn(t), np.eye(2**m)) @ psi
    for j in range(t):                                               # 계수 큐비트 q[j] 가 제어: U^{2^{t-1-j}} (q[0] 가 최상위 비트)
        psi = controlled_U_pow(U, j, t, m, 2**(t-1-j)) @ psi
    psi = kron(iqft(t), np.eye(2**m)) @ psi
    return psi
def qpe_probs(U, u, t):
    psi = qpe_state(U, u, t); m = int(np.log2(len(u)))
    return probs(psi).reshape(2**t, 2**m).sum(axis=1)                # 계수 레지스터 분포

# 1a. 정확한 경우: U = P(2π·1/8) = T 게이트, φ = 1/8 = 0.001₂
U = Pgate(2*np.pi/8); u = np.array([0, 1], dtype=complex)
for t in (3, 4):
    p = qpe_probs(U, u, t); k = int(np.argmax(p))
    print(f"1a) φ = 1/8, t = {t}: 측정 {bits(k, t)} = {k}/{2**t} = {k/2**t}, 확률 {p[k]:.4f}  (정확한 2진 소수 → 확률 1)")

# 1b. 근사 경우: φ = 1/3 = 0.010101...₂
U = Pgate(2*np.pi/3)
print("\n1b) φ = 1/3 (2진 무한소수)")
for t in (3, 4, 6, 8):
    p = qpe_probs(U, u, t); k = int(np.argmax(p)); best = round(2**t/3)
    # 이론: 최적 정수 k 에서 확률 = |sin(π 2^t δ)/(2^t sin(π δ))|² , δ = φ - k/2^t
    delta = 1/3 - best/2**t; p_th = (np.sin(np.pi*2**t*delta)/(2**t*np.sin(np.pi*delta)))**2
    within1 = p[(best-1) % 2**t] + p[best] + p[(best+1) % 2**t]
    print(f"   t = {t}: 최빈 {bits(k, t)} = {k/2**t:.5f} (오차 {abs(k/2**t - 1/3):.5f}), P(최빈) = {p[k]:.4f} (이론 {p_th:.4f}), P(±1 이내) = {within1:.4f}")
p3 = qpe_probs(U, u, 3); p6 = qpe_probs(U, u, 6)
print("   t = 3 전체 분포:", np.round(p3, 4), "→ 4/π² = 0.405 가 하한 (연습문제 1)")

# 1c. 정밀도 vs t: |추정 - φ| 의 기대값과 성공률 (|오차| ≤ 2^-n 를 성공으로, n = t - 2)
print("\n1c) t 에 따른 오차 (φ = 1/3): t = n + 2 비트로 n 비트 정확도를 확률 ≥ 1 - 1/(2(2^2-2)) = 0.75 로")
for t in (4, 6, 8, 10):
    p = qpe_probs(U, u, t); ks = np.arange(2**t)
    err = np.minimum(abs(ks/2**t - 1/3), 1 - abs(ks/2**t - 1/3))
    n_bits = t - 2
    print(f"   t = {t}: E|오차| = {np.sum(p*err):.5f},  P(|오차| ≤ 2^-{n_bits} = {2**-n_bits:.4f}) = {np.sum(p[err <= 2**-n_bits]):.4f}")

# ── 2. 고유상태 중첩 입력: 고유값들을 확률적으로 샘플링 ───────────────────
# U = diag(e^{2πi·1/8}, e^{2πi·3/8}) 가 아닌, 비대각 U: U = H Z H ... 간단히 2큐비트 U 에 고유상태가 섞인 입력
phis = [1/8, 3/8, 5/8, 3/4]                                         # 4개 고유값 (2큐비트 대각 U), 전부 t=3 에서 정확
D = np.diag([np.exp(2j*np.pi*f) for f in phis])
V = np.linalg.qr(rng.normal(size=(4, 4)) + 1j*rng.normal(size=(4, 4)))[0]  # 난수 고유벡터
U2 = V @ D @ V.conj().T
amps = np.array([0.6, 0.5, 0.4, np.sqrt(1 - 0.36 - 0.25 - 0.16)])     # 고유상태 중첩 계수
u_mix = V @ amps.astype(complex)
p = qpe_probs(U2, u_mix, 3)
print("\n2) 고유상태 중첩 입력 (t=3): 측정 분포 =", np.round(p, 4))
print("   고유값 위상 1/8, 3/8, 5/8, 3/4 에 해당하는 k = 1, 3, 5, 6 의 확률 =", np.round(p[[1, 3, 5, 6]], 4), " vs |계수|² =", np.round(amps**2, 4))
# 측정 후 고유상태 레지스터가 해당 고유상태로 붕괴하는지
psi = qpe_state(U2, u_mix, 3).reshape(8, 4); post = psi[3] / np.linalg.norm(psi[3])
print("   k=3 측정 후 레지스터 상태와 고유벡터 v_1 의 충실도:", abs(np.vdot(V[:, 1], post))**2)

# ── 3. 반복 위상 추정 (Kitaev / IPE): 계수 큐비트 1개, 비트를 최하위부터 ──────
def ipe(U, u, n_bits, shots=1):
    """반복 위상 추정. 각 단계 k: H, 제어 U^{2^{n-1-k}}, 이전 비트로 위상 보정 Rz, H, 측정. 반환 φ 추정(비트열)"""
    m = int(np.log2(len(u))); est_bits = []
    for k in range(n_bits):
        power = 2**(n_bits - 1 - k)
        # 보정 각: 이미 알아낸 하위 비트들 -2π Σ b_j / 2^{j+1} (j 는 상대 위치)
        corr = -2*np.pi*sum(b / 2**(j+2) for j, b in enumerate(est_bits))
        psi = np.kron(np.array([1, 0], dtype=complex), u)
        psi = kron(H, np.eye(2**m)) @ psi
        psi = controlled_U_pow(U, 0, 1, m, power) @ psi
        psi = kron(Pgate(corr), np.eye(2**m)) @ psi
        psi = kron(H, np.eye(2**m)) @ psi
        p1 = probs(psi).reshape(2, 2**m).sum(axis=1)[1]
        # shots 번 재서 다수결
        votes = (rng.random(shots) < p1).sum(); b = int(votes > shots/2) if shots > 1 else int(rng.random() < p1)
        est_bits.insert(0, b)                                        # 최하위 비트부터 얻으므로 앞에 삽입
    return est_bits
U = Pgate(2*np.pi/8); u = np.array([0, 1], dtype=complex)
print("\n3) 반복 위상 추정 (계수 큐비트 1개)")
b = ipe(U, u, 4); print(f"   φ = 1/8, 4비트: 추정 0.{''.join(map(str, b))}₂ = {sum(x/2**(i+1) for i, x in enumerate(b))}  (큐비트 2개로 QPE 와 같은 답)")
U = Pgate(2*np.pi/3)
print("   φ = 1/3: 샷 수에 따른 정확도 (6비트, 200회 반복)")
for shots in (1, 5, 21):
    errs = []
    for _ in range(200):
        b = ipe(U, u, 6, shots); est = sum(x/2**(i+1) for i, x in enumerate(b)); errs.append(min(abs(est - 1/3), 1 - abs(est - 1/3)))
    errs = np.array(errs)
    print(f"   shots/비트 = {shots:>2d}: 중앙 오차 {np.median(errs):.4f}, P(오차 ≤ 1/64) = {(errs <= 1/64).mean():.3f}")

# ── 4. 노이즈: 제어 U 마다 계수 큐비트에 위상 감쇠 → 분포 퍼짐 ──────────────
def phase_damp(rho, lam, q, N):
    K0 = np.array([[1, 0], [0, np.sqrt(1-lam)]]); K1 = np.array([[0, 0], [0, np.sqrt(lam)]])
    A_ = single(K0, q, N); B_ = single(K1, q, N); return A_ @ rho @ A_.conj().T + B_ @ rho @ B_.conj().T
def qpe_noisy(U, u, t, lam):
    m = int(np.log2(len(u))); N = t + m
    psi = np.kron(np.eye(2**t)[0], u); psi = kron(Hn(t), np.eye(2**m)) @ psi
    rho = np.outer(psi, psi.conj())
    for j in range(t):
        C = controlled_U_pow(U, j, t, m, 2**(t-1-j)); rho = C @ rho @ C.conj().T
        for q in range(t): rho = phase_damp(rho, lam, q, N)             # 제어 단계마다 계수 레지스터 전체에
    Q = kron(iqft(t), np.eye(2**m)); rho = Q @ rho @ Q.conj().T
    return np.real(np.diag(rho)).reshape(2**t, 2**m).sum(axis=1)
U = Pgate(2*np.pi/8); u = np.array([0, 1], dtype=complex)
print("\n4) 노이즈 (φ = 1/8, t = 4): 제어 단계마다 계수 큐비트 위상 감쇠 λ → P(정답 0010)")
noisy = {}
for lam in (0, 0.02, 0.05, 0.1, 0.2):
    p = qpe_noisy(U, u, 4, lam); noisy[lam] = p
    print(f"   λ = {lam:.2f}: P(정답) = {p[2]:.4f}, 최빈 = {bits(int(np.argmax(p)), 4)}")

# ── 5. 그림 ───────────────────────────────────────────────────────
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Noto Sans CJK JP"; plt.rcParams["axes.unicode_minus"] = False
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 3.4))
U = Pgate(2*np.pi/3); p_exact = qpe_probs(Pgate(2*np.pi/8), u, 6)
ax1.bar(np.arange(64)/64, p6, width=1/64, color="#4C6EF5", label="φ = 1/3 (무한소수)"); ax1.bar(np.arange(64)/64, p_exact, width=1/64, color="#F03E3E", alpha=0.6, label="φ = 1/8 (정확)")
ax1.axvline(1/3, ls=":", color="#868E96"); ax1.set_xlabel("추정 k/2^6"); ax1.set_ylabel("확률"); ax1.set_title("QPE t=6: 정확 vs 근사 위상"); ax1.legend(frameon=False, fontsize=8)
ts = [4, 6, 8, 10]; errs_t = []
for t in ts:
    p = qpe_probs(U, u, t); ks = np.arange(2**t); err = np.minimum(abs(ks/2**t - 1/3), 1 - abs(ks/2**t - 1/3)); errs_t.append(np.sum(p*err))
ax2.semilogy(ts, errs_t, 'o-', color="#4C6EF5", label="E|오차|"); ax2.semilogy(ts, [2.0**-t for t in ts], '--', color="#868E96", label="2^-t")
ax2.set_xlabel("계수 비트 t"); ax2.set_title("정밀도는 비트당 2배"); ax2.legend(frameon=False, fontsize=8)
for lam, p in noisy.items(): ax3.plot(np.arange(16), p, 'o-', ms=3, label=f"λ={lam}")
ax3.set_xlabel("측정값 k (정답 2)"); ax3.set_ylabel("확률"); ax3.set_title("위상 감쇠가 피크를 퍼뜨린다 (t=4)"); ax3.legend(frameon=False, fontsize=7)
plt.tight_layout(); plt.savefig("ch10_results.png", dpi=150); print("5) saved ch10_results.png")

# ── 게이트 비용 기록 (2큐비트 게이트 / 깊이 / T-count) ──
# 본문 회로를 게이트 리스트로 다시 짜서 세고, 위에서 만든 행렬과 같은지(전역 위상 제외) 확인한다. 난수 미사용.
# 규약(tools/gate_cost.py): 기저 CX + 1큐비트 게이트, 전결합. SWAP = 3 CX, CP(θ) = 2 CX + 위상 회전 3, CP(π) = CZ, CP(0) 생략.
# 깊이 = 게이트를 가능한 한 앞에 놓은(ASAP) 층 수. π/4 의 배수가 아닌 회전은 R 로 세고 T ≈ 3·log2(1/ε)·R (ε = 1e-6).
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tools"))
import gate_cost as gc
_g = gc.G
def qft_list(n, swaps=True):                                     # qft_gates(n) 과 같은 순서
    g = []
    for j in range(n):
        g.append(_g("h", j))
        for k in range(j+1, n): g.append(_g("cp", k, j, angle=2*np.pi/2**(k-j+1)))
    if swaps: g += [_g("swap", j, n-1-j) for j in range(n//2)]
    return g
def qpe_list(lam, t):                                            # U = P(lam) 의 QPE: H^t, 제어 U^{2^{t-1-j}} = CP(lam·2^{t-1-j}), QFT†
    return ([_g("h", j) for j in range(t)] + [_g("cp", j, t, angle=lam*2**(t-1-j)) for j in range(t)]
            + gc.inverse(qft_list(t)))
def qpe_matrix(U, t):                                            # qpe_state 가 쓰는 행렬들의 곱
    return kron(iqft(t), np.eye(2)) @ seq([controlled_U_pow(U, j, t, 1, 2**(t-1-j)) for j in range(t)]) @ kron(Hn(t), np.eye(2))
print("\n6) 게이트 비용 기록 — 기저 CX + 1큐비트, 전결합; SWAP = 3 CX, CP = 2 CX + 위상 회전 3 (CP(π) = CZ); 깊이 = ASAP 층 수; T ≈ 3·log2(1/ε)·R, ε = 1e-6")
for name, lam, tt in (("U = T (φ = 1/8)", 2*np.pi/8, 3), ("U = T (φ = 1/8)", 2*np.pi/8, 4), ("φ = 1/3", 2*np.pi/3, 4), ("φ = 1/3", 2*np.pi/3, 6), ("φ = 1/3", 2*np.pi/3, 8)):
    body = qpe_list(lam, tt); okU = gc.same_up_to_phase(gc.unitary(body, tt+1), qpe_matrix(Pgate(lam), tt)); assert okU
    full = [_g("x", tt)] + body                                   # 고유상태 |1> 준비 포함
    okS = np.allclose(abs(gc.apply(full, np.eye(2**(tt+1))[0], tt+1))**2, probs(qpe_state(Pgate(lam), np.array([0, 1], dtype=complex), tt))); assert okS
    if "1/3" in name: assert gc.two_qubit_count(full) == 2*tt + tt*(tt-1) + 3*(tt//2)
    print(f"   QPE t = {tt}, {name}: {gc.fmt(full)}  [유니터리 = 본문 행렬: {okU}, 분포 = qpe_state: {okS}]")
print("   QPE 일반식(U = 위상 게이트, 각도 일반): 2큐비트 = 2t + t(t−1) + 3⌊t/2⌋  (제어 U^{2^j} t개 + 역 QFT 의 CP t(t−1)/2 개 + SWAP ⌊t/2⌋ 개) →",
      {tt: 2*tt + tt*(tt-1) + 3*(tt//2) for tt in (4, 6, 8)})
print("   QFT n 큐비트: CP n(n−1)/2 개 → CX n(n−1) (+ SWAP ⌊n/2⌋ 개 = 3⌊n/2⌋ CX); 그중 CP(π/2) n−1 개는 T 3개씩(정확), 나머지는 회전 3개씩")
for nn in (3, 4, 6, 8):
    a = qft_list(nn, swaps=True); b_ = qft_list(nn, swaps=False)
    ok_a = gc.same_up_to_phase(gc.unitary(a, nn), qft(nn)); ok_b = gc.same_up_to_phase(gc.unitary(b_, nn), seq(qft_gates(nn)[:len(qft_gates(nn)) - nn//2])); assert ok_a and ok_b
    assert gc.two_qubit_count(b_) == nn*(nn-1) and gc.two_qubit_count(a) == nn*(nn-1) + 3*(nn//2) and gc.t_count(b_) == (3*(nn-1), 3*(nn-1)*(nn-2)//2)
    print(f"     n = {nn}: SWAP 포함 — {gc.fmt(a)};  SWAP 없이 — {gc.fmt(b_)}  [= qft({nn}): {ok_a and ok_b}]")
# 반복 위상 추정 한 라운드: H, 제어 U^{2^m} = CP, 보정 P(corr), H  (계수 큐비트 0, 고유상태 큐비트 1)
def ipe_round_list(lam, power, corr): return [_g("h", 0), _g("cp", 0, 1, angle=lam*power), _g("p", 0, angle=corr), _g("h", 0)]
def ipe_round_matrix(U, power, corr): return kron(H, np.eye(2)) @ kron(Pgate(corr), np.eye(2)) @ controlled_U_pow(U, 0, 1, 1, power) @ kron(H, np.eye(2))
est = []; tot2 = 0; print("   반복 위상 추정 한 라운드 (큐비트 2개), U = T, 4비트 — 라운드 k 는 제어 T^{2^{3−k}} 와 보정 P(corr):")
for k in range(4):
    power = 2**(3-k); corr = -2*np.pi*sum(b / 2**(j+2) for j, b in enumerate(est))
    rd = ipe_round_list(2*np.pi/8, power, corr); okR = gc.same_up_to_phase(gc.unitary(rd, 2), ipe_round_matrix(Pgate(2*np.pi/8), power, corr)); assert okR
    p1 = abs(gc.apply([_g("x", 1)] + rd, np.eye(4)[0], 2).reshape(2, 2)[1, 1])**2; est.insert(0, int(round(p1))); tot2 += gc.two_qubit_count(rd)
    print(f"     k = {k} (T^{power}, corr = {corr/np.pi + 0:.2f}π): {gc.fmt(rd)}  [= 본문 행렬: {okR}, 읽은 비트 {est[0]}]")
print(f"     네 라운드 합: 2큐비트 {tot2}, 읽은 값 0.{''.join(map(str, est))}₂")
corr = -2*np.pi*sum(b / 2**(j+2) for j, b in enumerate([1, 0, 1, 0, 1]))
rd = ipe_round_list(2*np.pi/3, 1, corr); okR = gc.same_up_to_phase(gc.unitary(rd, 2), ipe_round_matrix(Pgate(2*np.pi/3), 1, corr)); assert okR
print(f"   반복 위상 추정 한 라운드, φ = 1/3 (6비트의 마지막 라운드, 보정 있음): {gc.fmt(rd)}  [= 본문 행렬: {okR}]  — 각도가 일반이면 라운드마다 같은 값(첫 라운드만 보정이 없어 R 3)")
