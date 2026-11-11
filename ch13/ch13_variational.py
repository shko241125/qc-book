"""13장 변분 알고리즘 — numpy 검증 스크립트
VQE(TFIM n=4, 하드웨어 효율 ansatz), 파라미터 이동 규칙, 샷 잡음, QAOA MaxCut, 불모 고원(barren plateau)
규약: big-endian (q0 = 최상위 비트), seed=7
"""
import numpy as np
from scipy.optimize import minimize
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
matplotlib.rcParams["font.family"] = ["Noto Sans CJK JP", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

rng = np.random.default_rng(7)
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]])
Z = np.diag([1.0, -1.0]).astype(complex)
H1 = np.array([[1, 1], [1, -1]]) / np.sqrt(2)


# ---------- 상태벡터 유틸 (big-endian: 큐비트 q 는 텐서 축 q) ----------
def apply_1q(psi, g, q, n):
    psi = psi.reshape([2] * n)
    psi = np.moveaxis(np.tensordot(g, psi, axes=([1], [q])), 0, q)
    return psi.reshape(-1)


def apply_cz(psi, a, b, n):
    psi = psi.reshape([2] * n).copy()
    idx = [slice(None)] * n
    idx[a] = 1
    idx[b] = 1
    psi[tuple(idx)] *= -1
    return psi.reshape(-1)


def ry(t):
    c, s = np.cos(t / 2), np.sin(t / 2)
    return np.array([[c, -s], [s, c]], complex)


def rx(t):
    c, s = np.cos(t / 2), np.sin(t / 2)
    return np.array([[c, -1j * s], [-1j * s, c]])


def pauli_matrix(s):
    m = {"I": I2, "X": X, "Y": Y, "Z": Z}
    out = np.array([[1.0 + 0j]])
    for ch in s:
        out = np.kron(out, m[ch])
    return out


def tfim_terms(n, J=1.0, h=1.0):
    terms = []
    for i in range(n - 1):
        s = ["I"] * n; s[i] = "Z"; s[i + 1] = "Z"
        terms.append((-J, "".join(s)))
    for i in range(n):
        s = ["I"] * n; s[i] = "X"
        terms.append((-h, "".join(s)))
    return terms


def ham_matrix(terms, n):
    Hm = np.zeros((2 ** n, 2 ** n), complex)
    for c, s in terms:
        Hm += c * pauli_matrix(s)
    return Hm


# ---------- 하드웨어 효율 ansatz: L 층 × [Ry(θ) 전 큐비트 → CZ 사슬] + 마지막 Ry 층 ----------
def n_params(n, L):
    return n * (L + 1)


def ansatz_state(theta, n, L):
    psi = np.zeros(2 ** n, complex); psi[0] = 1.0
    k = 0
    for layer in range(L + 1):
        for q in range(n):
            psi = apply_1q(psi, ry(theta[k]), q, n); k += 1
        if layer < L:
            for q in range(n - 1):
                psi = apply_cz(psi, q, q + 1, n)
    return psi


def apply_cx(psi, c, t, n):
    psi = apply_1q(psi, H1, t, n)
    psi = apply_cz(psi, c, t, n)
    return apply_1q(psi, H1, t, n)


def ansatz_state_cx(theta, n, L):
    """Qiskit RealAmplitudes(entanglement='linear') 와 같은 구조: Ry 층 + CX 사슬."""
    psi = np.zeros(2 ** n, complex); psi[0] = 1.0
    k = 0
    for layer in range(L + 1):
        for q in range(n):
            psi = apply_1q(psi, ry(theta[k]), q, n); k += 1
        if layer < L:
            for q in range(n - 1):
                psi = apply_cx(psi, q, q + 1, n)
    return psi


def hva_state(theta, n, p, zz_diag):
    """해밀토니안 변분 ansatz (HVA): |+>^n 에서 e^{-iγ H_ZZ} e^{-iβ H_X} 를 p 층 (TFIM 전용, 2p 파라미터)."""
    psi = np.full(2 ** n, 1 / np.sqrt(2 ** n), complex)
    for l in range(p):
        psi = np.exp(-1j * theta[2 * l] * zz_diag) * psi
        for q in range(n):
            psi = apply_1q(psi, rx(-2 * theta[2 * l + 1]), q, n)   # e^{-iβ(-X)} = Rx(-2β)
    return psi


ANSATZ = {"Ry+CZ": ansatz_state, "Ry+CX": ansatz_state_cx}


def energy(theta, Hm, n, L, kind="Ry+CZ"):
    psi = ANSATZ[kind](theta, n, L)
    return float(np.real(psi.conj() @ Hm @ psi))


def grad_param_shift(theta, Hm, n, L, kind="Ry+CZ"):
    """Ry(θ) = e^{-iθY/2}: 생성원 고유값 ±1/2 → shift π/2, 계수 1/2 (Schuld et al. 2019)."""
    g = np.zeros_like(theta)
    for k in range(len(theta)):
        tp = theta.copy(); tp[k] += np.pi / 2
        tm = theta.copy(); tm[k] -= np.pi / 2
        g[k] = 0.5 * (energy(tp, Hm, n, L, kind) - energy(tm, Hm, n, L, kind))
    return g


def grad_finite_diff(theta, Hm, n, L, kind="Ry+CZ", eps=1e-5):
    g = np.zeros_like(theta)
    for k in range(len(theta)):
        tp = theta.copy(); tp[k] += eps
        tm = theta.copy(); tm[k] -= eps
        g[k] = (energy(tp, Hm, n, L, kind) - energy(tm, Hm, n, L, kind)) / (2 * eps)
    return g


# ======================================================================
n = 4
terms = tfim_terms(n)
Hm = ham_matrix(terms, n)
evals, evecs = np.linalg.eigh(Hm)
E0, E1 = evals[0], evals[1]
gs = evecs[:, 0]
print(f"1) TFIM n={n}, J=h=1: E0 = {E0:.6f}, E1 = {E1:.6f}  (12장과 동일)")

# ---------- 2) 파라미터 이동 규칙 vs 유한차분 ----------
L = 2
theta0 = rng.uniform(0, 2 * np.pi, n_params(n, L))
g_ps = grad_param_shift(theta0, Hm, n, L, "Ry+CX")
g_fd = grad_finite_diff(theta0, Hm, n, L, "Ry+CX")
print(f"2) 파라미터 이동 vs 유한차분(eps=1e-5): max|차이| = {np.max(np.abs(g_ps - g_fd)):.2e}, "
      f"|grad| = {np.linalg.norm(g_ps):.4f}, 파라미터 {len(theta0)}개 → 회로 평가 {2*len(theta0)}회/그래디언트")

# ---------- 3) VQE: ansatz 종류 × 층 수 ----------
print("3) VQE (L-BFGS-B, 파라미터 이동 그래디언트, 무작위 초기값 5개 중 최선)")
zz_diag = np.zeros(2 ** n)
_bits = ((np.arange(2 ** n)[:, None] >> (n - 1 - np.arange(n))) & 1); _zv = 1 - 2 * _bits
for i in range(n - 1):
    zz_diag += -1.0 * _zv[:, i] * _zv[:, i + 1]


def run_vqe(fun, jac, npar, starts=5):
    best = None
    for trial in range(starts):
        th = rng.uniform(0, 2 * np.pi, npar)
        hist = []
        res = minimize(fun, th, jac=jac, method="L-BFGS-B",
                       callback=lambda t: hist.append(fun(t)), options={"maxiter": 300, "gtol": 1e-9})
        if best is None or res.fun < best[0].fun:
            best = (res, hist)
    return best


vqe_hist = {}
vqe_summary = {}
for kind in ("Ry+CZ", "Ry+CX"):
    for L in (1, 2, 3, 4):
        res, hist = run_vqe(lambda t: energy(t, Hm, n, L, kind), lambda t: grad_param_shift(t, Hm, n, L, kind),
                            n_params(n, L))
        psi = ANSATZ[kind](res.x, n, L)
        fid = abs(gs.conj() @ psi) ** 2
        vqe_hist[(kind, L)] = hist
        vqe_summary[(kind, L)] = (res.fun, fid, res.nfev, res.x)
        print(f"   {kind} L={L}: 파라미터 {n_params(n, L):2d}, 2큐비트 {(n-1)*L:2d}, 깊이≈{2*L+1:2d} | "
              f"E-E0 = {res.fun-E0:.2e}, 충실도 {fid:.6f}, 에너지 평가 {res.nfev} 회")
for p in (1, 2, 3, 4):
    res, hist = run_vqe(lambda t: float(np.real((lambda s: s.conj() @ Hm @ s)(hva_state(t, n, p, zz_diag)))),
                        None, 2 * p, starts=10)
    psi = hva_state(res.x, n, p, zz_diag)
    fid = abs(gs.conj() @ psi) ** 2
    vqe_hist[("HVA", p)] = hist
    vqe_summary[("HVA", p)] = (res.fun, fid, res.nfev, res.x)
    print(f"   HVA p={p}: 파라미터 {2*p:2d}, 2큐비트 {2*(n-1)*p:2d} (e^{{-iγZZ}} 당 CNOT 2) | "
          f"E-E0 = {res.fun-E0:.2e}, 충실도 {fid:.6f}, 에너지 평가 {res.nfev} 회")


# 동역학 리 대수(DLA) 차원: 회로가 닿을 수 있는 방향의 수 (Larocca et al. 2022)
def layer_unitary(apply_layer):
    U = np.eye(2 ** n, dtype=complex)
    return np.array([apply_layer(U[:, c]) for c in range(2 ** n)]).T


def dla_dim(gens, return_basis=False):
    basis = []
    def add(M):
        v = M.reshape(-1)
        if basis:
            B = np.array(basis).T
            coef = np.linalg.lstsq(B, v, rcond=None)[0]
            if np.linalg.norm(v - B @ coef) < 1e-8:
                return False
        basis.append(v / np.linalg.norm(v)); return True
    queue = [g for g in gens if add(g)]
    i = 0
    while i < len(queue):
        A = queue[i]; i += 1
        for Bm in list(queue):
            C = A @ Bm - Bm @ A
            if np.linalg.norm(C) > 1e-10 and add(C):
                queue.append(C)
    return basis if return_basis else len(basis)


def cz_chain(psi):
    for q in range(n - 1): psi = apply_cz(psi, q, q + 1, n)
    return psi


def cx_chain(psi):
    for q in range(n - 1): psi = apply_cx(psi, q, q + 1, n)
    return psi


Ys = [pauli_matrix("".join("Y" if k == q else "I" for k in range(n))).imag for q in range(n)]   # Y = i·A, A 실수 반대칭
dla_basis = {}
for name, chain in (("Ry+CZ", cz_chain), ("Ry+CX", cx_chain)):
    Uc = layer_unitary(chain).real
    order, M = 1, Uc.copy()
    while not np.allclose(M, np.eye(2 ** n)): M = M @ Uc; order += 1          # 얽힘층의 차수: CZ 사슬 2, CX 사슬 4
    gens, Wk = [], np.eye(2 ** n)
    for k in range(order):                                                     # 층별 회전은 W^{-k} Y W^k 로 나타난다
        gens += [Wk.T @ A @ Wk for A in Ys]; Wk = Wk @ Uc
    dla_basis[name] = dla_dim(gens, return_basis=True)
    print(f"   DLA 차원 {name}: {len(dla_basis[name])}  (얽힘층 차수 {order}; 실수 상태 전체 so(16) = 120)")
# 도달 가능성은 대수 차원이 아니라 궤도 G|0...0> 로 판정: 궤도의 접공간 차원 = rank{A|ψ> : A ∈ 대수}
A_cz = [b.reshape(2 ** n, 2 ** n) for b in dla_basis["Ry+CZ"]]
for lab, psi in (("|0000>", np.eye(2 ** n)[0]), ("TFIM 바닥상태", gs.real)):
    T = np.array([a @ psi for a in A_cz])
    print(f"   Ry+CZ 궤도 접공간 차원 at {lab}: {np.linalg.matrix_rank(T, tol=1e-8)}  (실수 단위구 = 15)")
# 무작위 실수 목표 상태를 Ry+CZ L=6 으로 직접 맞추기 (충실도 최대화)
infid = []
for trial in range(5):
    tgt = rng.normal(size=2 ** n); tgt /= np.linalg.norm(tgt)
    best = min(minimize(lambda t: 1 - abs(tgt @ ansatz_state(t, n, 6)) ** 2, rng.uniform(0, 2 * np.pi, n_params(n, 6)),
                        method="BFGS").fun for _ in range(4))
    infid.append(best)
print(f"   Ry+CZ L=6, 무작위 실수 목표 5개의 최소 비충실도: {np.round(infid, 3)}  (so(16) 전체면 0)")

# ---------- 4) 샷 잡음: 두 측정 기저(ZZ 그룹, X 그룹)로 에너지 추정 ----------
def energy_shots(theta, n, L, shots, rng):
    psi = ansatz_state_cx(theta, n, L)
    # Z 기저: ZZ 항 (모두 교환) 한 번에
    p = np.abs(psi) ** 2
    samp = rng.choice(2 ** n, size=shots, p=p)
    bits = ((samp[:, None] >> (n - 1 - np.arange(n))) & 1)      # big-endian
    z = 1 - 2 * bits
    e_zz = sum(-1.0 * np.mean(z[:, i] * z[:, i + 1]) for i in range(n - 1))
    # X 기저: H^{⊗n} 후 Z 측정
    psix = psi
    for q in range(n):
        psix = apply_1q(psix, H1, q, n)
    px = np.abs(psix) ** 2
    sampx = rng.choice(2 ** n, size=shots, p=px)
    bitsx = ((sampx[:, None] >> (n - 1 - np.arange(n))) & 1)
    xz = 1 - 2 * bitsx
    e_x = sum(-1.0 * np.mean(xz[:, i]) for i in range(n))
    return e_zz + e_x


L = 3
theta_opt = vqe_summary[("Ry+CX", 3)][3]
E_exact_opt = energy(theta_opt, Hm, n, L, "Ry+CX")
print(f"4) 샷 잡음 (Ry+CX L=3 최적점, E_exact = {E_exact_opt:.4f}; 기저 2개 × shots)")
shot_list = [100, 400, 1600, 6400, 25600]
shot_std = []
for s in shot_list:
    ests = np.array([energy_shots(theta_opt, n, L, s, rng) for _ in range(200)])
    shot_std.append(ests.std())
    print(f"   shots={s:6d}: 평균 {ests.mean():.4f}, 표준편차 {ests.std():.4f}, std·√shots = {ests.std()*np.sqrt(s):.3f}")
ratio = shot_std[0] / shot_std[-1]
print(f"   shots 256배 → std {ratio:.1f}배 감소 (1/√N 이면 16)")

# 샷 잡음 하에서 그래디언트 하강 (파라미터 이동, 1000 shots/평가)
L = 2
th = rng.uniform(0, 2 * np.pi, n_params(n, L))
lr = 0.1
noisy_hist = []
for it in range(150):
    g = np.zeros_like(th)
    for k in range(len(th)):
        tp = th.copy(); tp[k] += np.pi / 2
        tm = th.copy(); tm[k] -= np.pi / 2
        g[k] = 0.5 * (energy_shots(tp, n, L, 1000, rng) - energy_shots(tm, n, L, 1000, rng))
    th = th - lr * g
    noisy_hist.append(energy(th, Hm, n, L, "Ry+CX"))
print(f"   1000 shots/평가 그래디언트 하강 150 스텝 (Ry+CX L=2): 최종 정확 에너지 {noisy_hist[-1]:.4f} "
      f"(L=2 무잡음 최적 E-E0 {vqe_summary[("Ry+CX", 2)][0]-E0:.1e}, E0 {E0:.4f}); 총 회로 실행 {150*2*len(th)*2*1000:,} 샷")

# ---------- 5) QAOA MaxCut: 6-노드 고리 ----------
nq = 6
edges = [(i, (i + 1) % nq) for i in range(nq)]
maxcut = nq  # 짝수 고리: 교대로 자르면 모든 변
# 비용 해밀토니안 C = Σ_(i,j) (1 - Z_i Z_j)/2 — 대각
z_bits = ((np.arange(2 ** nq)[:, None] >> (nq - 1 - np.arange(nq))) & 1)
zval = 1 - 2 * z_bits
cost_diag = np.zeros(2 ** nq)
for (i, j) in edges:
    cost_diag += 0.5 * (1 - zval[:, i] * zval[:, j])


def qaoa_state(gammas, betas):
    psi = np.full(2 ** nq, 1 / np.sqrt(2 ** nq), complex)
    for g, b in zip(gammas, betas):
        psi = np.exp(-1j * g * cost_diag) * psi
        for q in range(nq):
            psi = apply_1q(psi, rx(2 * b), q, nq)
    return psi


def qaoa_expect(params):
    p = len(params) // 2
    psi = qaoa_state(params[:p], params[p:])
    return float(np.real(np.sum(np.abs(psi) ** 2 * cost_diag)))


print(f"5) QAOA MaxCut, {nq}-노드 고리 (최대 컷 {maxcut})")
# p=1 지형
G = np.linspace(0, np.pi, 61)
B = np.linspace(0, np.pi / 2, 31)
land = np.array([[qaoa_expect(np.array([g, b])) for g in G] for b in B])
ib, ig = np.unravel_index(np.argmax(land), land.shape)
print(f"   p=1 격자 최대 <C> = {land.max():.4f} at γ={G[ig]:.3f}, β={B[ib]:.3f}; 근사비 {land.max()/maxcut:.4f}")
qaoa_summary = {}
for p in (1, 2, 3):
    best = None
    for trial in range(20):
        x0 = np.concatenate([rng.uniform(0, np.pi, p), rng.uniform(0, np.pi / 2, p)])
        res = minimize(lambda x: -qaoa_expect(x), x0, method="L-BFGS-B")
        if best is None or res.fun < best.fun:
            best = res
    val = -best.fun
    psi = qaoa_state(best.x[:p], best.x[p:])
    p_opt = float(np.sum(np.abs(psi) ** 2 * (cost_diag == maxcut)))
    qaoa_summary[p] = (val, val / maxcut, p_opt)
    print(f"   p={p}: 최적 <C> = {val:.4f}, 근사비 {val/maxcut:.4f} (고리 공식 (2p+1)/(2p+2) = {(2*p+1)/(2*p+2):.4f}), "
          f"최적해 측정 확률 {p_opt:.3f}, 2큐비트 게이트 = {2*len(edges)*p} (e^{{-iγZZ}} 당 CNOT 2)")

# ---------- 6) 불모 고원: 그래디언트 분산 vs n ----------
def ansatz_state_bp(theta, n, L, rng_local):
    """무작위 회전축(Rx/Ry/Rz) 하드웨어 효율 ansatz (McClean et al. 2018 스타일)."""
    psi = np.zeros(2 ** n, complex); psi[0] = 1.0
    # 초기 Ry(π/4) 로 축 편향 제거
    for q in range(n):
        psi = apply_1q(psi, ry(np.pi / 4), q, n)
    k = 0
    for layer in range(L):
        for q in range(n):
            ax = theta[1][k]
            g = {0: rx, 1: ry, 2: lambda t: np.diag([np.exp(-1j * t / 2), np.exp(1j * t / 2)])}[ax](theta[0][k])
            psi = apply_1q(psi, g, q, n); k += 1
        for q in range(n - 1):
            psi = apply_cz(psi, q, q + 1, n)
    return psi


def grad_var(n, L, obs_diag_fn, samples=100):
    """첫 파라미터(θ_0, Ry 로 고정)에 대한 ∂<O>/∂θ 의 분산."""
    grads = []
    O = obs_diag_fn(n)
    for _ in range(samples):
        th = rng.uniform(0, 2 * np.pi, n * L)
        axes = rng.integers(0, 3, n * L); axes[0] = 1
        def e(t0):
            t = th.copy(); t[0] = t0
            psi = ansatz_state_bp((t, axes), n, L, rng)
            return float(np.real(psi.conj() @ (O * psi)))
        grads.append(0.5 * (e(th[0] + np.pi / 2) - e(th[0] - np.pi / 2)))
    return np.var(grads)


def obs_global(n):   # Z⊗Z⊗...⊗Z (대각)
    bits = ((np.arange(2 ** n)[:, None] >> (n - 1 - np.arange(n))) & 1)
    return np.prod(1 - 2 * bits, axis=1).astype(complex)


def obs_local(n):    # Z_0 Z_1
    bits = ((np.arange(2 ** n)[:, None] >> (n - 1 - np.arange(n))) & 1)
    return ((1 - 2 * bits[:, 0]) * (1 - 2 * bits[:, 1])).astype(complex)


print("6) 불모 고원: Var[∂<O>/∂θ_0] (무작위 파라미터 100개)")
ns = [2, 3, 4, 5, 6, 7, 8]
var_deep_g, var_deep_l, var_shallow_l = [], [], []
for nn in ns:
    vg = grad_var(nn, 2 * nn, obs_global)     # 깊은 회로(L=2n), 전역 관측량
    vl = grad_var(nn, 2 * nn, obs_local)      # 깊은 회로, 국소 관측량
    vs = grad_var(nn, 1, obs_local)           # 얕은 회로(L=1), 국소 관측량
    var_deep_g.append(vg); var_deep_l.append(vl); var_shallow_l.append(vs)
    print(f"   n={nn}: 깊은(L={2*nn})·전역 {vg:.2e} | 깊은·국소 {vl:.2e} | 얕은(L=1)·국소 {vs:.2e}")
slope_g = np.polyfit(ns, np.log2(var_deep_g), 1)[0]
slope_l = np.polyfit(ns, np.log2(var_deep_l), 1)[0]
print(f"   log2 분산 기울기: 깊은·전역 {slope_g:.2f}/큐비트 (2^(-n) 이면 -1), 깊은·국소 {slope_l:.2f}/큐비트; "
      f"얕은·국소는 n 에 거의 무관")

# ---------- 7) 그림 ----------
fig, axs = plt.subplots(2, 2, figsize=(11, 8))
ax = axs[0, 0]
for key in (("Ry+CZ", 3), ("Ry+CX", 1), ("Ry+CX", 2), ("Ry+CX", 3), ("HVA", 2), ("HVA", 4)):
    hist = vqe_hist[key]
    ax.semilogy(np.arange(1, len(hist) + 1), np.array(hist) - E0 + 1e-12, label=(f"{key[0]} p={key[1]}" if key[0] == "HVA" else f"{key[0]} L={key[1]}"))
ax.set_xlabel("L-BFGS 반복"); ax.set_ylabel("E - E0"); ax.set_title("VQE 수렴 (TFIM n=4)"); ax.legend()
ax = axs[0, 1]
im = ax.imshow(land, origin="lower", extent=[G[0], G[-1], B[0], B[-1]], aspect="auto", cmap="viridis")
ax.plot(G[ig], B[ib], "r*", ms=12)
ax.set_xlabel("gamma"); ax.set_ylabel("beta"); ax.set_title("QAOA p=1 지형 <C> (6-노드 고리)")
fig.colorbar(im, ax=ax)
ax = axs[1, 0]
ax.loglog(shot_list, shot_std, "o-", label="측정 std")
ax.loglog(shot_list, shot_std[0] * np.sqrt(shot_list[0] / np.array(shot_list)), "--", label="1/sqrt(N)")
ax.set_xlabel("shots (기저당)"); ax.set_ylabel("에너지 추정 std"); ax.set_title("샷 잡음"); ax.legend()
ax = axs[1, 1]
ax.semilogy(ns, var_deep_g, "o-", label="깊은(L=2n)·전역 Z^n")
ax.semilogy(ns, var_deep_l, "s-", label="깊은(L=2n)·국소 Z0Z1")
ax.semilogy(ns, var_shallow_l, "^-", label="얕은(L=1)·국소 Z0Z1")
ax.semilogy(ns, [2.0 ** (-k) for k in ns], "k:", label="2^(-n)")
ax.set_xlabel("n"); ax.set_ylabel("Var[dE/dθ0]"); ax.set_title("불모 고원"); ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig("ch13_results.png", dpi=130)
print("7) saved ch13_results.png")

# ── 게이트 비용 기록 (2큐비트 게이트 / 깊이 / T-count) ──
# 위에서 상태벡터에 직접 건 ansatz 들을 게이트 리스트로 다시 짜서 세고, 같은 상태를 내는지(전역 위상 제외) 확인한다.
# 전역 rng 를 건드리지 않도록 자체 난수 발생기를 쓴다.
# 규약(tools/gate_cost.py): 기저 CX + 1큐비트 게이트, 전결합. CZ = 2큐비트 게이트 1개, exp(−iθ Z⊗Z) = 2 CX + Rz(2θ).
# 깊이 = 게이트를 가능한 한 앞에 놓은(ASAP) 층 수. π/4 의 배수가 아닌 회전은 개수 R 로 세고 T ≈ 3·log2(1/ε)·R (ε = 1e-6).
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tools"))
import gate_cost as gc
_g = gc.G
_rng = np.random.default_rng(13)
def hea_list(theta, n_, L_, ent, order=None):                      # ansatz_state(ent='cz') / ansatz_state_cx(ent='cx') 와 같은 순서
    g, k = [], 0
    for layer in range(L_ + 1):
        for q in range(n_): g.append(_g("ry", q, angle=theta[k])); k += 1
        if layer < L_: g += [_g(ent, q, q + 1) for q in (order or range(n_ - 1))]
    return g
def hva_list(theta, n_, p_, order):                                # |+>^n 뒤 [e^{-iγ H_ZZ} (ZZ 결합을 order 순으로), Rx(-2β)] × p
    g = [_g("h", q) for q in range(n_)]
    for l in range(p_):
        g += [_g("rzz", i, i + 1, angle=-theta[2 * l]) for i in order] + [_g("rx", q, angle=-2 * theta[2 * l + 1]) for q in range(n_)]
    return g
def qaoa_list(gammas, betas, edge_order):                          # |+>^n 뒤 [e^{-iγC}: 변마다 e^{+iγ Z_i Z_j/2} (전역 위상 제외), Rx(2β)] × p
    g = [_g("h", q) for q in range(nq)]
    for ga, be in zip(gammas, betas):
        g += [_g("rzz", i, j, angle=-ga / 2) for (i, j) in edge_order] + [_g("rx", q, angle=2 * be) for q in range(nq)]
    return g
def _state(g, n_): return gc.apply(g, np.eye(2 ** n_, dtype=complex)[0], n_)
print("8) 게이트 비용 기록 — 기저 CX + 1큐비트, 전결합; CZ = 1, e^{−iθZZ} = 2 CX + Rz; 깊이 = ASAP 층 수; 각도는 일반값 → 회전 수 R, T ≈ 3·log2(1/ε)·R, ε = 1e-6")
for kind, ent, ref in (("Ry+CZ", "cz", ansatz_state), ("Ry+CX", "cx", ansatz_state_cx)):
    for L_ in (1, 2, 3, 4):
        th_ = _rng.uniform(0, 2 * np.pi, n_params(n, L_)); gl = hea_list(th_, n, L_, ent)
        ok_ = gc.same_up_to_phase(_state(gl, n), ref(th_, n, L_)) and gc.same_up_to_phase(_state(hea_list(vqe_summary[(kind, L_)][3], n, L_, ent), n), ref(vqe_summary[(kind, L_)][3], n, L_)); assert ok_
        assert gc.two_qubit_count(gl) == (n - 1) * L_ and gc.depth(gl) == 3 * L_ + 2 and gc.t_count(gl) == (0, n * (L_ + 1))
        extra = ""
        if ent == "cz":                                            # CZ 끼리는 교환 → 01·23·12 순으로 놓아도 같은 상태
            gb = hea_list(th_, n, L_, ent, order=[0, 2, 1]); assert gc.same_up_to_phase(_state(gb, n), ref(th_, n, L_)); extra = f"; CZ 를 01·23·12 순으로 놓으면 깊이 {gc.depth(gb)}"
        print(f"   {kind} L={L_} (사슬 01·12·23): {gc.fmt(gl)}{extra}  [상태 = {ref.__name__}: {ok_}]")
assert all(gc.depth(hea_list(_rng.uniform(0, 2 * np.pi, n_ * (L_ + 1)), n_, L_, "cx")) == n_ + 3 * L_ - 2 for n_ in (3, 4, 5, 6, 8) for L_ in (1, 2, 3, 5))
print("     → 3) 의 '깊이≈2L+1' 은 얽힘 사슬 하나를 한 층으로 센 층 수다. 게이트 단위 ASAP 깊이는 사슬이 계단으로 놓여 3L+2 (n = 4; 일반 n 은 n + 3L − 2).")
for p_ in (1, 2, 3, 4):
    th_ = _rng.uniform(0, 2 * np.pi, 2 * p_); gb = hva_list(th_, n, p_, [0, 2, 1]); gch = hva_list(th_, n, p_, [0, 1, 2])
    ok_ = all(gc.same_up_to_phase(_state(x, n), hva_state(th_, n, p_, zz_diag)) for x in (gb, gch)) and gc.same_up_to_phase(_state(hva_list(vqe_summary[("HVA", p_)][3], n, p_, [0, 2, 1]), n), hva_state(vqe_summary[("HVA", p_)][3], n, p_, zz_diag)); assert ok_
    assert gc.two_qubit_count(gb) == 2 * (n - 1) * p_ and gc.t_count(gb) == (0, (2 * n - 1) * p_)
    print(f"   HVA p={p_} (H 4 + [ZZ 01·23·12, Rx 4] × p): {gc.fmt(gb)}; ZZ 를 01·12·23 순으로 놓으면 깊이 {gc.depth(gch)}  [상태 = hva_state: {ok_}]")
ring_brick = [(0, 1), (2, 3), (4, 5), (1, 2), (3, 4), (5, 0)]
for p_ in (1, 2, 3):
    ga_, be_ = _rng.uniform(0, np.pi, p_), _rng.uniform(0, np.pi / 2, p_); gb = qaoa_list(ga_, be_, ring_brick); gch = qaoa_list(ga_, be_, edges)
    ok_ = all(gc.same_up_to_phase(_state(x, nq), qaoa_state(ga_, be_)) for x in (gb, gch)); assert ok_
    assert gc.two_qubit_count(gb) == 2 * len(edges) * p_ and gc.t_count(gb) == (0, 2 * nq * p_)
    print(f"   QAOA 6-노드 고리 p={p_} (H 6 + [ZZ 6: 01·23·45 뒤 12·34·50, Rx 6] × p): {gc.fmt(gb)}; 변을 01·12·…·50 순으로 놓으면 깊이 {gc.depth(gch)}  [상태 = qaoa_state: {ok_}]")
