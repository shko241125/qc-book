"""16장 원시 요소 사전 — numpy 검증 스크립트
블록 인코딩(LCU) · 큐비트화(qubitization) · QSP/QSVT · 진폭 추정(QAE, MLAE) · Szegedy 양자 워크를
행렬 수준에서 직접 만들어 정의대로 동작하는지 확인한다.
규약: big-endian (보조 레지스터가 앞, 시스템이 뒤; 각 레지스터 안에서 q0 = 최상위 비트), seed=7.
"""
import os, time
import numpy as np
from numpy.polynomial import chebyshev as C
from scipy.linalg import expm
from scipy.optimize import least_squares
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
matplotlib.rcParams["font.family"] = ["Noto Sans CJK JP", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

rng = np.random.default_rng(7)
T0 = time.time()
I2 = np.eye(2); X = np.array([[0, 1], [1, 0]], float); Z = np.diag([1.0, -1.0])
PAULI = {"I": I2, "X": X, "Z": Z}


def pauli_string(s):
    M = np.array([[1.0]])
    for c in s: M = np.kron(M, PAULI[c])
    return M


# ---------- 1) 블록 인코딩: Pauli 합의 LCU ----------
def tfim_terms(n, J=1.0, h=1.0):
    """열린 사슬 TFIM H = −J Σ Z_i Z_{i+1} − h Σ X_i (12~14장과 같은 모형). (계수, Pauli 문자열) 목록."""
    terms = []
    for i in range(n - 1):
        s = ["I"] * n; s[i] = s[i + 1] = "Z"; terms.append((-J, "".join(s)))
    for i in range(n):
        s = ["I"] * n; s[i] = "X"; terms.append((-h, "".join(s)))
    return terms


def householder_first_column(v):
    """첫 열이 v(실수, 단위 벡터)인 실직교 행렬 — PREP 의 행렬 표현(실제 회로는 상태 준비 회로)."""
    d = len(v); e = np.zeros(d); e[0] = 1.0
    u = v - e
    if np.linalg.norm(u) < 1e-14: return np.eye(d)
    u /= np.linalg.norm(u)
    return np.eye(d) - 2 * np.outer(u, u)


def lcu_block_encoding(terms, n):
    """U = (PREP† ⊗ I) · SELECT · (PREP ⊗ I).  ⟨0|_a U |0⟩_a = H/α,  α = Σ|c_i|.
    PREP|0⟩ = Σ √(|c_i|/α)|i⟩, SELECT = Σ |i⟩⟨i| ⊗ sgn(c_i) P_i (쓰지 않는 색인은 I)."""
    L = len(terms); m = int(np.ceil(np.log2(L))); da = 2 ** m; ds = 2 ** n
    coef = np.array([c for c, _ in terms]); alpha = np.abs(coef).sum()
    amp = np.zeros(da); amp[:L] = np.sqrt(np.abs(coef) / alpha)
    PREP = householder_first_column(amp)
    SEL = np.zeros((da * ds, da * ds))
    for i in range(da):
        Pi = np.sign(coef[i]) * pauli_string(terms[i][1]) if i < L else np.eye(ds)
        SEL[i * ds:(i + 1) * ds, i * ds:(i + 1) * ds] = Pi
    U = np.kron(PREP.T, np.eye(ds)) @ SEL @ np.kron(PREP, np.eye(ds))
    return U, alpha, m


def hamiltonian(terms, n):
    return sum(c * pauli_string(s) for c, s in terms)


n_sys = 4
terms = tfim_terms(n_sys)
Hm = hamiltonian(terms, n_sys)
U, alpha, m_anc = lcu_block_encoding(terms, n_sys)
ds = 2 ** n_sys; da = 2 ** m_anc
evals, evecs = np.linalg.eigh(Hm)
print("1) 블록 인코딩 (LCU), TFIM n = 4, J = h = 1")
print(f"   항 L = {len(terms)}, 보조 큐비트 m = {m_anc}, α = Σ|c_i| = {alpha:.1f}, ‖H‖ = {np.abs(evals).max():.4f}, E0 = {evals[0]:.6f}")
print(f"   U 유니터리 오차 {np.abs(U.T @ U - np.eye(da * ds)).max():.1e}, U 에르미트 오차 {np.abs(U - U.T).max():.1e}")
blk = U[:ds, :ds]
print(f"   ‖⟨0|U|0⟩ − H/α‖ = {np.abs(blk - Hm / alpha).max():.1e}  (부분 정규화 α/‖H‖ = {alpha/np.abs(evals).max():.3f})")
psi0 = evecs[:, 0]
succ = np.linalg.norm(blk @ psi0) ** 2
print(f"   바닥상태에 U 를 걸고 보조 = 0 을 얻을 확률 = (E0/α)² = {succ:.4f}  (= {(evals[0]/alpha)**2:.4f})")

# ---------- 2) 큐비트화: W = (2Π − I)·U 의 고유위상 = ±arccos(λ/α) ----------
Pi = np.zeros(da * ds); Pi[:ds] = 1.0; R = np.diag(2 * Pi - 1)       # 2Π − I (보조가 |0…0⟩ 일 때 +1)
W = R @ U
wev = np.linalg.eigvals(W)
ph = np.angle(wev)
pred = np.concatenate([np.arccos(evals / alpha), -np.arccos(evals / alpha)])
# W 의 고유값 중 시스템 고유값 λ_k 에 해당하는 2ds 개가 ±arccos(λ/α); 나머지(da−2)·ds 개는 ±1 근처의 '쓰지 않는' 공간
match = [np.min(np.abs(ph - p)) for p in pred]
print("\n2) 큐비트화 (qubitization)")
print(f"   W = (2Π−I)U 유니터리, 고유값 {len(wev)}개. 예측 위상 ±arccos(λ_k/α) {len(pred)}개 각각과 가장 가까운 고유위상의 최대 차 = {max(match):.1e}")
th0 = np.arccos(evals[0] / alpha)
print(f"   바닥상태: λ0/α = {evals[0]/alpha:.4f} → 위상 ±{th0:.4f} rad. W 로 QPE 를 하면 E0 = α·cos θ")
for mbits in (6, 8, 10, 12):
    th_hat = np.round(th0 / (2 * np.pi) * 2 ** mbits) / 2 ** mbits * 2 * np.pi
    print(f"     QPE {mbits:2d}비트(W 호출 {2**mbits - 1:5d}번): θ̂ = {th_hat:.5f} → Ê0 = {alpha*np.cos(th_hat):.5f}, 오차 {abs(alpha*np.cos(th_hat)-evals[0]):.1e} (상한 α·π/2^m·sin θ ≈ {alpha*np.pi/2**mbits*np.sin(th0):.1e})")
# 큐비트화 vs U 그대로: U 는 에르미트·유니터리라 고유값이 ±1 뿐 — 위상 정보가 없다
print(f"   비교: U 자체의 고유값은 {np.unique(np.round(np.linalg.eigvalsh(U), 8))} 뿐 (반사) — 반사 R 을 곱해야 회전이 된다")

# ---------- 3) QSP 와 QSVT ----------
def Rsig(x):
    s = np.sqrt(max(0.0, 1 - x * x)); return np.array([[x, s], [s, -x]])


def ez(phi): return np.diag([np.exp(1j * phi), np.exp(-1j * phi)])


def qsp_poly(phis, x):
    """반사 규약 QSP: P(x) = ⟨0| e^{iφ0 Z} Π_k [R(x) e^{iφk Z}] |0⟩ (d = len(phis) − 1 차). x 는 배열 가능."""
    x = np.atleast_1d(np.asarray(x, float)); s = np.sqrt(np.clip(1 - x * x, 0, None))
    # 행 벡터 ⟨0| 를 오른쪽으로 밀며 곱한다: v ← v·e^{iφZ}, v ← v·R(x)
    v0 = np.exp(1j * phis[0]) * np.ones_like(x, complex); v1 = np.zeros_like(x, complex)
    for p in phis[1:]:
        v0, v1 = v0 * x + v1 * s, v0 * s - v1 * x
        v0, v1 = v0 * np.exp(1j * p), v1 * np.exp(-1j * p)
    return v0


xs = np.linspace(-1, 1, 401)


def fit_qsp(f_target, d):
    """Re P(x) = f(x) (f 는 x 의 함수, 차수 d·패리티 d 인 다항식) 를 만족하는 위상 φ_0..φ_d 를 찾는다.
    Dong–Meng–Whaley–Lin 2021: 대칭 위상(φ_k = φ_{d−k}), 초기값 (π/4, 0, …, 0, π/4), 양의 체비쇼프 점에서 Levenberg–Marquardt.
    점 수 = 미지수 수로 먼저 풀고, 수렴이 덜 되면 점을 두 배로 늘려 다시 푼다."""
    nfree = (d + 2) // 2
    full = lambda v: np.concatenate([v, v[::-1][(d + 1) % 2:]])
    best = None
    for k in (nfree, 2 * nfree):
        nodes = cheb_nodes_pos(k); fv = f_target(nodes)
        v0 = np.zeros(nfree); v0[0] = np.pi / 4
        sol = least_squares(lambda v: qsp_poly(full(v), nodes).real - fv, v0, method="lm",
                            xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=20000)
        xx = np.linspace(-1, 1, 2001)
        err = np.abs(qsp_poly(full(sol.x), xx).real - f_target(xx)).max()
        if best is None or err < best[1]: best = (full(sol.x), err)
        if err < 1e-12: break
    return best


def cheb_nodes_pos(k): return np.cos((2 * np.arange(1, k + 1) - 1) * np.pi / (4 * k))   # (0,1] 의 체비쇼프 점


# (b) 부호 함수 근사: f(x) = 0.8·erf(κx) 의 체비쇼프 전개를 홀수 차수 d 에서 자른 다항식
def odd_cheb_target(fun, d):
    c = C.chebinterpolate(fun, d)
    c[0::2] = 0.0                       # 홀수 함수 → 짝수 차 계수 0
    return c


print("\n3) QSP / QSVT")
kappa = 6.0
sign_rows = []
for d in (11, 21, 31):
    c = odd_cheb_target(lambda x: 0.8 * np.vectorize(__import__("math").erf)(kappa * x), d)
    phis, err_fit = fit_qsp(lambda x, c=c: C.chebval(x, c), d)
    Pre = qsp_poly(phis, xs).real
    gap = 0.25
    out = np.abs(xs) >= gap
    err_sign = np.abs(Pre[out] - 0.8 * np.sign(xs[out])).max()
    sign_rows.append((d, phis, Pre, err_fit, err_sign))
    print(f"   부호 근사 d = {d:2d}: 위상 {len(phis)}개 맞춤 잔차 {err_fit:.1e}, |x| ≥ {gap} 에서 |Re P − 0.8·sign| 최대 {err_sign:.3f}")


# (c) QSVT: 같은 위상을 블록 인코딩 U 에 쓰면 Π U_Φ Π = P(H/α)  (U 가 에르미트 → 각 고유값 λ 의 2차원 부분공간에서 U = R(λ))
def qsvt(U, phis, Pdiag):
    """U_Φ = e^{iφ0(2Π−I)} Π_k [U e^{iφk(2Π−I)}]. 실제 회로에서 e^{iφ(2Π−I)} 는 보조 1개 + 다중 제어 NOT 2개 + Rz."""
    D = lambda p: np.diag(np.exp(1j * p * Pdiag))
    M = D(phis[0])
    for p in phis[1:]: M = M @ U @ D(p)
    return M


Pd = 2 * Pi - 1
d_sig, ph_sig, _, _, _ = sign_rows[1]
UP = qsvt(U, ph_sig, Pd); UM = qsvt(U, -ph_sig, Pd)
lam = evals / alpha
target = evecs @ np.diag(qsp_poly(ph_sig, lam)) @ evecs.T
print(f"   QSVT (d = {d_sig}, 부호 근사 위상) on H/α: ‖Π U_Φ Π − P(H/α)‖ = {np.abs(UP[:ds, :ds] - target).max():.1e}")
re_blk = 0.5 * (UP + UM)[:ds, :ds]
print(f"   Re 부분 = (U_Φ + U_−Φ)/2 의 블록 (보조 1개 더, LCU): ‖블록 − Re P(H/α)‖ = {np.abs(re_blk - evecs @ np.diag(qsp_poly(ph_sig, lam).real) @ evecs.T).max():.1e}")

# (c') 규약 변환 검산: R(x) = −i e^{iπ/4 Z} W_x e^{iπ/4 Z} → 반사 규약 위상 Φ 를 W_x 규약으로: 양 끝 +π/4, 안쪽 +π/2, 전역 (−i)^d
def qsp_wx(phis, x):
    s_ = np.sqrt(1 - x * x); Wx = np.array([[x, 1j * s_], [1j * s_, x]])
    Mw = ez(phis[0])
    for p_ in phis[1:]: Mw = Mw @ Wx @ ez(p_)
    return Mw[0, 0]
ph_r = np.random.default_rng(99).normal(size=8); dconv = len(ph_r) - 1
ph_w = ph_r.copy(); ph_w[[0, -1]] += np.pi / 4; ph_w[1:-1] += np.pi / 2
conv_err = max(abs(qsp_poly(ph_r, x)[0] - (-1j) ** dconv * qsp_wx(ph_w, x)) for x in np.linspace(-1, 1, 41))
print(f"   규약 변환 검산(무작위 위상 8개, d = {dconv}): 반사 규약 P 와 (−i)^d × W_x 규약 P 의 최대 차 {conv_err:.1e}")

# (d) 해밀토니안 시뮬레이션: cos(τx), sin(τx) (τ = αt) 의 Jacobi–Anger 전개 → QSVT
from scipy.special import jv


def ja_degree(tau, eps):
    """Jacobi–Anger 를 차수 d 에서 자른 꼬리 2Σ_{k>d}|J_k(τ)| < ε 인 최소 d (cos·sin 공통)."""
    ks = np.arange(0, int(2 * tau + 60)); J = np.abs(jv(ks, tau))
    tail = 2 * np.cumsum(J[::-1])[::-1]                # tail[j] = 2Σ_{k≥j}|J_k|
    return int(np.argmax(tail < eps)) - 1               # 첫 j 에서 하나 뺀 것이 최소 차수


def cos_sin_cheb(tau, d):
    """cos(τx) = J0 + 2Σ(−1)^k J_{2k} T_{2k},  sin(τx) = 2Σ(−1)^k J_{2k+1} T_{2k+1}."""
    cc = np.zeros(d + 1); cs = np.zeros(d + 1)
    for k in range(d + 1):
        if k % 2 == 0: cc[k] = (1 if k == 0 else 2) * (-1) ** (k // 2) * jv(k, tau)
        else: cs[k] = 2 * (-1) ** (k // 2) * jv(k, tau)
    return cc, cs


t_sim = 0.5; tau = alpha * t_sim; eps = 1e-8
dd = ja_degree(tau, eps); dc = dd - dd % 2; dsn = dd - 1 + dd % 2            # cos 는 짝수 차, sin 은 홀수 차 (둘 다 ≤ dd)
cc, _ = cos_sin_cheb(tau, dc); _, cs = cos_sin_cheb(tau, dsn)
ph_c, ec = fit_qsp(lambda x: 0.5 * C.chebval(x, cc), dc)
ph_s, es = fit_qsp(lambda x: 0.5 * C.chebval(x, cs), dsn)
blk_c = 0.5 * (qsvt(U, ph_c, Pd) + qsvt(U, -ph_c, Pd))[:ds, :ds]
blk_s = 0.5 * (qsvt(U, ph_s, Pd) + qsvt(U, -ph_s, Pd))[:ds, :ds]
Uexact = expm(-1j * Hm * t_sim)
err_sim = np.abs(2 * (blk_c - 1j * blk_s) - Uexact).max()
print(f"   해밀토니안 시뮬레이션 t = {t_sim}: τ = αt = {tau:.1f}, ε = {eps:.0e} → 차수 cos {dc}, sin {dsn} (위상 맞춤 잔차 {ec:.1e}, {es:.1e})")
print(f"   2·(블록_cos − i·블록_sin) 과 e^(−iHt) 의 최대 차 = {err_sim:.1e}  (블록 인코딩 U 호출 {dc} + {dsn}회, 부분 정규화 1/2 → OAA 로 회복)")

print("   비용 비교 (TFIM n=4): QSVT 차수(블록 인코딩 호출) vs 2차 Trotter 단계 수")


def trotter2_steps(t, eps):
    """2차(대칭) Trotter: e^{−iA t/2r} e^{−iB t/r} e^{−iA t/2r} 의 r 거듭제곱, 연산자 노름 오차 ≤ ε 인 최소 r (2배 탐색 + 이분)."""
    A = sum(c * pauli_string(s) for c, s in terms if "Z" in s); B = sum(c * pauli_string(s) for c, s in terms if "X" in s)
    Ue = expm(-1j * Hm * t)
    def err(r):
        S = expm(-1j * A * t / (2 * r)) @ expm(-1j * B * t / r) @ expm(-1j * A * t / (2 * r))
        return np.linalg.norm(np.linalg.matrix_power(S, r) - Ue, 2)
    hi = 1
    while err(hi) > eps: hi *= 2
    lo = hi // 2
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if err(mid) > eps: lo = mid
        else: hi = mid
    return hi


nA = sum(1 for _, s_ in terms if "Z" in s_)
cost_rows = []
for eps_c in (1e-3, 1e-6):
    for t in (0.5, 1, 2, 4, 8, 16):
        dq = ja_degree(alpha * t, eps_c)
        r2 = trotter2_steps(t, eps_c)
        cost_rows.append((eps_c, t, alpha * t, dq, r2))
        print(f"     ε = {eps_c:.0e}, t = {t:4}: αt = {alpha*t:5.1f}, QSVT 차수 d = {dq:3d} (cos+sin 블록 호출 ≈ {2*dq:4d}) | Trotter-2 단계 r = {r2:5d} (Pauli 지수 회전 {r2*len(terms)+nA:6d})")

# ---------- 4) 진폭 추정 (QAE) ----------
print("\n4) 진폭 추정: a = ⟨ψ|Π_good|ψ⟩ = sin²θ")


def qae_dist(theta, m):
    """정준 QAE (Brassard–Høyer–Mosca–Tapp 2002): Q 의 고유위상 ±2θ 에 m 비트 QPE. 출력 y 의 확률 (해석식)."""
    M = 2 ** m; y = np.arange(M)
    def F(phi):                      # |(1/M) Σ_k e^{2πik(φ − y/M)}|²
        delta = phi - y / M
        num = np.sin(np.pi * M * delta) ** 2; den = (M * np.sin(np.pi * delta)) ** 2
        out = np.where(np.abs(den) < 1e-30, 1.0, num / np.where(den == 0, 1, den))
        return out
    return 0.5 * (F(theta / np.pi) + F(-theta / np.pi))


# (a) 회로 수준 검증: 3큐비트 A (무작위 실직교), '좋은' 상태 = 마지막 큐비트가 1
nq = 3; dA = 2 ** nq
Q_, _ = np.linalg.qr(rng.normal(size=(dA, dA))); A = Q_
good = np.array([(k & 1) == 1 for k in range(dA)], float)
psi = A[:, 0]; a_true = float(np.sum(good * psi ** 2)); theta = np.arcsin(np.sqrt(a_true))
S0 = np.eye(dA) - 2 * np.outer(np.eye(dA)[0], np.eye(dA)[0]); Schi = np.diag(1 - 2 * good)
Qop = -A @ S0 @ A.T @ Schi                                             # Grover 반복 Q = −A S0 A† Sχ
m = 5; M = 2 ** m
# 게이트 수준 QPE: 계수 레지스터 q0..q_{m−1}(빅 엔디안, q0 = 최상위) + 시스템. 상태를 (2,)*m + (dA,) 텐서로 다룬다.
def apply_1q(T, G, q):                                 # 계수 큐비트 q 에 2×2 게이트
    return np.moveaxis(np.tensordot(G, T, axes=([1], [q])), 0, q)
def apply_cphase(T, phi, a, b):                        # 계수 큐비트 a, b 사이 제어 위상 diag(1,1,1,e^{iφ})
    T = T.copy(); idx = [slice(None)] * T.ndim; idx[a] = 1; idx[b] = 1
    T[tuple(idx)] *= np.exp(1j * phi); return T
def apply_cU(T, Umat, c):                              # 계수 큐비트 c 가 1 일 때 시스템(마지막 축)에 Umat
    T = T.copy(); idx = [slice(None)] * T.ndim; idx[c] = 1
    T[tuple(idx)] = np.tensordot(T[tuple(idx)], Umat.T, axes=([T.ndim - 2], [0])); return T
def apply_swap(T, a, b): return np.swapaxes(T, a, b)
Hg = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
T = np.zeros((2,) * m + (dA,), complex); T[(0,) * m] = psi
for q in range(m): T = apply_1q(T, Hg, q)
for jq in range(m):                                    # q_j 가 Q^{2^{m−1−j}} 를 제어 (N&C 5.2 그림)
    T = apply_cU(T, np.linalg.matrix_power(Qop, 2 ** (m - 1 - jq)), jq)
# 교과서 역 QFT: 스왑 먼저, 그다음 QFT 게이트를 역순·켤레 위상으로
for q in range(m // 2): T = apply_swap(T, q, m - 1 - q)
for jq in reversed(range(m)):
    for k in reversed(range(jq + 1, m)):
        T = apply_cphase(T, -2 * np.pi / 2 ** (k - jq + 1), jq, k)
    T = apply_1q(T, Hg, jq)
p_circ = np.sum(np.abs(T.reshape(M, dA)) ** 2, axis=1)
# 스왑을 빼먹으면? (비트 반전 — 흔한 실수)
T2 = np.zeros((2,) * m + (dA,), complex); T2[(0,) * m] = psi
for q in range(m): T2 = apply_1q(T2, Hg, q)
for jq in range(m): T2 = apply_cU(T2, np.linalg.matrix_power(Qop, 2 ** (m - 1 - jq)), jq)
for jq in reversed(range(m)):
    for k in reversed(range(jq + 1, m)):
        T2 = apply_cphase(T2, -2 * np.pi / 2 ** (k - jq + 1), jq, k)
    T2 = apply_1q(T2, Hg, jq)
p_noswap = np.sum(np.abs(T2.reshape(M, dA)) ** 2, axis=1)
p_anal = qae_dist(theta, m)
print(f"   게이트 수준(3큐비트 A, m = {m}, 제어 Q^(2^(m−1−j)) + 교과서 역 QFT): a = {a_true:.4f}, θ = {theta:.4f}; 해석식과 최대 차 {np.abs(p_circ - p_anal).max():.1e}")
print(f"   역 QFT 의 스왑을 빼면: 해석식과 최대 차 {np.abs(p_noswap - p_anal).max():.2f}, 최빈 y = {np.argmax(p_noswap)} (비트 반전이 아니라 분포가 망가짐)")
yb = np.argmax(p_anal); print(f"   최빈 출력 y = {yb} 와 거울상 y = {M - yb} (같은 â = sin²(πy/M) = {np.sin(np.pi*yb/M)**2:.4f}), 두 값의 확률 합 {p_anal[yb] + p_anal[M - yb]:.3f}")

# (b) 오차 스케일링: 정준 QAE (QPE) vs 고전 몬테카를로 vs MLAE (Suzuki et al. 2020)
a0 = 0.3; th0a = np.arcsin(np.sqrt(a0))
ae_rows = []
q81 = 8 / np.pi ** 2                                                       # BHMT 정리 12 의 성공 확률
z81 = 1.3126                                                               # 정규분포 양측 81% 분위 (Φ⁻¹(0.905))
for m in range(3, 13):
    M = 2 ** m; pd_ = qae_dist(th0a, m); est = np.sin(np.pi * np.arange(M) / M) ** 2
    err = np.abs(est - a0); order = np.argsort(err); cum = np.cumsum(pd_[order])
    e81 = err[order][np.searchsorted(cum, q81)]                             # 81% 분위 오차
    bound = 2 * np.pi * np.sqrt(a0 * (1 - a0)) / M + np.pi ** 2 / M ** 2   # BHMT 상한
    p_in = pd_[err <= bound + 1e-15].sum()
    rmse_q = np.sqrt(np.sum(pd_ * err ** 2))
    nAc = 2 * (M - 1) + 1                                                  # A 호출 수 (Q 하나 = A 둘)
    e81_mc = z81 * np.sqrt(a0 * (1 - a0) / nAc)
    ae_rows.append((m, nAc, e81, bound, p_in, rmse_q, e81_mc))
    print(f"   QAE m = {m:2d}: A 호출 {nAc:5d} | 81% 분위 오차 {e81:.2e} (BHMT 상한 {bound:.2e}, 상한 안 확률 {p_in:.3f}) | RMSE {rmse_q:.2e} | 같은 A 호출 몬테카를로 81% 오차 {e81_mc:.2e}")
print("   → QAE 한 번의 RMSE 는 QPE 분포의 1/k² 꼬리 때문에 1/M 로 줄지 않는다: 여러 번 돌려 중앙값을 쓰는 것이 표준")
# (연습문제 3) a 를 (0,1) 에서 고르게 평균한 한 번 실행 RMSE × √M — 1/√M 이면 m 에 무관한 상수
avec = np.linspace(0.005, 0.995, 400)
for m in (4, 6, 8, 10, 12):
    M = 2 ** m; est = np.sin(np.pi * np.arange(M) / M) ** 2
    mse = np.mean([np.sum(qae_dist(np.arcsin(np.sqrt(a_)), m) * (est - a_) ** 2) for a_ in avec])
    print(f"   (연습문제 3) m = {m:2d}: a 평균 RMSE·√M = {np.sqrt(mse) * np.sqrt(M):.3f}")
fr = [(m, (2 ** m * th0a / np.pi) % 1) for m in (4, 12)]
print("   a = 0.3 에서 f = Mθ/π 의 소수 부분: " + ", ".join(f"m={m_}: f={f_:.3f}, sin²(πf)={np.sin(np.pi*f_)**2:.3f}" for m_, f_ in fr))
m7 = 7; p7 = qae_dist(th0a, m7)


def mlae(theta, Ks, shots, rng, grid=np.linspace(1e-4, np.pi / 2 - 1e-4, 20001)):
    """MLAE: Q^k A|0⟩ 를 k ∈ Ks 마다 shots 번 측정해 좋은 결과 수 h_k 를 얻고, 우도 최대 θ 를 격자 탐색."""
    ll = np.zeros_like(grid)
    for k in Ks:
        pk = np.sin((2 * k + 1) * theta) ** 2
        h = rng.binomial(shots, pk)
        pg = np.clip(np.sin((2 * k + 1) * grid) ** 2, 1e-15, 1 - 1e-15)
        ll += h * np.log(pg) + (shots - h) * np.log(1 - pg)
    return np.sin(grid[np.argmax(ll)]) ** 2


mlae_rows = []
shots = 100
for Kmax in range(0, 8):
    Ks = [0] + [2 ** j for j in range(Kmax)]
    calls = shots * sum(Ks)                                               # Q 호출 총수 (측정당 k 번)
    ests = np.array([mlae(th0a, Ks, shots, rng) for _ in range(200)])
    rmse = np.sqrt(np.mean((ests - a0) ** 2))
    Nq = shots * sum(2 * k + 1 for k in Ks)                               # A 호출 총수 (k 번 Q = 2k+1 번 A)
    mlae_rows.append((len(Ks), calls, Nq, rmse))
    print(f"   MLAE k = {','.join(map(str, Ks)):>22s} × {shots}샷: A 호출 {Nq:6d}, RMSE {rmse:.2e} (같은 A 호출의 몬테카를로 {np.sqrt(a0*(1-a0)/Nq):.2e}, 비 {np.sqrt(a0*(1-a0)/Nq)/rmse:.1f})")

# ---------- 5) Szegedy 양자 워크: 위상 간격 ≈ 2√(2δ) (고전 스펙트럼 간격 δ 의 제곱근) ----------
print("\n5) Szegedy 양자 워크")


def szegedy(P):
    """P: 행 확률 행렬 (가역). |ψ_x⟩ = |x⟩ ⊗ Σ_y √P_xy |y⟩, Π_A = Σ|ψ_x⟩⟨ψ_x|, S = SWAP, W = S(2Π_A − I)S(2Π_A − I)."""
    N = P.shape[0]
    Psi = np.zeros((N * N, N))
    for x in range(N): Psi[x * N:(x + 1) * N, x] = np.sqrt(P[x])
    RA = 2 * Psi @ Psi.T - np.eye(N * N)
    S = np.zeros((N * N, N * N))
    for x in range(N):
        for y in range(N): S[y * N + x, x * N + y] = 1
    return S @ RA @ S @ RA, Psi


sz_rows = []
for N in (8, 12, 16, 24, 32):
    P = np.zeros((N, N))
    for x in range(N):                                                  # 게으른(lazy) 고리 위 랜덤 워크
        P[x, x] = 0.5; P[x, (x + 1) % N] += 0.25; P[x, (x - 1) % N] += 0.25
    lamP = np.sort(np.linalg.eigvalsh(P))[::-1]                          # 대칭이라 판별 행렬 D = P
    delta = 1 - lamP[1]
    Wz, Psi = szegedy(P)
    ev = np.linalg.eigvals(Wz); ang = np.abs(np.angle(ev))
    # '바쁜' 부분공간(span{Psi, S Psi})의 위상: ±2 arccos(λ_k(D)); λ=1 은 위상 0
    pred = 2 * np.arccos(np.clip(lamP, -1, 1))
    nz = np.sort(ang[ang > 1e-9])
    Delta = nz[0]
    sz_rows.append((N, delta, Delta, 2 * np.sqrt(2 * delta)))
    ok = max(np.min(np.abs(ang - p)) for p in pred)
    print(f"   고리 N = {N:2d}: 고전 간격 δ = {delta:.5f} (혼합 ~1/δ = {1/delta:7.1f} 걸음) | 양자 위상 간격 Δ = {Delta:.5f} (~1/Δ = {1/Delta:5.1f}), 2√(2δ) = {2*np.sqrt(2*delta):.5f}, Δ/√δ = {Delta/np.sqrt(delta):.3f} | 예측 위상 일치 {ok:.1e}")
print(f"   N 이 커지면 Δ/√δ → 2√2 = {2*np.sqrt(2):.3f}: 양자 워크는 고전 혼합 시간의 제곱근 만큼의 단계로 위상 간격을 분해한다")

# ---------- 그림 ----------
fig, axs = plt.subplots(2, 3, figsize=(15, 8.2))
ax = axs[0, 0]; xx = np.linspace(-1, 1, 400)
ax.plot(xx, np.arccos(xx), color="#868E96", lw=1, label="θ = arccos(λ/α)")
ax.plot(evals / alpha, np.arccos(evals / alpha), "o", color="#4C6EF5", ms=5, label="TFIM n=4 고유값 (16개)")
ax.axvline(-np.abs(evals).max() / alpha, ls=":", color="#F03E3E", lw=1); ax.axvline(np.abs(evals).max() / alpha, ls=":", color="#F03E3E", lw=1, label="±‖H‖/α (부분 정규화 1.47)")
ax.set_xlabel("λ/α"); ax.set_ylabel("W 의 고유위상 θ (rad)"); ax.set_title("(1) 큐비트화: W = (2Π−I)U 의 위상"); ax.legend(fontsize=7)

ax = axs[0, 1]
ax.plot(xs, 0.8 * np.sign(xs), color="#868E96", lw=1, label="0.8·sign(x)")
for (d_, ph_, Pre_, _, _), c_ in zip(sign_rows, ("#74C0FC", "#4C6EF5", "#1C2D7A")):
    ax.plot(xs, Pre_, color=c_, lw=1.3, label=f"QSP Re P(x), d = {d_}")
ax.axvspan(-0.25, 0.25, color="#F1F3F5"); ax.set_xlabel("x"); ax.set_ylabel("Re P(x)"); ax.set_title("(2) QSP: 위상 d+1 개로 만든 부호 근사"); ax.legend(fontsize=7)

ax = axs[0, 2]
for eps_c, mk in ((1e-3, "o"), (1e-6, "s")):
    rows = [r for r in cost_rows if r[0] == eps_c]
    ts = [r[1] for r in rows]
    ax.loglog(ts, [2 * r[3] for r in rows], mk + "-", color="#4C6EF5", label=f"QSVT 블록 호출 2d, ε={eps_c:.0e}")
    ax.loglog(ts, [r[4] * len(terms) + nA for r in rows], mk + "--", color="#F03E3E", label=f"Trotter-2 Pauli 회전, ε={eps_c:.0e}")
ax.set_xlabel("시뮬레이션 시간 t"); ax.set_ylabel("호출 / 회전 수"); ax.set_title(r"(3) $e^{-iHt}$: QSVT vs Trotter (TFIM n=4)"); ax.legend(fontsize=7)

ax = axs[1, 0]
nAq = np.array([r[1] for r in ae_rows])
ax.loglog(nAq, [r[2] for r in ae_rows], "o-", color="#4C6EF5", label="정준 QAE 81% 분위 오차")
ax.loglog(nAq, [r[5] for r in ae_rows], "o:", color="#74C0FC", label="정준 QAE RMSE (꼬리)")
ax.loglog([r[2] for r in mlae_rows], [r[3] for r in mlae_rows], "s-", color="#2F9E44", label="MLAE RMSE (200회)")
NN = np.logspace(1, 4.5, 50)
ax.loglog(NN, np.sqrt(a0 * (1 - a0) / NN), "k--", lw=1, label="몬테카를로 √(a(1−a)/N)")
ax.loglog(NN, 2 * np.pi * np.sqrt(a0 * (1 - a0)) / NN, color="#868E96", lw=1, label="1/N 기울기 기준선")
ax.set_xlabel("A 호출 수 N"); ax.set_ylabel("추정 오차"); ax.set_title("(4) 진폭 추정 a = 0.3: 1/N vs 1/√N"); ax.legend(fontsize=7)

ax = axs[1, 1]
ax.semilogy(np.arange(len(p7)), np.maximum(p7, 1e-8), "o-", color="#4C6EF5", ms=3)
yt = th0a / np.pi * len(p7); ax.axvline(yt, color="#F03E3E", ls=":", label=f"참값 Mθ/π = {yt:.2f} (a = 0.3)")
ax.axvline(len(p7) - yt, color="#F03E3E", ls=":")
yy = np.arange(int(np.ceil(yt)), int(np.ceil(yt)) + 49); ax.semilogy(yy, np.sin(np.pi * (yt % 1)) ** 2 / (2 * np.pi ** 2 * (yy - yt) ** 2), "k--", lw=0.8, label="sin²(πf)/(2π²(y−Mθ/π)²)")
ax.set_xlabel("QPE 출력 y (m = 7)"); ax.set_ylabel("확률 (로그)"); ax.set_title("(4b) QAE 출력 분포: 1/k² 꼬리"); ax.legend(fontsize=7)

ax = axs[1, 2]
Ns = np.array([r[0] for r in sz_rows])
ax.loglog(Ns, [r[1] for r in sz_rows], "o-", color="#868E96", label="고전 스펙트럼 간격 δ (∝ N⁻²)")
ax.loglog(Ns, [r[2] for r in sz_rows], "s-", color="#4C6EF5", label="Szegedy 위상 간격 Δ (∝ N⁻¹)")
ax.set_xticks(Ns); ax.set_xticklabels([str(v) for v in Ns]); ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_xlabel("고리 크기 N"); ax.set_ylabel("간격"); ax.set_title("(5) 양자 워크: Δ ≈ 2√(2δ)"); ax.legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)), "ch16_results.png"), dpi=120)
print("\n6) saved ch16_results.png")
print(f"\n   (실행 시간 {time.time()-T0:.1f} s)")
