"""18장 탐색 도구와 반례 — numpy 검증 스크립트
17장 체크리스트의 각 단계를 자동화하는 세 도구를 작은 크기에서 직접 만들어 돌린다.
(A) 자동 회로 합성: 구조(CNOT 배치) 탐색 + 연속 파라미터 최적화로 목표 유니터리를 만든다.
    2큐비트 표적의 최소 CNOT 수(0~3)를 Shende–Bullock–Markov 불변량과 대조하고,
    CCZ(= Toffoli)가 전결합에서 6개, 선형 연결에서 8개를 필요로 함을 전수 탐색으로 확인한다.
(B) 변분 탐색: ADAPT-VQE가 연산자 풀에서 회로를 한 게이트씩 키운다. 풀을 잘못 고르면 멈춘다(기울기 0 함정).
(C) dequantization 검사: 표본·질의 접근(ℓ² 표본)만으로 저계수 행렬의 특이 부분공간과 추천 행을 얻는다
    (Frieze–Kannan–Vempala 표본화, Tang 2019). 읽는 성분 수가 행렬 크기에 무관함과, 안정 계수가 크면 깨짐을 본다.
(D) 자유 페르미온 검사: 파울리 항을 Jordan–Wigner로 옮겨 Majorana 차수를 세고, 2차면 2n×2n 행렬로 동역학을 정확히 푼다.
    17장의 TFIM 급랭이 "MPS로는 어렵지만 고전으로는 쉬운" 예임을 n = 1000에서 확인한다.
규약: big-endian(큐비트 0 = 가장 왼쪽 텐서 인자). 난수는 seed 고정.
"""
import itertools, os, time
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import svds, expm_multiply
from scipy.optimize import minimize
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
matplotlib.rcParams["font.family"] = ["Noto Sans CJK JP", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

T0 = time.time()
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1, -1]).astype(complex)


def kron_all(ms):
    out = np.eye(1, dtype=complex)
    for m in ms: out = np.kron(out, m)
    return out


def op1(n, q, m):
    return kron_all([m if i == q else I2 for i in range(n)])


def cx(n, c, t):
    P0, P1 = np.diag([1, 0]).astype(complex), np.diag([0, 1]).astype(complex)
    return op1(n, c, P0) + op1(n, c, P1) @ op1(n, t, X)


# =====================================================================
# (A) 자동 회로 합성
# =====================================================================
class Template:
    """회로 틀: 모든 큐비트에 U3, 그 뒤 층마다 CNOT(a,b) + U3(a) + U3(b).  U3 = Rz·Ry·Rz (각 3개).
    R_P(θ) = exp(−iθP/2).  CNOT의 방향은 주변 U3가 흡수하므로(H⊗H 켤레) 쌍 {a,b}만 구조를 정한다."""

    def __init__(self, n, pairs):
        self.n, self.d, self.pairs = n, 2 ** n, list(pairs)
        self.ops = []                                   # ('r', 생성원 행렬) 또는 ('c', CNOT 행렬)
        def u3(q):
            for g in (Z, Y, Z): self.ops.append(("r", op1(n, q, g)))
        for q in range(n): u3(q)
        for a, b in self.pairs:
            self.ops.append(("c", cx(n, a, b))); u3(a); u3(b)
        self.npar = sum(1 for kind, _ in self.ops if kind == "r")
        self.Id = np.eye(self.d, dtype=complex)

    def unitary(self, th):
        V, k = self.Id.copy(), 0
        for kind, m in self.ops:
            if kind == "c": V = m @ V
            else:
                V = np.cos(th[k] / 2) * V - 1j * np.sin(th[k] / 2) * (m @ V); k += 1
        return V

    def cost_grad(self, th, Ud):
        """f = 1 − |Tr(U†V)|²/d² (전역 위상 무시 공정 비충실도)와 해석적 기울기(앞·뒤 부분곱)."""
        pre, V, k = [], self.Id.copy(), 0
        for kind, m in self.ops:
            if kind == "c": V = m @ V
            else:
                V = np.cos(th[k] / 2) * V - 1j * np.sin(th[k] / 2) * (m @ V); k += 1
                pre.append(V)
        t = np.trace(Ud @ V); f = 1 - abs(t) ** 2 / self.d ** 2
        g, B, k = np.zeros(self.npar), Ud.copy(), self.npar - 1
        for kind, m in reversed(self.ops):
            if kind == "c": B = B @ m
            else:
                gk = -0.5j * np.trace(B @ m @ pre[k])    # ∂Tr/∂θ_k,  B = U†·(k 뒤의 게이트들)
                g[k] = -2 * np.real(np.conj(t) * gk) / self.d ** 2
                B = np.cos(th[k] / 2) * B - 1j * np.sin(th[k] / 2) * (B @ m); k -= 1
        return f, g


def instantiate(tmpl, U, rng, starts=4, tol=1e-12):
    """구조가 정해진 틀의 각도를 L-BFGS로 맞춘다(무작위 시작 여러 번). 반환: 최소 비충실도, 각도."""
    Ud, best = U.conj().T, (1.0, None)
    for _ in range(starts):
        r = minimize(tmpl.cost_grad, rng.uniform(0, 2 * np.pi, tmpl.npar), args=(Ud,), jac=True,
                     method="L-BFGS-B", options=dict(maxiter=3000, ftol=1e-16, gtol=1e-10))
        if r.fun < best[0]: best = (r.fun, r.x)
        if best[0] < tol: break
    return max(best[0], 0.0), best[1]


def haar(d, rng):
    G = (rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))) / np.sqrt(2)
    Q, R = np.linalg.qr(G)
    return Q * (np.diag(R) / abs(np.diag(R)))


def sbm_cnot_count(U):
    """Shende–Bullock–Markov(2004) 판정: γ(U) = U (Y⊗Y) Uᵀ (Y⊗Y), U ∈ SU(4).
    γ = ±I → 0개,  tr γ = 0 이고 γ² = −I → 1개,  tr γ 가 실수 → 2개,  그 밖 → 3개."""
    U = U / np.linalg.det(U) ** 0.25
    YY = np.kron(Y, Y); g = U @ YY @ U.T @ YY; tr = np.trace(g)
    if np.allclose(g, np.eye(4), atol=1e-9) or np.allclose(g, -np.eye(4), atol=1e-9): return 0, tr
    if abs(tr) < 1e-9 and np.allclose(g @ g, -np.eye(4), atol=1e-9): return 1, tr
    if abs(tr.imag) < 1e-9: return 2, tr
    return 3, tr


rng = np.random.default_rng(18)
print("=" * 78); print("(A) 자동 회로 합성 — 구조 탐색 + 각도 최적화"); print("=" * 78)

# A0. 기울기 검산
T_ = Template(3, [(0, 1), (1, 2), (0, 2)]); U_ = haar(8, rng); th_ = rng.uniform(0, 6, T_.npar)
f_, g_ = T_.cost_grad(th_, U_.conj().T); e_ = 1e-6
gn_ = np.array([(T_.cost_grad(th_ + e_ * np.eye(T_.npar)[i], U_.conj().T)[0]
                 - T_.cost_grad(th_ - e_ * np.eye(T_.npar)[i], U_.conj().T)[0]) / (2 * e_) for i in range(T_.npar)])
print(f"A0) 해석적 기울기 vs 중심차분 (3큐비트, 파라미터 {T_.npar}개): 최대 차 {np.abs(g_-gn_).max():.1e}")

# A1. 2큐비트 표적: CNOT k개로 닿는 최소 비충실도
SWAP = np.eye(4)[[0, 2, 1, 3]].astype(complex)
CZ = np.diag([1, 1, 1, -1]).astype(complex)
CP = np.diag([1, 1, 1, 1j])                                            # 제어 위상 π/2
ISWAP = np.array([[1, 0, 0, 0], [0, 0, 1j, 0], [0, 1j, 0, 0], [0, 0, 0, 1]])
w_, v_ = np.linalg.eigh(SWAP); SQSWAP = (v_ * np.sqrt(w_.astype(complex))) @ v_.conj().T
QFT2 = np.array([[1, 1, 1, 1], [1, 1j, -1, -1j], [1, -1, 1, -1], [1, -1j, -1, 1j]]) / 2
targets2 = [("CZ", CZ), ("CP(π/2)", CP), ("iSWAP", ISWAP), ("√SWAP", SQSWAP), ("SWAP", SWAP), ("QFT₂", QFT2)]
print("\nA1) 2큐비트 표적: CNOT k개 틀로 닿는 최소 비충실도 1 − |Tr(U†V)|²/16   (시작점 6개)")
print("    표적       k=0      k=1      k=2      k=3     수치 최소 k   SBM 판정")
tab2 = {}
for name, U in targets2:
    row = [instantiate(Template(2, [(0, 1)] * k), U, rng, starts=6)[0] for k in range(4)]
    kmin = next(k for k, f in enumerate(row) if f < 1e-10)
    sbm, _ = sbm_cnot_count(U); tab2[name] = row
    print(f"    {name:8s} " + " ".join(f"{f:8.4f}" for f in row) + f"      {kmin}            {sbm}" + ("" if kmin == sbm else "   ← 불일치"))

# A2. Haar 무작위 SU(4) 20개
print("\nA2) Haar 무작위 2큐비트 유니터리 20개")
fl2, fl3, imtr = [], [], []
for _ in range(20):
    U = haar(4, rng)
    fl2.append(instantiate(Template(2, [(0, 1)] * 2), U, rng, starts=6)[0])
    fl3.append(instantiate(Template(2, [(0, 1)] * 3), U, rng, starts=6)[0])
    c_, tr_ = sbm_cnot_count(U); imtr.append(abs(tr_.imag))
fl2, fl3, imtr = np.array(fl2), np.array(fl3), np.array(imtr)
print(f"    CNOT 2개: 비충실도 최소 {fl2.min():.1e}, 중앙값 {np.median(fl2):.1e}, 최대 {fl2.max():.1e}  (닿은 표적 {int((fl2<1e-10).sum())}/20)")
print(f"    CNOT 3개: 비충실도 최대 {fl3.max():.1e}  (닿은 표적 {int((fl3<1e-10).sum())}/20)")
print(f"    k=2 바닥과 |Im tr γ(U)| 의 순위 상관(Spearman): {np.corrcoef(np.argsort(np.argsort(fl2)), np.argsort(np.argsort(imtr)))[0,1]:.3f}"
      f"   (최소 바닥 표적의 |Im tr γ| = {imtr[np.argmin(fl2)]:.1e})")
print("    파라미터 세기: SU(4) 15차원, 틀의 유효 자유도 6 + 4k → k ≥ ⌈9/4⌉ = 3;  일반 n: ⌈(4ⁿ−3n−1)/4⌉ = "
      + ", ".join(f"n={n}: {int(np.ceil((4**n-3*n-1)/4))}" for n in (2, 3, 4, 5)))

# A3. 3큐비트 CCZ(= H 켤레로 Toffoli): 구조 전수 탐색
CCZ = np.diag([1] * 7 + [-1]).astype(complex)


def canon(seq, perms):
    """표적이 큐비트 순열과 전치에 대칭이면(CCZ) 순열·뒤집기로 같은 구조는 한 번만 본다."""
    cands = []
    for p in perms:
        s = tuple(tuple(sorted((p[a], p[b]))) for a, b in seq)
        cands += [s, s[::-1]]
    return min(cands)


def exhaustive(U, pairs, perms, kmax, starts, rng):
    out = []
    for k in range(kmax + 1):
        classes = sorted({canon(s, perms) for s in itertools.product(pairs, repeat=k)})
        t0 = time.time(); fs = []
        for s in classes:
            fs.append(instantiate(Template(3, list(s)), U, rng, starts=starts)[0])
        fs = np.array(fs); solved = [classes[i] for i in np.where(fs < 1e-10)[0]]
        out.append((k, len(classes), fs.min(), solved, time.time() - t0))
        print(f"    CNOT {k}개: 구조 {3 if len(pairs)==3 else 2}^{k} = {len(pairs)**k:4d}가지 → 대칭 묶음 {len(classes):3d}개, "
              f"최소 비충실도 {fs.min():.3e}, 닿은 구조 {len(solved)}개  ({time.time()-t0:.1f} s)", flush=True)
        if solved: break
    return out


print("\nA3) CCZ 합성 — 전결합(쌍 3종): CNOT 수를 늘려 가며 모든 구조를 시도 (구조당 시작점 4개)")
res_all = exhaustive(CCZ, [(0, 1), (0, 2), (1, 2)], list(itertools.permutations(range(3))), 6, 4, rng)
for s in res_all[-1][3]: print("      닿은 구조:", " ".join(f"{a}{b}" for a, b in s))
print("    선형 연결 0–1–2(쌍 2종):")
res_line = exhaustive(CCZ, [(0, 1), (1, 2)], [(0, 1, 2), (2, 1, 0)], 8, 4, rng)
for s in res_line[-1][3]: print("      닿은 구조:", " ".join(f"{a}{b}" for a, b in s))
print(f"    바닥값 검산: 1 − (6/8)² = {1-(6/8)**2:.4f} (항등), 3/8 = {3/8:.4f}, (2−√2)/4 = {(2-np.sqrt(2))/4:.4f}")
tA = time.time() - T0

# =====================================================================
# (B) 변분 탐색 — ADAPT-VQE
# =====================================================================
print("\n" + "=" * 78); print("(B) 변분 탐색 — ADAPT-VQE로 회로를 한 게이트씩 키우기"); print("=" * 78)
PM = {"I": sp.identity(2, dtype=complex, format="csr"), "X": sp.csr_matrix(X), "Y": sp.csr_matrix(Y), "Z": sp.csr_matrix(Z)}


def pauli(n, spec):
    """spec = {큐비트: 'X'|'Y'|'Z'} 인 파울리 문자열의 희소 행렬 (big-endian)."""
    out = sp.identity(1, dtype=complex, format="csr")
    for q in range(n): out = sp.kron(out, PM[spec.get(q, "I")], format="csr")
    return out


def tfim(n, J=1.0, h=1.0):
    """열린 사슬 TFIM H = −J Σ Z_i Z_{i+1} − h Σ X_i (12~17장과 같은 모형)."""
    H = sp.csr_matrix((2 ** n, 2 ** n), dtype=complex)
    for i in range(n - 1): H = H - J * pauli(n, {i: "Z", i + 1: "Z"})
    for i in range(n): H = H - h * pauli(n, {i: "X"})
    return H


def apply_rots(th, ops, psi0):
    psi = psi0.copy()
    for t, P in zip(th, ops): psi = np.cos(t / 2) * psi - 1j * np.sin(t / 2) * (P @ psi)
    return psi


def energy_grad(th, ops, psi0, H):
    """E(θ) = ⟨ψ|H|ψ⟩, |ψ⟩ = Π_k exp(−iθ_k P_k/2)|ψ0⟩ 와 기울기(수반 방법, 게이트 수에 선형)."""
    psi = apply_rots(th, ops, psi0); lam = H @ psi
    E, g = np.real(np.vdot(psi, lam)), np.zeros(len(th))
    for k in range(len(th) - 1, -1, -1):
        P, t = ops[k], th[k]
        g[k] = 2 * np.real(np.vdot(lam, -0.5j * (P @ psi)))
        psi = np.cos(t / 2) * psi + 1j * np.sin(t / 2) * (P @ psi)
        lam = np.cos(t / 2) * lam + 1j * np.sin(t / 2) * (P @ lam)
    return E, g


def make_pool(n, kind):
    """풀 원소 = (이름, 행렬, CNOT 수). Y가 홀수 개인 파울리 문자열은 실수 진폭을 보존한다(qubit-ADAPT 풀)."""
    pool = []
    if kind in ("Y", "full"):
        pool += [(f"Y{i}", pauli(n, {i: "Y"}), 0) for i in range(n)]
    if kind in ("YZ", "full"):
        pool += [(f"Y{i}Z{j}", pauli(n, {i: "Y", j: "Z"}), 2) for i in range(n) for j in range(n) if i != j]
    if kind == "full":
        pool += [(f"Y{i}X{j}", pauli(n, {i: "Y", j: "X"}), 2) for i in range(n) for j in range(n) if i != j]
    return pool


def adapt(H, psi0, pool, E0, max_ops, gtol=1e-6, etol=1e-10):
    """ADAPT-VQE: (1) 풀 전체의 기울기 |∂E/∂θ_P|_{θ=0} = |⟨ψ|i[H,P]|ψ⟩|/… 측정, (2) 최대인 것을 붙임,
    (3) 모든 각도를 다시 최적화. 기울기 최대값 < gtol 이면 멈춘다."""
    chosen, th, psi, hist, nfev, ngrad = [], np.zeros(0), psi0.copy(), [], 0, 0
    while True:
        lam = H @ psi; E = np.real(np.vdot(psi, lam))
        grads = np.array([2 * np.real(np.vdot(lam, -0.5j * (P @ psi))) for _, P, _ in pool])
        j = int(np.argmax(abs(grads))); gmax = abs(grads[j]); ngrad += len(pool)
        hist.append((len(chosen), E - E0, gmax, sum(pool[c][2] for c in chosen)))
        if gmax < gtol or E - E0 < etol or len(chosen) >= max_ops: break
        chosen.append(j); th = np.append(th, 0.0)
        ops = [pool[c][1] for c in chosen]
        r = minimize(energy_grad, th, args=(ops, psi0, H), jac=True, method="BFGS", options=dict(gtol=1e-9, maxiter=2000))
        th, psi, nfev = r.x, apply_rots(r.x, ops, psi0), nfev + r.nfev
    return hist, [pool[c][0] for c in chosen], (ngrad, nfev)


def first_below(hist, tol):
    for nops, err, _, ncx in hist:
        if err < tol: return nops, ncx
    return None


def hva_setup(n):
    """HVA 층 exp(−iβ ΣX)·exp(−iγ ΣZZ) 를 빠르게: ΣZZ는 대각, ΣX = H^{⊗n}·(대각)·H^{⊗n}."""
    bits = (np.arange(2 ** n)[:, None] >> np.arange(n - 1, -1, -1)) & 1; s = 1 - 2 * bits
    dz = (s[:, :-1] * s[:, 1:]).sum(1).astype(float); dx = s.sum(1).astype(float)
    Hn = kron_all([(X + Z) / np.sqrt(2)] * n).real
    return dz, dx, Hn


def hva_energy_grad(par, dz, dx, Hn, psi0, Hd):
    p = len(par) // 2; psi = psi0.copy(); states = []
    for l in range(p):
        psi = np.exp(-1j * par[2 * l] * dz) * psi
        psi = Hn @ (np.exp(-1j * par[2 * l + 1] * dx) * (Hn @ psi)); states.append(psi)
    lam = Hd @ psi; E = np.real(np.vdot(psi, lam)); g = np.zeros(2 * p)
    for l in range(p - 1, -1, -1):
        psi = states[l]
        g[2 * l + 1] = 2 * np.real(np.vdot(lam, -1j * (Hn @ (dx * (Hn @ psi)))))
        psi = Hn @ (np.exp(1j * par[2 * l + 1] * dx) * (Hn @ psi)); lam = Hn @ (np.exp(1j * par[2 * l + 1] * dx) * (Hn @ lam))
        g[2 * l] = 2 * np.real(np.vdot(lam, -1j * dz * psi))
        lam = np.exp(1j * par[2 * l] * dz) * lam
    return E, g


def hva(Hd, n, p, psi0, rng, warm=None, starts=20):
    """13장의 HVA(파라미터 2p개). 시작점: 이전 층 수의 최적값 + 0(warm), 그리고 [−π, π] 무작위 starts개.
    반환: 최소 에너지, 그 각도, 시작점별 최종 에너지, 에너지 평가 횟수."""
    dz, dx, Hn = hva_setup(n); inits = [rng.uniform(-np.pi, np.pi, 2 * p) for _ in range(starts)]
    if warm is not None: inits.insert(0, np.concatenate([warm, [0.0, 0.0]]))
    finals, best, nfev = [], (np.inf, None), 0
    for x0 in inits:
        r = minimize(hva_energy_grad, x0, args=(dz, dx, Hn, psi0, Hd), jac=True, method="BFGS", options=dict(gtol=1e-10, maxiter=5000))
        finals.append(r.fun); nfev += r.nfev
        if r.fun < best[0]: best = (r.fun, r.x)
    return best[0], best[1], np.array(finals), nfev


adapt_curves, adapt_cost, hva_pts, trap_rows = {}, {}, {}, []
print("B1) TFIM(J = h = 1) 바닥상태, 시작 |+…+⟩, 풀 {Y_i Z_j} (n(n−1)개, 각 CNOT 2개)")
print("    n   풀 크기  E0           오차<1e-3: (게이트, CNOT)   오차<1e-6          오차<1e-8          마지막 (게이트, 오차)   탐색 비용")
for n in (4, 6, 8):
    H = tfim(n); E0 = np.linalg.eigvalsh(H.toarray())[0]
    plus = np.ones(2 ** n, dtype=complex) / np.sqrt(2 ** n)
    hist, names, (ngrad, nfev) = adapt(H, plus, make_pool(n, "YZ"), E0, max_ops=60)
    adapt_curves[n] = hist; adapt_cost[n] = (ngrad, nfev)
    print(f"    {n}   {n*(n-1):4d}    {E0:.6f}   {str(first_below(hist,1e-3)):22s}  {str(first_below(hist,1e-6)):17s}  {str(first_below(hist,1e-8)):17s}"
          f"  ({hist[-1][0]}, {abs(hist[-1][1]):.1e})   풀 기울기 측정 {ngrad}회, 에너지·기울기 평가 {nfev}회")
    if n == 4: print("        n = 4에서 고른 순서:", " ".join(names))

print("\nB2) 풀과 시작 상태에 따른 함정 (n = 4, 6): 멈춘 뒤의 오차와 최대 기울기")
print("    n  시작      풀                    풀 크기  붙인 게이트  멈춘 오차 E−E0   최대 기울기   설명되는 값")
for n in (4, 6):
    H = tfim(n); ev = np.linalg.eigvalsh(H.toarray()); E0 = ev[0]
    plus = np.ones(2 ** n, dtype=complex) / np.sqrt(2 ** n); zero = np.zeros(2 ** n, dtype=complex); zero[0] = 1
    # 설명용 기준값: (a) |+…+⟩의 에너지, (b) 최선의 실수 곱상태(평균장), (c) 패리티 ΠX 두 구역 바닥의 평균
    e_plus = np.real(np.vdot(plus, H @ plus)) - E0
    ry = [pauli(n, {i: "Y"}) for i in range(n)]
    e_mf = min(minimize(energy_grad, rng.uniform(0, np.pi, n), args=(ry, zero, H), jac=True, method="BFGS").fun for _ in range(10)) - E0
    Px = pauli(n, {i: "X" for i in range(n)}).toarray().real; w, v = np.linalg.eigh(H.toarray())
    par = np.real(np.einsum("ij,ij->j", v.conj(), Px @ v)); e_odd = w[np.where(par < -0.5)[0][0]]
    e_sym = (e_odd - E0) / 2
    refs = {("|+…+⟩", "Y"): f"|+…+⟩ 자체 {e_plus:.3f}", ("|0…0⟩", "Y"): f"평균장 최적 {e_mf:.3f}",
            ("|0…0⟩", "YZ"): f"(E₁−E₀)/2 = {e_sym:.3f}"}
    for iname, psi0, kind in (("|+…+⟩", plus, "Y"), ("|+…+⟩", plus, "YZ"), ("|0…0⟩", zero, "Y"), ("|0…0⟩", zero, "YZ"), ("|0…0⟩", zero, "full")):
        pool = make_pool(n, kind); hist, _, _ = adapt(H, psi0, pool, E0, max_ops=60)
        trap_rows.append((n, iname, kind, len(pool), hist[-1][0], hist[-1][1], hist[-1][2], hist))
        label = {"Y": "{Y_i}", "YZ": "{Y_iZ_j}", "full": "{Y_i, Y_iZ_j, Y_iX_j}"}[kind]
        print(f"    {n}  {iname}   {label:22s} {len(pool):3d}      {hist[-1][0]:3d}        {abs(hist[-1][1]):9.2e}      {hist[-1][2]:.1e}     {refs.get((iname, kind), '바닥상태 도달')}")
    print(f"       (n = {n}: 패리티 ΠX 홀 구역의 최저 에너지 E₁ = {e_odd:.6f}, 간극 E₁−E₀ = {e_odd-E0:.4f} = 4·sin(π/(4n+2)) = {4*np.sin(np.pi/(4*n+2)):.4f})")

print("\nB3) 비교 기준 — 구조를 아는 틀 HVA(13장): 층 p, 파라미터 2p, CNOT 2p(n−1). 시작점 = 이전 최적 + 무작위 20개")
print("    n    p   파라미터   CNOT    최선 오차 E−E0   최선의 2배 이내로 간 시작점 비율   에너지·기울기 평가(누적)")
for n in (4, 6):
    Hd = tfim(n).toarray(); E0 = np.linalg.eigvalsh(Hd)[0]
    plus = np.ones(2 ** n, dtype=complex) / np.sqrt(2 ** n); hva_pts[n] = []; warm, tot = None, 0
    for p in range(1, 17):
        Eb, warm, finals, nfev = hva(Hd, n, p, plus, rng, warm=warm); err = max(Eb - E0, 0.0); tot += nfev
        frac = np.mean(finals - E0 < 2 * err + 1e-9)
        hva_pts[n].append((2 * p, max(err, 1e-16), 2 * p * (n - 1)))
        print(f"    {n}   {p:2d}     {2*p:3d}     {2*p*(n-1):3d}      {err:9.2e}            {frac:4.2f}                         {tot}")
        if err < 1e-9: break
tB = time.time() - T0 - tA

# =====================================================================
# (C) dequantization 검사 — ℓ² 표본만으로 저계수 근사
# =====================================================================
print("\n" + "=" * 78); print("(C) dequantization 검사 — 표본·질의 접근으로 저계수 근사 (FKV 표본화)"); print("=" * 78)


def lowrank_matrix(m, n, k, noise, rng, spectrum=None):
    """계수 k의 '선호 행렬' + 작은 조밀 잡음. 특이값은 1 → 0.5로 선형 감소(기본)."""
    U = np.linalg.qr(rng.normal(size=(m, k)))[0]; V = np.linalg.qr(rng.normal(size=(n, k)))[0]
    s = np.linspace(1.0, 0.5, k) if spectrum is None else spectrum
    A = (U * s) @ V.T
    if noise > 0: A = A + noise * np.linalg.norm(s) * rng.normal(size=(m, n)) / np.sqrt(m * n)
    return A


class SQAccess:
    """표본·질의(sample-and-query) 접근. 전처리 O(mn) 한 번(QRAM에 데이터를 싣는 비용과 같은 급).
    이후: 행 i를 ‖A_i‖²/‖A‖_F²로, 행 안의 열 j를 A_ij²/‖A_i‖²로 뽑고, 성분·노름을 O(1)에 읽는다. 읽은 성분 수를 센다."""

    def __init__(self, A):
        self.A = A; self.row2 = (A ** 2).sum(1); self.F2 = self.row2.sum(); self.p = self.row2 / self.F2
        self.cum = np.cumsum(A ** 2, axis=1); self.reads = 0

    def sample_rows(self, r, rng): return rng.choice(len(self.p), size=r, p=self.p)

    def sample_cols_in_row(self, i, s, rng):
        return np.minimum(np.searchsorted(self.cum[i], rng.uniform(0, self.row2[i], s)), self.A.shape[1] - 1)

    def entries(self, rows, cols):
        self.reads += len(rows) * len(cols); return self.A[np.ix_(rows, cols)]


def fkv(sq, r, c, k, rng):
    """FKV: 행 r개를 노름² 비례로 뽑아 S(r×n, E[SᵀS] = AᵀA), 그 열 c개를 다시 뽑아 W(r×c). W의 SVD만 계산.
    열 분포는 '표본 행을 균등하게 하나 고르고 그 안에서 성분² 비례'의 혼합 — S 전체를 읽지 않는다.
    반환: 표본 행 번호, 계수 coef(r×k) — 근사 오른쪽 특이벡터 v̂_t = Sᵀ u_t/σ_t = Aᵀ[idx]·coef[:,t]."""
    idx = sq.sample_rows(r, rng); scale = 1 / np.sqrt(r * sq.p[idx])
    pick = rng.integers(0, r, c)
    jdx = np.array([sq.sample_cols_in_row(idx[l], 1, rng)[0] for l in pick])
    Wraw = sq.entries(idx, jdx) * scale[:, None]                         # S의 r×c 부분
    q = (Wraw ** 2 / (sq.row2[idx] * scale ** 2)[:, None]).mean(0)        # 혼합 분포에서의 열 확률 q_j
    W = Wraw / np.sqrt(c * q)[None, :]
    u, sg, _ = np.linalg.svd(W, full_matrices=False)
    return idx, (u[:, :k] / sg[:k]) * scale[:, None], sg[:k]


k = 5
print("C1) 계수 5 + 잡음 행렬: 표본 행·열 수 r(= c)에 따른 초과 오차")
print("    초과 오차 = (‖A − A·V̂V̂ᵀ‖_F² − ‖A − A_k‖_F²)/‖A‖_F²  (5회 평균; 검증을 위해서만 A 전체를 씀)")
r_list = [25, 50, 100, 200, 400, 800]
print("    크기 m×n        안정 계수  최적 꼬리 | " + " ".join(f"r={r:<5d}" for r in r_list) + "| r=800에서 읽은 성분 / mn")
fkv_curves, mats, top = {}, {}, {}
for m, n in [(250, 250), (1000, 1000), (4000, 4000)]:
    A = lowrank_matrix(m, n, k, 0.3, rng); sq = SQAccess(A); mats[m] = (A, sq)
    _, sv, vt = svds(A, k=k); top[m] = vt[np.argsort(-sv)].T              # 위 k개 특이쌍만 (검증용)
    opt = (sq.F2 - (sv ** 2).sum()) / sq.F2; srank = sq.F2 / sv.max() ** 2
    row = []
    for r in r_list:
        ex = []
        for _ in range(5):
            sq.reads = 0
            idx, coef, sg = fkv(sq, r, r, k, rng)
            Q = np.linalg.qr(A[idx].T @ coef)[0]                         # 검증용으로만 v̂를 통째로 만든다
            ex.append((sq.F2 - np.linalg.norm(A @ Q) ** 2) / sq.F2 - opt)
        row.append(np.mean(ex))
    fkv_curves[m] = row
    print(f"    {m:4d}×{n:<4d}      {srank:5.2f}     {opt:.3f}    | " + " ".join(f"{e:7.4f}" for e in row) + f"| {sq.reads/(m*n):.4f}")
sl = np.polyfit(np.log(r_list), np.log(fkv_curves[1000]), 1)[0]
print(f"    1000×1000 곡선의 로그-로그 기울기 {sl:.2f} (초과 오차 ∝ 1/r 이면 −1)")

print("\nC2) 추천 질의: 사용자 i의 계수-5 근사 행 (A_k)_i 를 표본만으로 — 4000×4000, r = c = 400, 내적 표본 s = 400")
A, sq = mats[4000]; m, n = A.shape
Vk = top[4000]
sq.reads = 0
idx, coef, sg = fkv(sq, 400, 400, k, rng); stage1 = sq.reads
users = rng.choice(m, 20, replace=False); rel, prec, qreads = [], [], []
for i in users:
    sq.reads = 0
    jj = sq.sample_cols_in_row(i, 400, rng)                              # j ∝ A_ij²
    Bij = sq.entries(idx, jj); ai = sq.entries([i], jj)[0]
    ip = sq.row2[i] * (Bij / ai[None, :]).mean(1)                        # ⟨A_i, A_l⟩ 의 비편향 추정 (l = 표본 행)
    lam = coef.T @ ip                                                    # λ_t ≈ ⟨A_i, v̂_t⟩
    qreads.append(sq.reads)
    xhat = (A[idx].T @ (coef @ lam))                                     # 검증용: 예측 행 전체 (실제로는 필요한 성분만 r개씩 읽음)
    xtrue = (A[i] @ Vk) @ Vk.T
    rel.append(np.linalg.norm(xhat - xtrue) / np.linalg.norm(xtrue))
    prec.append(len(set(np.argsort(-xhat)[:10]) & set(np.argsort(-xtrue)[:10])) / 10)
print(f"    1단계(부분공간): 성분 {stage1:,}개 읽음 = 전체 {m*n:,}개의 {stage1/(m*n):.2%}")
print(f"    2단계(질의당): 성분 {int(np.mean(qreads)):,}개 읽음 = 전체의 {np.mean(qreads)/(m*n):.2%}")
print(f"    사용자 20명: 상대 오차 ‖x̂ − (A_k)_i‖/‖(A_k)_i‖ 평균 {np.mean(rel):.3f} (최대 {np.max(rel):.3f}), 상위 10개 추천 일치율 평균 {np.mean(prec):.2f}")
# 정확한 내적을 쓰면(표본 오차 제거) 남는 오차 = 부분공간 오차
rel_exact = []
for i in users:
    lam = coef.T @ (A[idx] @ A[i]); xhat = A[idx].T @ (coef @ lam); xtrue = (A[i] @ Vk) @ Vk.T
    rel_exact.append(np.linalg.norm(xhat - xtrue) / np.linalg.norm(xtrue))
print(f"    (내적을 정확히 계산하면 평균 {np.mean(rel_exact):.3f} — 나머지는 내적 표본 오차)")

print("\nC3) 어디서 깨지는가 — 안정 계수 ‖A‖_F²/‖A‖² 의존성 (1000×1000, r = 200, 5회 평균)")
print("    Gram 오차 = ‖AᵀA − SᵀS‖₂/‖A‖²,  대략 (안정 계수)/√r;  10% 오차에 필요한 r ≈ (안정 계수/0.1)²")
print("    행렬                         안정 계수   Gram 오차(측정)  안정계수/√r   10%에 필요한 표본 행 수(상한)   행 수 m = 1000 과 비교")
N = 1000
lap = (2 * np.eye(N) - np.eye(N, k=1) - np.eye(N, k=-1))
cases = [("계수 5 (추천 행렬형)", lowrank_matrix(N, N, 5, 0.0, rng)),
         ("계수 40, 평평한 스펙트럼", lowrank_matrix(N, N, 40, 0.0, rng, spectrum=np.ones(40))),
         ("계수 300, 평평한 스펙트럼", lowrank_matrix(N, N, 300, 0.0, rng, spectrum=np.ones(300))),
         ("3중 대각(희소·최대 계수)", lap)]
srank_rows = []
for name, A_ in cases:
    sq_ = SQAccess(A_); G = A_.T @ A_; nrm2 = np.linalg.eigvalsh(G)[-1]; srank = sq_.F2 / nrm2; errs = []
    for _ in range(5):
        idx_ = sq_.sample_rows(200, rng); S_ = A_[idx_] / np.sqrt(200 * sq_.p[idx_])[:, None]
        errs.append(abs(np.linalg.eigvalsh(G - S_.T @ S_)[[0, -1]]).max() / nrm2)
    need = (srank / 0.1) ** 2
    srank_rows.append((name, srank, np.mean(errs)))
    print(f"    {name:24s}   {srank:7.1f}      {np.mean(errs):6.3f}          {srank/np.sqrt(200):6.2f}        {need:9.2e}        {'행렬의 행 수보다 적음' if need < N else '행렬 전체(1000행)보다 많음'}")
tC = time.time() - T0 - tA - tB

# =====================================================================
# (D) 자유 페르미온 검사 — 얽힘 검사가 놓치는 고전 경계
# =====================================================================
print("\n" + "=" * 78); print("(D) 자유 페르미온(matchgate) 검사 — Majorana 차수 세기와 2n×2n 동역학"); print("=" * 78)


def majorana_degree(spec, n, axis="X"):
    """파울리 문자열(spec = {자리: 'X'|'Y'|'Z'})을 Jordan–Wigner로 옮겼을 때의 Majorana 단항식 차수.
    axis = 점유수 역할을 하는 파울리(끈). axis = 'X'이면 c_{2j} = (Π_{k<j} X_k) Z_j, c_{2j+1} = (Π_{k<j} X_k) Y_j,
    X_j = i c_{2j} c_{2j+1}.  끈이 아닌 파울리(홀)가 오른쪽에 홀수 개 있으면 그 자리의 연산자에 끈이 한 번 곱해진다."""
    deg, odd_right = 0, 0
    for j in range(n - 1, -1, -1):
        pj = spec.get(j, "I")
        is_axis, is_id = (pj == axis), (pj == "I")
        if odd_right:                       # 끈(axis)을 곱한다: I ↔ axis, 나머지 둘은 서로 바뀜(여전히 홀)
            is_axis, is_id = is_id, is_axis
        deg += 0 if is_id else (2 if is_axis else 1)
        if pj not in ("I", axis): odd_right ^= 1
    return deg


def model_terms(name, n):
    T = []
    if name.startswith("TFIM"):
        T += [{i: "Z", i + 1: "Z"} for i in range(n - 1)] + [{i: "X"} for i in range(n)]
        if "세로장" in name: T += [{i: "Z"} for i in range(n)]
        if "다음 이웃" in name: T += [{i: "Z", i + 2: "Z"} for i in range(n - 2)]
    if name == "XY 사슬": T += [{i: a, i + 1: a} for i in range(n - 1) for a in "XY"]
    if name == "Heisenberg 사슬": T += [{i: a, i + 1: a} for i in range(n - 1) for a in "XYZ"]
    if name == "2D TFIM 2×3 (뱀 순서)":
        idx = lambda r, c: r * 3 + (c if r % 2 == 0 else 2 - c)
        bonds = [(idx(r, c), idx(r, c + 1)) for r in range(2) for c in range(2)] + [(idx(0, c), idx(1, c)) for c in range(3)]
        T += [{a: "Z", b: "Z"} for a, b in bonds] + [{i: "X"} for i in range(6)]
    return T


print("D1) 모형별 최대 Majorana 차수 (끈 축 X/Y/Z 중 최선). 2 이하면 자유 페르미온 → 고전 다항 시간")
print("    모형                          항 수   축 X   축 Y   축 Z    판정")
ff_rows = []
for name in ("TFIM 1D", "TFIM 1D + 세로장", "TFIM 1D + 다음 이웃 ZZ", "XY 사슬", "Heisenberg 사슬", "2D TFIM 2×3 (뱀 순서)"):
    T = model_terms(name, 6); degs = [max(majorana_degree(s, 6, ax) for s in T) for ax in "XYZ"]
    ok = min(degs) <= 2; ff_rows.append((name, degs, ok))
    print(f"    {name:26s}   {len(T):3d}     {degs[0]}      {degs[1]}      {degs[2]}     {'자유 페르미온' if ok else '상호작용(이 변환으로는 2차 아님)'}")


def majorana_generator(n, J=1.0, h=1.0):
    """TFIM = (i/4) Σ A_ab c_a c_b 의 실 반대칭 A (2n×2n).  Heisenberg 그림: dc/dt = A c."""
    A = np.zeros((2 * n, 2 * n))
    for j in range(n):     A[2 * j, 2 * j + 1], A[2 * j + 1, 2 * j] = -2 * h, 2 * h           # −h X_j = −h·i c_{2j} c_{2j+1}
    for j in range(n - 1): A[2 * j + 1, 2 * j + 2], A[2 * j + 2, 2 * j + 1] = -2 * J, 2 * J   # −J Z_j Z_{j+1} = −J·i c_{2j+1} c_{2j+2}
    return A


def gaussian_entropy(G, l):
    """공분산 Γ_ab = (i/2)⟨[c_a, c_b]⟩ 의 앞 2l×2l 블록에서 블록 엔트로피(비트).
    고유값은 ±ν 쌍으로 나오고 쌍마다 모드 하나다: 정렬한 위쪽 l개만 센다(ν = 0 인 모드도 정확히 한 번)."""
    nu = np.clip(np.sort(np.linalg.eigvalsh(1j * G[:2 * l, :2 * l]))[l:], 0, 1)
    q = np.clip((1 + nu) / 2, 1e-15, 1 - 1e-15)
    return float(-(q * np.log2(q) + (1 - q) * np.log2(1 - q)).sum())


def ff_quench(n, times):
    """|+…+⟩ (모든 X_j = +1, 가우스 상태)에서 임계 TFIM으로 급랭. 고유분해 한 번 → 임의 시각의 Γ(t) = O Γ0 Oᵀ."""
    A = majorana_generator(n); lam, V = np.linalg.eigh(1j * A)            # A = V (−i lam) V†
    G0 = np.zeros((2 * n, 2 * n))
    for j in range(n): G0[2 * j, 2 * j + 1], G0[2 * j + 1, 2 * j] = 1, -1
    out = []
    for t_ in times:
        O = ((V * np.exp(-1j * lam * t_)) @ V.conj().T).real
        out.append(O @ G0 @ O.T)
    return out


def half_entropy_sv(psi, n):
    s = np.linalg.svd(psi.reshape(2 ** (n // 2), -1), compute_uv=False) ** 2; s = s[s > 1e-15]
    return float(-(s * np.log2(s)).sum())


n = 12; times = [0.5, 1.0, 2.0, 3.0, 4.0]
H12 = tfim(n).tocsc(); plus = np.ones(2 ** n, dtype=complex) / np.sqrt(2 ** n)
Gs = ff_quench(n, times); dS, dX = [], []
for t_, G in zip(times, Gs):
    psi = expm_multiply(-1j * t_ * H12, plus)
    dS.append(abs(half_entropy_sv(psi, n) - gaussian_entropy(G, n // 2)))
    xs = [np.real(np.vdot(psi, pauli(n, {j: "X"}) @ psi)) for j in range(n)]
    dX.append(max(abs(xs[j] - G[2 * j, 2 * j + 1]) for j in range(n)))
print(f"\nD2) 검산 (n = 12, |+…+⟩ → 임계 TFIM, t = {times}): 상태벡터 vs 24×24 행렬")
print(f"    반쪽 엔트로피 최대 차 {max(dS):.1e} 비트,  ⟨X_j⟩ 최대 차 {max(dX):.1e}")

print("\nD3) 같은 급랭을 n = 1000에서 — 2000×2000 행렬만으로")
t1 = time.time(); big_t = [2, 5, 10, 20, 40]; Gb = ff_quench(1000, big_t); ff_big = []
print("    t      반쪽 엔트로피 S (비트)    MPS가 필요로 하는 결합 차원 하한 2^S")
for t_, G in zip(big_t, Gb):
    S_ = gaussian_entropy(G, 500); ff_big.append((t_, S_))
    print(f"    {t_:3d}        {S_:7.2f}                 {2.0**S_:.1e}")
rate = (ff_big[-1][1] - ff_big[1][1]) / (big_t[-1] - big_t[1])
print(f"    성장률 {rate:.3f} 비트/단위시간 (t = 5~40),  계산 시간 {time.time()-t1:.1f} s  — 상태벡터라면 2^1000 차원")

print("\nD4) Wick(가우스성) 검사: 2점 상관만으로 만든 엔트로피 S_가우스 와 실제 S (n = 12, t = 2, 반쪽 블록)")
maj = []
for j in range(n // 2):
    string = {k: "X" for k in range(j)}
    maj += [pauli(n, {**string, j: "Z"}), pauli(n, {**string, j: "Y"})]
print("    모형                       S_실제     S_가우스    차이(비가우스성)")
wick_rows = []
for name, g in (("TFIM (자유)", 0.0), ("TFIM + 세로장 0.5", 0.5), ("TFIM + 세로장 1.0", 1.0)):
    Hg = H12.copy()
    for i in range(n): Hg = Hg - g * pauli(n, {i: "Z"}).tocsc()
    psi = expm_multiply(-2.0j * Hg, plus)
    Gm = np.zeros((n, n))
    for a in range(n):
        va = maj[a] @ psi
        for b in range(a + 1, n):
            Gm[a, b] = np.real(1j * np.vdot(va, maj[b] @ psi)); Gm[b, a] = -Gm[a, b]
    St, Sg = half_entropy_sv(psi, n), gaussian_entropy(Gm, n // 2); wick_rows.append((name, St, Sg))
    print(f"    {name:22s}   {St:7.4f}    {Sg:7.4f}     {Sg-St:8.1e}")

print("\nD5) stoquastic 검사 (n = 6): 계산 기저에서 H의 비대각 원소가 모두 ≤ 0 이면 양자 몬테카를로에 부호 문제가 없다")
print("    기저에 따라 달라지므로, 한쪽 부분격자(홀수 자리)에 Z를 곱하는 회전도 함께 본다")


def ham_from(terms, n):
    H = sp.csr_matrix((2 ** n, 2 ** n), dtype=complex)
    for c, spec in terms: H = H + c * pauli(n, spec)
    return H


def max_offdiag(H):
    M = H.toarray(); off = M - np.diag(np.diag(M))
    return float(off.real.max()), float(abs(off.imag).max())


n6 = 6; Zodd = pauli(n6, {i: "Z" for i in range(1, n6, 2)})
heis = lambda pairs, J: [(J, {a: s, b: s}) for a, b in pairs for s in "XYZ"]
nn = [(i, i + 1) for i in range(n6 - 1)]; nnn = [(i, i + 2) for i in range(n6 - 2)]
sto_cases = [("TFIM 1D (−ZZ − X)", [(-1.0, s) for s in model_terms("TFIM 1D", n6)]),
             ("2D TFIM 2×3", [(-1.0, s) for s in model_terms("2D TFIM 2×3 (뱀 순서)", n6)]),
             ("반강자성 Heisenberg 사슬", heis(nn, 1.0)),
             ("J1–J2 Heisenberg (J2 = 0.5, 좌절)", heis(nn, 1.0) + heis(nnn, 0.5))]
print("    모형                               비대각 최대(원래 기저)   비대각 최대(부분격자 회전 뒤)   판정")
sto_rows = []
for name, terms in sto_cases:
    H_ = ham_from(terms, n6); a, ai = max_offdiag(H_); b, bi = max_offdiag(Zodd @ H_ @ Zodd)
    ok = (a <= 1e-12 and ai < 1e-12) or (b <= 1e-12 and bi < 1e-12); sto_rows.append((name, a, b, ok))
    print(f"    {name:32s}      {a:+.1f}                    {b:+.1f}                 {'stoquastic → QMC 가능' if ok else '부호 문제(이 두 기저에서는)'}")
tD = time.time() - T0 - tA - tB - tC

# =====================================================================
# (E) 연습문제 검산
# =====================================================================
print("\n" + "=" * 78); print("(E) 연습문제 검산"); print("=" * 78)
_, trg = sbm_cnot_count(SWAP); Usw = SWAP / np.linalg.det(SWAP) ** 0.25
gsw = Usw @ np.kron(Y, Y) @ Usw.T @ np.kron(Y, Y)
print(f"E1) γ(SWAP) = ({gsw[0,0]:.3f})·I 인가: {np.allclose(gsw, gsw[0,0]*np.eye(4))},  tr γ = {trg:.3f} → 실수가 아니므로 CNOT 3개")
print("    파라미터 하한 ⌈(4ⁿ−3n−1)/4⌉: " + ", ".join(f"n={n_}: {int(np.ceil((4**n_-3*n_-1)/4))}" for n_ in (2, 3, 4)) + "  (CCZ의 6개는 '일반' 유니터리가 아니라서 모순이 아님)")
for n_ in (4, 6):
    H_ = tfim(n_).toarray(); w_, v_ = np.linalg.eigh(H_)
    Px_ = pauli(n_, {i: "X" for i in range(n_)}).toarray().real
    zero_ = np.zeros(2 ** n_); zero_[0] = 1
    w_even = (1 + zero_ @ Px_ @ zero_) / 2
    par_ = np.real(np.einsum("ij,ij->j", v_.conj(), Px_ @ v_)); e_odd_ = w_[np.where(par_ < -0.5)[0][0]]
    print(f"E2) n = {n_}: |0…0⟩의 패리티 짝/홀 구역 가중치 = {w_even:.2f}/{1-w_even:.2f} → 도달 가능한 최저 에너지 = (E₀+E₁)/2, 오차 하한 (E₁−E₀)/2 = {(e_odd_-w_[0])/2:.4f}")
c_fit = srank_rows[0][2] * np.sqrt(200) / srank_rows[0][1]
t_svd = time.time(); _ = np.linalg.svd(rng.normal(size=(10000, 192)), compute_uv=False); t_svd = time.time() - t_svd
print(f"E3) Gram 오차 ≈ c·(안정 계수)/√r 의 c = {c_fit:.2f} (계수 5 행렬, r = 200 측정) → 안정 계수 3, 오차 0.05에 필요한 r ≈ {(c_fit*3/0.05)**2:.0f}")
r_e3 = (c_fit * 3 / 0.05) ** 2
print(f"    10,000×192 행렬: 표본 행 {r_e3:.0f}개를 통째로 읽으면 성분 {r_e3*192:,.0f}개 = 전체 {10000*192:,}개의 {r_e3*192/(10000*192):.0%}")
print(f"    10,000×192 조밀 행렬의 전체 SVD: {t_svd:.2f} s — 통째로 읽고 푸는 편이 빠르다")

# =====================================================================
# (F) 보강 — 본문 주장의 근거를 따로 확인 (난수는 별도 seed, 위 결과에 영향 없음)
# =====================================================================
print("\n" + "=" * 78); print("(F) 보강 확인"); print("=" * 78)
rngF = np.random.default_rng(1818); tF = time.time()

# F1. 시작점 수에 따른 합성 결과: 닿는 구조를 시작점 4개가 놓치는가, 선형 7개는 시작점을 늘려도 안 닿는가
def count_solved(pairs, perms, k, starts):
    classes = sorted({canon(s, perms) for s in itertools.product(pairs, repeat=k)})
    fs = np.array([instantiate(Template(3, list(s)), CCZ, rngF, starts=starts)[0] for s in classes])
    count_solved.last = [classes[i] for i in np.where(fs < 1e-10)[0]]
    return len(classes), int((fs < 1e-10).sum()), fs.min(), np.sort(np.round(fs, 4))
nc, ns, fmin, fl_ = count_solved([(0, 1), (0, 2), (1, 2)], list(itertools.permutations(range(3))), 6, 12)
print(f"F1) 전결합 CNOT 6개, 구조당 시작점 12개: {nc}묶음 중 {ns}묶음이 닿음 (시작점 4개일 때 {len(res_all[-1][3])}묶음)")
solved12 = list(count_solved.last)
for s_ in solved12: print("      닿은 구조:", " ".join(f"{a}{b}" for a, b in s_))
both = sorted(set(solved12) | set(res_all[-1][3]))
print(f"    두 실행(시작점 4개, 12개)을 합치면 닿은 묶음 {len(both)}개; 4개 실행만 찾은 것 {len(set(res_all[-1][3]) - set(solved12))}개, 12개 실행만 찾은 것 {len(set(solved12) - set(res_all[-1][3]))}개")
nc5, ns5, fmin5, fl5 = count_solved([(0, 1), (0, 2), (1, 2)], list(itertools.permutations(range(3))), 5, 4)
print(f"    전결합 CNOT 5개 구조별 바닥값의 종류: {sorted(set(fl5.tolist()))}  (구조 전체의 최솟값 {fmin5:.4f})")
nc7, ns7, fmin7, _ = count_solved([(0, 1), (1, 2)], [(0, 1, 2), (2, 1, 0)], 7, 12)
print(f"    선형 CNOT 7개, 구조당 시작점 12개: {nc7}묶음 중 {ns7}묶음이 닿음, 최소 비충실도 {fmin7:.4f}")

# F2. Gram 오차 대 표본 행 수 (1000×1000, 3회 평균): 10% 에 닿는 r 을 직접 잰다
print("F2) Gram 오차 ‖AᵀA − SᵀS‖₂/‖A‖² 대 표본 행 수 r (행렬의 행 수는 1000; r > 1000 은 전체를 읽는 것보다 많다)")
r_gram = [200, 1000, 4000, 16000]
print("    행렬                         " + "  ".join(f"r={r:<6d}" for r in r_gram))
gram_tab = []
for name, A_ in cases:
    sq_ = SQAccess(A_); G = A_.T @ A_; nrm2 = np.linalg.eigvalsh(G)[-1]; row = []
    for r in r_gram:
        e = []
        for _ in range(3):
            idx_ = sq_.sample_rows(r, rngF); S_ = A_[idx_] / np.sqrt(r * sq_.p[idx_])[:, None]
            e.append(abs(np.linalg.eigvalsh(G - S_.T @ S_)[[0, -1]]).max() / nrm2)
        row.append(np.mean(e))
    gram_tab.append((name, row)); print(f"    {name:24s}   " + "  ".join(f"{x:8.3f}" for x in row))
print(f"    3중 대각 행렬의 조건수: {np.linalg.cond(lap):.2e}  (N² 에 비례)")

# F3. 동역학 리 대수(DLA) 차원: HVA 생성원 {ΣZZ, ΣX} 대 ADAPT 풀 {Y_iZ_j}
def dla_dim(gens, tol=1e-7, cap=600):
    """생성원 iP 들의 교환자 닫힘의 실수 차원. 복소 행렬을 실벡터로 펴서 정규화한 뒤 Gram–Schmidt(두 번)로 독립성을 본다."""
    basis, mats, queue = [], [], [1j * g for g in gens]
    d = gens[0].shape[0]
    while queue and len(basis) < cap:
        M = queue.pop(0); v = np.concatenate([M.real.ravel(), M.imag.ravel()]); nv0 = np.linalg.norm(v)
        if nv0 < 1e-12: continue
        v = v / nv0
        for _ in range(2):
            for b in basis: v = v - (b @ v) * b
        if np.linalg.norm(v) > tol:
            v = v / np.linalg.norm(v); basis.append(v)
            Mn = v[:d * d].reshape(d, d) + 1j * v[d * d:].reshape(d, d)
            for Bm in mats: queue.append(Mn @ Bm - Bm @ Mn)
            mats.append(Mn)
    return len(basis)
print("F3) 동역학 리 대수 차원")
for n_ in (3, 4, 5):
    zz_ = sum(pauli(n_, {i: "Z", i + 1: "Z"}) for i in range(n_ - 1)).toarray(); xs_ = sum(pauli(n_, {i: "X"}) for i in range(n_)).toarray()
    msg = f"    n = {n_}: HVA 생성원 {{ΣZZ, ΣX}} → {dla_dim([zz_, xs_])} (= n²)"
    if n_ < 5: msg += f",  ADAPT 풀 {{Y_iZ_j}} → {dla_dim([P.toarray() for _, P, _ in make_pool(n_, 'YZ')])}"
    print(msg)
deg_pool = sorted({(abs(i - j), majorana_degree({i: "Y", j: "Z"}, 4, "X")) for i in range(4) for j in range(4) if i != j})
print(f"    풀 원소 Y_iZ_j 의 Majorana 차수 (|i−j|, 차수): {deg_pool}")

# F4. 17장의 급랭(|0…0⟩ 시작)은 이 Jordan–Wigner 에서 가우스 상태 하나가 아니다 — 패리티 성분 GHZ± 는 각각 가우스
print("F4) Wick 검사, 임계 TFIM 급랭 n = 12, t = 2: 시작 상태별 실제 S 와 가우스 S")
n = 12
zero12 = np.zeros(2 ** n, dtype=complex); zero12[0] = 1; one12 = np.zeros(2 ** n, dtype=complex); one12[-1] = 1
for name, psi0 in (("|0…0⟩", zero12), ("GHZ₊ = (|0…0⟩+|1…1⟩)/√2", (zero12 + one12) / np.sqrt(2)), ("GHZ₋ = (|0…0⟩−|1…1⟩)/√2", (zero12 - one12) / np.sqrt(2))):
    psi = expm_multiply(-2.0j * H12, psi0); Gm = np.zeros((n, n))
    for a in range(n):
        va = maj[a] @ psi
        for b in range(a + 1, n):
            Gm[a, b] = np.real(1j * np.vdot(va, maj[b] @ psi)); Gm[b, a] = -Gm[a, b]
    c1 = max(abs(np.vdot(psi, maj[a] @ psi)) for a in range(n))
    s_sv, s_g = half_entropy_sv(psi, n), gaussian_entropy(Gm, n // 2)
    print(f"    {name:28s}  S_실제 = {s_sv:.4f}   S_가우스 = {s_g:.4f}   차이 = {abs(s_g - s_sv):.1e}   max|⟨c_a⟩| = {c1:.3f}")

# F5. "시작점을 늘리면" 의 근거: (a) CCZ 를 6개로 만드는 구조별 시작점 하나의 성공률, (b) HVA 의 최선 오차
print("F5) 시작점 수 의존성")
six = sorted(set(both) | {canon(((0, 1), (0, 2), (1, 2), (1, 2), (0, 1), (0, 2)), list(itertools.permutations(range(3))))})   # 마지막 것은 두 실행이 모두 놓친 구조
print(f"    (a) 전결합 CNOT 6개 구조 {len(six)}개 각각에 시작점 하나짜리 시도 40번: 닿은 비율")
rates6 = []
for s_ in six:
    tm = Template(3, list(s_)); ok = sum(instantiate(tm, CCZ, rngF, starts=1)[0] < 1e-10 for _ in range(40)); rates6.append(ok / 40)
    print(f"        {' '.join(f'{a}{b}' for a, b in s_)}:  {ok}/40 = {ok/40:.2f}")
print(f"        한 번이라도 닿은 구조 {sum(r > 0 for r in rates6)}개, 성공률 범위 {min(rates6):.2f}~{max(rates6):.2f}")
Hd6 = tfim(6).toarray(); E06 = np.linalg.eigvalsh(Hd6)[0]; plus6 = np.ones(64, dtype=complex) / 8
hva_more = {}
for p_ in (6, 8):
    t5 = time.time(); Eb, _, fin_, _ = hva(Hd6, 6, p_, plus6, rngF, starts=250); hva_more[p_] = Eb - E06
    old = dict((pp // 2, e) for pp, e, _ in hva_pts[6])[p_]
    print(f"    (b) HVA n = 6, p = {p_}: 무작위 시작점 250개 중 최선 오차 {Eb-E06:.2e} (B3 의 21개 중 최선 {old:.2e}, {old/max(Eb-E06,1e-16):.0f}배)  ({time.time()-t5:.0f} s)", flush=True)
tF = time.time() - tF

# ---------- 그림 ----------
fig, axs = plt.subplots(3, 2, figsize=(12, 12.9))
FLOOR = 1e-16
ax = axs[0, 0]
cols = ["#4C6EF5", "#74C0FC", "#F59F00", "#F03E3E", "#868E96", "#2F9E44"]
for (name, _), c_ in zip(targets2, cols):
    ax.plot(range(4), [max(f, FLOOR) for f in tab2[name]], "o-", color=c_, label=name, lw=1.6, ms=5)
ax.plot([2] * 20, np.maximum(fl2, FLOOR), "x", color="#343A40", ms=5, alpha=0.6, label="Haar 무작위 20개 (k = 2, 3)")
ax.plot([3] * 20, np.maximum(fl3, FLOOR), "x", color="#343A40", ms=5, alpha=0.6)
ax.set_yscale("log"); ax.set_ylim(3e-17, 3); ax.set_xticks(range(4))
ax.set_xlabel("틀에 넣은 CNOT 수 k"); ax.set_ylabel("닿은 최소 비충실도 (로그)")
ax.set_title("(1) 2큐비트 합성: 필요한 CNOT 수는 0~3, 무작위 표적은 3"); ax.legend(fontsize=7.5, ncol=2, loc="lower left")

ax = axs[0, 1]
ax.step([r[0] for r in res_all], [max(r[2], FLOOR) for r in res_all], where="mid", color="#4C6EF5", lw=2, label="전결합 (쌍 3종)")
ax.plot([r[0] for r in res_all], [max(r[2], FLOOR) for r in res_all], "o", color="#4C6EF5")
ax.step([r[0] for r in res_line], [max(r[2], FLOOR) for r in res_line], where="mid", color="#F03E3E", lw=2, ls="--", label="선형 연결 0–1–2 (쌍 2종)")
ax.plot([r[0] for r in res_line], [max(r[2], FLOOR) for r in res_line], "s", color="#F03E3E")
ax.set_yscale("log"); ax.set_ylim(3e-17, 3); ax.set_xticks(range(9))
ax.annotate("6개에서 처음 닿음", (6, FLOOR), (3.2, 1e-12), fontsize=8, color="#4C6EF5", arrowprops=dict(arrowstyle="->", color="#4C6EF5"))
ax.annotate("선형은 8개", (8, FLOOR), (6.3, 1e-7), fontsize=8, color="#F03E3E", arrowprops=dict(arrowstyle="->", color="#F03E3E"))
ax.set_xlabel("CNOT 수 k (그 수의 모든 구조 중 최선)"); ax.set_ylabel("닿은 최소 비충실도 (로그)")
ax.set_title("(2) CCZ(Toffoli) 구조 전수 탐색: 연결이 줄면 6 → 8"); ax.legend(fontsize=8, loc="lower left")

ax = axs[1, 0]
for n_, c_ in zip((4, 6, 8), ("#74C0FC", "#4C6EF5", "#343A40")):
    h_ = adapt_curves[n_]
    ax.plot([x[0] for x in h_], [max(x[1], FLOOR) for x in h_], "-", color=c_, lw=1.8, label=f"ADAPT n = {n_} (풀 Y_iZ_j)")
    if n_ in hva_pts:
        ax.plot([x[0] for x in hva_pts[n_]], [max(x[1], 3e-11) for x in hva_pts[n_]], "^", color=c_, ms=7, mec="k", mew=0.5)
ax.plot([], [], "^", color="w", mec="k", label="HVA (같은 색 = 같은 n)")
trap_style = {("|+…+⟩", "Y"): ("-.", "함정 n = 6: 시작 |+…+⟩, 풀 {Y_i} (기울기 0)"),
              ("|0…0⟩", "Y"): ("--", "함정 n = 6: 시작 |0…0⟩, 풀 {Y_i} (평균장)"),
              ("|0…0⟩", "YZ"): (":", "함정 n = 6: 시작 |0…0⟩, 풀 {Y_iZ_j} (패리티)")}
for row in trap_rows:
    if row[0] == 6 and (row[1], row[2]) in trap_style:
        ls_, lab_ = trap_style[(row[1], row[2])]
        ax.axhline(abs(row[5]), color="#F03E3E", ls=ls_, lw=1.3, label=lab_)
ax.set_yscale("log"); ax.set_ylim(1e-11, 5); ax.set_xlabel("ADAPT: 게이트 수 / HVA: 파라미터 수 (10⁻¹¹ 아래는 바닥에 표시)"); ax.set_ylabel("에너지 오차 E − E₀ (로그)")
ax.set_title("(3) 변분 탐색: ADAPT는 구조 없이 찾지만 게이트 수가 빨리 는다"); ax.legend(fontsize=7, loc="lower right")

ax = axs[1, 1]
for m_, c_, mk in zip((250, 1000, 4000), ("#74C0FC", "#4C6EF5", "#343A40"), ("o", "s", "D")):
    ax.plot(r_list, fkv_curves[m_], mk + "-", color=c_, lw=1.5, ms=5, label=f"{m_}×{m_}")
ax.plot(r_list, fkv_curves[1000][0] * r_list[0] / np.array(r_list), "k:", lw=1, label="∝ 1/r")
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xticks(r_list); ax.set_xticklabels([str(r) for r in r_list])
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_xlabel("표본 행·열 수 r (= c)"); ax.set_ylabel("최적 계수-5 근사 대비 초과 오차 (로그)")
ax.set_title("(4) dequantization: 표본 수만 중요하고 행렬 크기는 무관"); ax.legend(fontsize=8)

ax = axs[2, 0]
xs_ = np.arange(len(srank_rows))
ax.bar(xs_, [r[2] for r in srank_rows], color=["#4C6EF5", "#74C0FC", "#F59F00", "#F03E3E"], width=0.6)
ax.axhline(1.0, color="#868E96", ls=":", lw=1); ax.text(-0.4, 1.15, "오차 1 = 추정이 쓸모없음", fontsize=8, color="#868E96")
for x_, r_ in zip(xs_, srank_rows): ax.text(x_, r_[2] * 1.12, f"안정 계수 {r_[1]:.0f}", ha="center", fontsize=8)
ax.set_xticks(xs_); ax.set_xticklabels(["계수 5\n(추천 행렬형)", "계수 40\n평평", "계수 300\n평평", "3중 대각\n(희소·최대 계수)"], fontsize=8)
ax.set_yscale("log"); ax.set_ylim(0.05, 20); ax.set_ylabel("Gram 오차 ‖AᵀA − SᵀS‖₂/‖A‖² (r = 200)")
ax.set_title("(5) dequantization이 깨지는 곳: 안정 계수가 크면 표본이 못 따라간다")

ax = axs[2, 1]
ax.plot([r[0] for r in ff_big], [r[1] for r in ff_big], "o-", color="#4C6EF5", lw=1.8, label=f"자유 페르미온 계산, n = 1000 (성장률 {rate:.2f} 비트/시간)")
for t_, S_ in ff_big[2:]:
    ax.annotate(f"MPS라면 χ ≥ 2^S ≈ {2.0**S_:.0e}", (t_, S_), xytext=(8, -16) if t_ < 40 else (-10, 8), textcoords="offset points",
                fontsize=8, color="#F03E3E", ha="left" if t_ < 40 else "right")
ax.set_ylim(0, 44)
ax.set_xlabel("시간 t"); ax.set_ylabel("반쪽 사슬 얽힘 엔트로피 S (비트)")
ax.set_title("(6) TFIM 급랭: MPS에는 지수 비용, 2000×2000 행렬로는 몇 초"); ax.legend(fontsize=8, loc="upper left")
fig.tight_layout(); fig.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)), "ch18_results.png"), dpi=120)

print(f"\n   (실행 시간: 합성 {tA:.0f} s + 변분 탐색 {tB:.0f} s + dequantization {tC:.0f} s + 자유 페르미온 {tD:.0f} s + 보강 {tF:.0f} s = 총 {time.time()-T0:.0f} s)")
