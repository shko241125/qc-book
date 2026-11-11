"""15장 양자 가속의 원천 — numpy 검증 스크립트
(1) 간섭: 양자 워크 vs 고전 랜덤 워크, 결어긋남으로 붕괴
(2) 얽힘만으로는 부족: Gottesman–Knill 안정자 테이블(CHP) 시뮬레이터
(3) 마법(T 게이트): 안정자 분기 합 비용 2^t
(4) 얽힘과 텐서망: 무작위 회로 vs 무작위 Clifford vs TFIM 바닥상태, MPS 절단 충실도
(5) 푸리에 샘플링의 최대 분리: Forrelation
규약: big-endian (q0 = 최상위 비트), seed=7
"""
import time
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
matplotlib.rcParams["font.family"] = ["Noto Sans CJK JP", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

rng = np.random.default_rng(7)
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]])
Z = np.diag([1.0, -1.0]).astype(complex)
H1 = np.array([[1, 1], [1, -1]], complex) / np.sqrt(2)
S1 = np.diag([1, 1j])
T1 = np.diag([1, np.exp(1j * np.pi / 4)])


def apply_1q(psi, g, q, n):
    psi = psi.reshape([2] * n)
    psi = np.moveaxis(np.tensordot(g, psi, axes=([1], [q])), 0, q)
    return psi.reshape(-1)


def apply_2q(psi, g4, a, b, n):
    psi = psi.reshape([2] * n)
    g = g4.reshape(2, 2, 2, 2)
    psi = np.tensordot(g, psi, axes=([2, 3], [a, b]))
    psi = np.moveaxis(psi, [0, 1], [a, b])
    return psi.reshape(-1)


CNOT4 = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]], complex)
SWAP4 = np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]], complex)


# ======================================================================
# (1) 간섭: 이산 시간 양자 워크(Hadamard 동전) vs 고전 랜덤 워크
# ======================================================================
def quantum_walk(steps, p_dephase=0.0, trials=1):
    """위치 -steps..steps, 동전 2차원. 동전 결어긋남은 매 스텝 확률 p 로 동전을 Z 기저에서 측정(혼합)하는 것으로
    모형화 — 밀도행렬 대신 궤적(trajectory) 평균."""
    L = 2 * steps + 1
    probs = np.zeros(L)
    for _ in range(trials):
        psi = np.zeros((L, 2), complex)
        psi[steps, 0] = 1 / np.sqrt(2); psi[steps, 1] = 1j / np.sqrt(2)   # 대칭 초기 동전
        for _ in range(steps):
            psi = psi @ H1.T                                  # 동전 던지기
            if p_dephase > 0 and rng.random() < p_dephase:    # 동전 측정 → 간섭 파괴
                pc = (np.abs(psi) ** 2).sum(axis=0)
                c = rng.choice(2, p=pc / pc.sum())
                keep = np.zeros_like(psi); keep[:, c] = psi[:, c]
                psi = keep / np.linalg.norm(keep)
            new = np.zeros_like(psi)
            new[:-1, 0] = psi[1:, 0]                          # 동전 0 → 왼쪽
            new[1:, 1] = psi[:-1, 1]                          # 동전 1 → 오른쪽
            psi = new
        probs += (np.abs(psi) ** 2).sum(axis=1)
    return probs / trials


def classical_walk(steps):
    from math import comb
    L = 2 * steps + 1
    p = np.zeros(L)
    for k in range(steps + 1):
        p[steps + 2 * k - steps] = comb(steps, k) / 2 ** steps
    return p


print("1) 간섭: 양자 워크 vs 고전 랜덤 워크")
walk = {}
for steps in (25, 50, 100, 200):
    xs = np.arange(-steps, steps + 1)
    pq = quantum_walk(steps); pc = classical_walk(steps)
    sq = np.sqrt((pq * xs ** 2).sum()); sc = np.sqrt((pc * xs ** 2).sum())
    walk[steps] = (sq, sc)
    print(f"   t={steps:3d}: σ_양자 = {sq:6.2f} (σ/t = {sq/steps:.3f}),  σ_고전 = {sc:6.2f} (σ/√t = {sc/np.sqrt(steps):.3f})")
steps = 100; xs = np.arange(-steps, steps + 1)
deph = {}
for p in (0.0, 0.05, 0.2, 1.0):
    pq = quantum_walk(steps, p, trials=1 if p == 0 else 2000)
    deph[p] = pq
    print(f"   t=100, 동전 결어긋남 확률 p={p:4.2f}: σ = {np.sqrt((pq * xs**2).sum()):6.2f}")


# ======================================================================
# (2) Gottesman–Knill: CHP 안정자 테이블 (Aaronson–Gottesman 2004)
#     행 0..n-1 = 파괴자(destabilizer), n..2n-1 = 안정자, 2n = 작업 행. 비트 x, z 와 부호 r.
# ======================================================================
class Tableau:
    def __init__(self, n):
        self.n = n
        self.x = np.zeros((2 * n + 1, n), dtype=np.uint8)
        self.z = np.zeros((2 * n + 1, n), dtype=np.uint8)
        self.r = np.zeros(2 * n + 1, dtype=np.uint8)
        for i in range(n):
            self.x[i, i] = 1          # 파괴자 X_i
            self.z[n + i, i] = 1      # 안정자 Z_i  (|0...0>)

    def h(self, a):
        x, z = self.x, self.z
        self.r ^= x[:, a] & z[:, a]
        x[:, a], z[:, a] = z[:, a].copy(), x[:, a].copy()

    def s(self, a):
        x, z = self.x, self.z
        self.r ^= x[:, a] & z[:, a]
        z[:, a] ^= x[:, a]

    def cnot(self, a, b):
        x, z = self.x, self.z
        self.r ^= x[:, a] & z[:, b] & (x[:, b] ^ z[:, a] ^ 1)
        x[:, b] ^= x[:, a]
        z[:, a] ^= z[:, b]

    @staticmethod
    def _g(x1, z1, x2, z2):
        """i 의 지수: 파울리 곱 (x1,z1)(x2,z2) 에서 생기는 위상 (원소별, 값 −1/0/1)."""
        x1 = x1.astype(np.int8); z1 = z1.astype(np.int8); x2 = x2.astype(np.int8); z2 = z2.astype(np.int8)
        return np.where((x1 == 0) & (z1 == 0), 0,
               np.where((x1 == 1) & (z1 == 1), z2 - x2,
               np.where((x1 == 1) & (z1 == 0), z2 * (2 * x2 - 1), x2 * (1 - 2 * z2))))

    def rowsum(self, hs, i):
        """행 집합 hs 각각에 행 i 를 곱한다 (벡터화)."""
        if len(hs) == 0: return
        hs = np.asarray(hs)
        tot = 2 * self.r[hs].astype(np.int64) + 2 * int(self.r[i]) + \
              self._g(self.x[i][None, :], self.z[i][None, :], self.x[hs], self.z[hs]).sum(axis=1)
        self.r[hs] = ((tot % 4) == 2).astype(np.uint8)
        self.x[hs] ^= self.x[i]; self.z[hs] ^= self.z[i]

    def measure(self, a, rng):
        n = self.n
        p_rows = np.nonzero(self.x[n:2 * n, a])[0]
        if len(p_rows):                                    # 무작위 결과
            p = n + p_rows[0]
            others = [i for i in range(2 * n) if i != p and self.x[i, a]]
            self.rowsum(others, p)
            self.x[p - n] = self.x[p]; self.z[p - n] = self.z[p]; self.r[p - n] = self.r[p]
            self.x[p] = 0; self.z[p] = 0; self.z[p, a] = 1
            self.r[p] = rng.integers(2)
            return int(self.r[p]), False
        w = 2 * n                                          # 결정적 결과: 작업 행에 누적
        self.x[w] = 0; self.z[w] = 0; self.r[w] = 0
        for i in np.nonzero(self.x[:n, a])[0]:
            self.rowsum([w], i + n)
        return int(self.r[w]), True

    def stabilizers(self):
        n = self.n
        return [(self.x[n + i].copy(), self.z[n + i].copy(), int(self.r[n + i])) for i in range(n)]


def pauli_from_bits(xb, zb, r):
    m = {(0, 0): I2, (1, 0): X, (0, 1): Z, (1, 1): Y}
    out = np.array([[1.0 + 0j]])
    for xi, zi in zip(xb, zb):
        out = np.kron(out, m[(int(xi), int(zi))])
    return (-1) ** r * out


def random_clifford_ops(n, ngates, rng):
    ops = []
    for _ in range(ngates):
        k = rng.integers(3)
        if k == 0: ops.append(("h", rng.integers(n)))
        elif k == 1: ops.append(("s", rng.integers(n)))
        else:
            a, b = rng.choice(n, 2, replace=False); ops.append(("cx", a, b))
    return ops


print("\n2) Gottesman–Knill: 안정자 테이블 시뮬레이터")
# 2a) 정확성: n=8 무작위 Clifford 회로 → 상태벡터가 테이블의 안정자 8개 모두의 +1 고유상태인가
ok_all = True
for trial in range(20):
    n = 8; ops = random_clifford_ops(n, 120, rng)
    tab = Tableau(n); psi = np.zeros(2 ** n, complex); psi[0] = 1
    for op in ops:
        if op[0] == "h": tab.h(op[1]); psi = apply_1q(psi, H1, op[1], n)
        elif op[0] == "s": tab.s(op[1]); psi = apply_1q(psi, S1, op[1], n)
        else: tab.cnot(op[1], op[2]); psi = apply_2q(psi, CNOT4, op[1], op[2], n)
    for xb, zb, r in tab.stabilizers():
        P = pauli_from_bits(xb, zb, r)
        if not np.allclose(P @ psi, psi, atol=1e-9): ok_all = False
print(f"   n=8, 무작위 Clifford 120 게이트 × 20 회로: 상태벡터가 테이블 안정자 8개의 +1 고유상태 — {ok_all}")

# 2b) 측정 분포 검증: n=5 회로, 테이블 샘플링 vs 상태벡터 |amp|^2
n = 5; ops = random_clifford_ops(n, 40, rng)
psi = np.zeros(2 ** n, complex); psi[0] = 1
for op in ops:
    if op[0] == "h": psi = apply_1q(psi, H1, op[1], n)
    elif op[0] == "s": psi = apply_1q(psi, S1, op[1], n)
    else: psi = apply_2q(psi, CNOT4, op[1], op[2], n)
p_sv = np.abs(psi) ** 2
counts = np.zeros(2 ** n)
for shot in range(4000):
    tab = Tableau(n)
    for op in ops:
        getattr(tab, {"h": "h", "s": "s", "cx": "cnot"}[op[0]])(*op[1:])
    bits = [tab.measure(q, rng)[0] for q in range(n)]
    counts[int("".join(map(str, bits)), 2)] += 1
tvd_tab = 0.5 * np.abs(counts / 4000 - p_sv).sum()
supp = int((p_sv > 1e-12).sum())
print(f"   n=5 측정 4000샷: 테이블 vs 상태벡터 TVD = {tvd_tab:.3f} (지지 {supp}개 균등 = 2^k; 샷 잡음 수준)")

# 2c) GHZ: 최대 얽힘이지만 테이블로 즉시
n = 1000; tab = Tableau(n); tab.h(0)
for i in range(n - 1): tab.cnot(i, i + 1)
first, det0 = tab.measure(0, rng)
rest = [tab.measure(q, rng) for q in (1, 500, 999)]
print(f"   GHZ n=1000: q0 측정 = {first} (무작위={not det0}), q1·q500·q999 = {[b for b, _ in rest]} (결정적={all(d for _, d in rest)})")

# 2d) 확장성: 게이트 20n + 전 큐비트 측정, n = 50..800
clif_times = []
for n in (50, 100, 200, 400, 800):
    ops = random_clifford_ops(n, 20 * n, rng)
    t0 = time.perf_counter(); tab = Tableau(n)
    for op in ops:
        if op[0] == "h": tab.h(op[1])
        elif op[0] == "s": tab.s(op[1])
        else: tab.cnot(op[1], op[2])
    t1 = time.perf_counter()
    for q in range(n): tab.measure(q, rng)
    t2 = time.perf_counter()
    clif_times.append((n, t1 - t0, t2 - t1))
    print(f"   n={n:4d}: 게이트 {20*n:6d}개 {t1-t0:6.2f} s, 측정 {n}개 {t2-t1:6.2f} s  (상태벡터라면 2^{n} 진폭)")


# ======================================================================
# (3) 마법: T = e^{iπ/8}[cos(π/8) I − i sin(π/8) Z] — Clifford 두 개의 합.
#     T 가 t 개인 회로 = Clifford 회로 2^t 개의 가중합 → 비용 2^t (n 이 아니라 t 에 지수)
# ======================================================================
def random_clifford_t(n, ngates, t_count, rng):
    ops = random_clifford_ops(n, ngates, rng)
    pos = sorted(rng.choice(len(ops) + t_count, t_count, replace=False))
    out, k = [], 0
    for i in range(len(ops) + t_count):
        if k < t_count and i == pos[k]: out.append(("t", rng.integers(n))); k += 1
        else: out.append(ops.pop(0))
    return out


def run_sv(ops, n, branch=None):
    """branch=None: 정확한 T. branch=비트열: k번째 T 를 I(0) 또는 Z(1) 로 바꾼 Clifford 회로."""
    psi = np.zeros(2 ** n, complex); psi[0] = 1; k = 0
    for op in ops:
        if op[0] == "h": psi = apply_1q(psi, H1, op[1], n)
        elif op[0] == "s": psi = apply_1q(psi, S1, op[1], n)
        elif op[0] == "cx": psi = apply_2q(psi, CNOT4, op[1], op[2], n)
        else:
            if branch is None: psi = apply_1q(psi, T1, op[1], n)
            else:
                if branch[k]: psi = apply_1q(psi, Z, op[1], n)
                k += 1
    return psi


print("\n3) 마법: T 게이트 하나 = Clifford 두 개의 합 → 분기 2^t")
c0 = np.exp(1j * np.pi / 8) * np.cos(np.pi / 8); c1 = np.exp(1j * np.pi / 8) * (-1j) * np.sin(np.pi / 8)
assert np.allclose(c0 * I2 + c1 * Z, T1)
t_rows = []
n = 6
for t in (0, 2, 4, 6, 8, 10):
    ops = random_clifford_t(n, 60, t, rng)
    exact = run_sv(ops, n)
    t0 = time.perf_counter()
    approx = np.zeros(2 ** n, complex)
    for b in range(2 ** t):
        bb = [(b >> (t - 1 - j)) & 1 for j in range(t)]
        w = np.prod([c1 if x else c0 for x in bb]) if t else 1.0
        approx += w * run_sv(ops, n, bb)
    dt = time.perf_counter() - t0
    l1 = (abs(c0) + abs(c1)) ** t
    t_rows.append((t, 2 ** t, dt, l1 ** 2))
    print(f"   n={n}, T {t:2d}개: 분기 {2**t:5d}개, 합 = 정확 상태? {np.allclose(approx, exact)}, 시간 {dt:6.3f} s, "
          f"‖c‖₁² = {l1**2:7.2f} (= 1.707^t, 샘플링 비용)")
print(f"   T 상태의 안정자 범위 ξ(T) = 1/cos²(π/8) = {1/np.cos(np.pi/8)**2:.4f} → 최적 분해 비용 ∝ ξ^t = 2^({np.log2(1/np.cos(np.pi/8)**2):.3f}·t)")


# ======================================================================
# (4) 얽힘과 텐서망: 반쪽 사슬 얽힘 엔트로피 S (비트), MPS 절단 충실도
# ======================================================================
def haar_u4(rng):
    A = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
    Q, R = np.linalg.qr(A)
    return Q * (np.diag(R) / np.abs(np.diag(R)))


def random_clifford_2q(rng):
    """H, S, CNOT(양방향) 40~41개를 무작위로 곱해 2큐비트 Clifford 군(위상 무시 11520원소)을 거의 균등하게 뽑는다.
    생성원 그래프가 이분 그래프라 길이가 늘 짝수면 절반(5760원소)에만 닿는다 → 길이의 홀짝도 무작위로."""
    U = np.eye(4, dtype=complex)
    for _ in range(40 + int(rng.integers(2))):
        k = rng.integers(5)
        if k == 0: U = np.kron(H1, I2) @ U
        elif k == 1: U = np.kron(I2, H1) @ U
        elif k == 2: U = np.kron(S1, I2) @ U
        elif k == 3: U = np.kron(I2, S1) @ U
        else: U = (CNOT4 if rng.random() < 0.5 else SWAP4 @ CNOT4 @ SWAP4) @ U
    return U


def schmidt(psi, n, cut):
    s = np.linalg.svd(psi.reshape(2 ** cut, 2 ** (n - cut)), compute_uv=False)
    return s ** 2


def entropy_bits(p):
    p = p[p > 1e-14]
    return float(-(p * np.log2(p)).sum())


print("\n4) 얽힘: 반쪽 사슬 엔트로피 S (비트) vs 벽돌 회로 깊이, n = 14")
n = 14; cut = n // 2
depths = list(range(0, 31, 3))
S_haar, S_cliff = [], []
psi_h = np.zeros(2 ** n, complex); psi_h[0] = 1
psi_c = psi_h.copy()
snap = {}
for d in range(max(depths) + 1):
    if d in depths:
        S_haar.append(entropy_bits(schmidt(psi_h, n, cut))); S_cliff.append(entropy_bits(schmidt(psi_c, n, cut)))
        if d in (3, 9, 18, 30): snap[d] = psi_h.copy()
    start = d % 2
    for q in range(start, n - 1, 2):
        psi_h = apply_2q(psi_h, haar_u4(rng), q, q + 1, n)
        psi_c = apply_2q(psi_c, random_clifford_2q(rng), q, q + 1, n)
page = n / 2 - 1 / (2 * np.log(2))
for d, sh, sc in zip(depths, S_haar, S_cliff):
    print(f"   깊이 {d:2d}: S_무작위 = {sh:5.2f}, S_Clifford = {sc:5.2f}")
print(f"   Page 값(무작위 상태의 평균) ≈ n/2 − 1/(2 ln 2) = {page:.2f} 비트, 최대 = {n//2}")
# 표본 하나는 요동이 크므로 Haar vs Clifford 성장 속도를 10회 평균으로 비교 (별도 시드 — 위의 rng 흐름은 그대로)
avg = {6: [[], []], 12: [[], []]}
for k in range(10):
    r2 = np.random.default_rng(1000 + k)
    ph_, pc_ = np.zeros(2 ** n, complex), np.zeros(2 ** n, complex); ph_[0] = pc_[0] = 1
    for d in range(12):
        for q in range(d % 2, n - 1, 2):
            ph_ = apply_2q(ph_, haar_u4(r2), q, q + 1, n); pc_ = apply_2q(pc_, random_clifford_2q(r2), q, q + 1, n)
        if d + 1 in avg:
            avg[d + 1][0].append(entropy_bits(schmidt(ph_, n, cut))); avg[d + 1][1].append(entropy_bits(schmidt(pc_, n, cut)))
for d in (6, 12):
    h_, c_ = np.array(avg[d][0]), np.array(avg[d][1])
    print(f"   10회 평균, 깊이 {d:2d}: Haar {h_.mean():.2f} ± {h_.std(ddof=1)/np.sqrt(10):.2f}, Clifford {c_.mean():.2f} ± {c_.std(ddof=1)/np.sqrt(10):.2f} 비트")

# TFIM 바닥상태 (12·13장의 모델, n=14) — 면적 법칙
from scipy.sparse import kron as skron, identity as sid, csr_matrix
from scipy.sparse.linalg import eigsh
def op_at(o, i, n):
    mats = [sid(2, format="csr") for _ in range(n)]; mats[i] = csr_matrix(o)
    out = mats[0]
    for m in mats[1:]: out = skron(out, m, format="csr")
    return out
Htf = sum(-op_at(Z.real, i, n) @ op_at(Z.real, i + 1, n) for i in range(n - 1)) + sum(-op_at(X.real, i, n) for i in range(n))
e0, v0 = eigsh(Htf, k=1, which="SA")
S_tfim = entropy_bits(schmidt(v0[:, 0].astype(complex), n, cut))
print(f"   TFIM 바닥상태 n=14 (임계점 J=h=1): S = {S_tfim:.3f} 비트 — 면적 법칙(임계라 log n 으로만 증가)")

print("   MPS 절단 충실도 (가운데 자름, 슈미트 상위 χ 개 유지):")
mps_rows = {}
for d, psi in sorted(snap.items()):
    p = np.sort(schmidt(psi, n, cut))[::-1]
    row = [p[:chi].sum() for chi in (2, 8, 32, 128)]
    mps_rows[d] = row
    print(f"   깊이 {d:2d}: χ=2 {row[0]:.3f}, χ=8 {row[1]:.3f}, χ=32 {row[2]:.3f}, χ=128 {row[3]:.4f}")
p = np.sort(schmidt(v0[:, 0].astype(complex), n, cut))[::-1]
row_tfim = [p[:chi].sum() for chi in (2, 8, 32, 128)]
print(f"   TFIM 바닥상태: χ=2 {row_tfim[0]:.4f}, χ=8 {row_tfim[1]:.6f}, χ=32 {row_tfim[2]:.8f}")


# ======================================================================
# (5) 푸리에 샘플링의 최대 분리: Forrelation (Aaronson 2010, Aaronson–Ambainis 2015)
#     Φ(f,g) = 2^{-3n/2} Σ_{x,y} f(x) (−1)^{x·y} g(y),  양자 2쿼리 회로의 P(0^n) = Φ²
# ======================================================================
def wht(v):
    """정규화된 Walsh–Hadamard 변환 (H^{⊗n} 을 벡터에) — O(N log N)."""
    v = v.astype(float).copy(); h = 1
    while h < len(v):
        for i in range(0, len(v), 2 * h):
            a = v[i:i + h].copy(); b = v[i + h:i + 2 * h].copy()
            v[i:i + h] = a + b; v[i + h:i + 2 * h] = a - b
        h *= 2
    return v / np.sqrt(len(v))


def forrelation(f, g):
    return float(np.dot(wht(f), g)) / len(f)


def quantum_p0(f, g):
    """|0^n> → H^n → U_f(위상) → H^n → U_g(위상) → H^n, 출력 |0^n> 확률."""
    N = len(f); psi = np.full(N, 1 / np.sqrt(N))
    psi = wht(f * psi); psi = wht(g * psi)
    return psi[0] ** 2


print("\n5) Forrelation: 1~2 쿼리 양자 vs Ω̃(√N) 고전")
for n in (6, 10, 14):
    N = 2 ** n
    phis_r, phis_f, p_r, p_f = [], [], [], []
    for trial in range(200):
        f = rng.choice([-1, 1], N)
        g_rand = rng.choice([-1, 1], N)
        g_forr = np.sign(wht(f)); g_forr[g_forr == 0] = 1
        phis_r.append(forrelation(f, g_rand)); phis_f.append(forrelation(f, g_forr))
        if trial < 20:
            p_r.append(quantum_p0(f, g_rand)); p_f.append(quantum_p0(f, g_forr))
    phis_r, phis_f = np.array(phis_r), np.array(phis_f)
    if n == 10: forr_hist = (phis_r, phis_f)
    print(f"   n={n:2d} (N={N:5d}): 무작위 쌍 |Φ| 평균 {np.mean(np.abs(phis_r)):.4f} (≈ √(2/π)·N^-1/2 = {np.sqrt(2/np.pi)*N**-0.5:.4f}), "
          f"forrelated Φ 평균 {phis_f.mean():.3f} (√(2/π) = {np.sqrt(2/np.pi):.3f}); "
          f"양자 P(0^n): 무작위 {np.mean(p_r):.1e}, forrelated {np.mean(p_f):.3f}; 고전 쿼리 하한 ~√N/log N = {np.sqrt(N)/np.log2(N):.0f}")
print("   양자: 쿼리 2번(제어형이면 1번)의 회로 몇 번 반복으로 구분; 고전 무작위 알고리즘은 Ω(√N / log N) 쿼리 필요(Aaronson–Ambainis 2015)")


# ---------- 그림 (2×3) ----------
ns = np.array([r[0] for r in clif_times]); ttot = np.array([r[1] + r[2] for r in clif_times])
slope = np.polyfit(np.log(ns[1:]), np.log(ttot[1:]), 1)[0]
print(f"   (그림용) Clifford 총 시간의 log-log 기울기 (n=100..800) = {slope:.2f}")
fig, axs = plt.subplots(2, 3, figsize=(15, 8))
ax = axs[0, 0]; steps = 100; xs = np.arange(-steps, steps + 1); m = (xs % 2 == 0)
for p, c in ((0.0, "#4C6EF5"), (0.05, "#7950F2"), (0.2, "#F59F00"), (1.0, "#868E96")):
    ax.plot(xs[m], deph[p][m], color=c, lw=1.2, label=f"양자 워크, 결어긋남 p={p}")
ax.plot(xs[m], classical_walk(steps)[m], "k:", lw=1.2, label="고전 랜덤 워크")
ax.set_xlabel("위치 x (t = 100, 짝수 자리)"); ax.set_ylabel("확률"); ax.set_title("(1) 간섭: σ ∝ t (양자) vs √t (고전)"); ax.legend(fontsize=7)

ax = axs[0, 1]
ax.loglog(ns, ttot, "o-", color="#4C6EF5", label="안정자 테이블 총 시간")
ax.loglog(ns, ttot[-1] * (ns / ns[-1]) ** slope, "k:", label=f"∝ n^{slope:.2f} (적합)")
ax.set_xticks(ns); ax.set_xticklabels([str(v) for v in ns]); ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_xlabel("큐비트 수 n"); ax.set_ylabel("시간 (s)"); ax.set_title("(2) Clifford 회로: n 에 다항"); ax.legend(fontsize=8)

ax = axs[0, 2]
ax.semilogy([r[0] for r in t_rows], [r[1] for r in t_rows], "s-", color="#F03E3E", label="정확한 분기 수 2^t")
ax.semilogy([r[0] for r in t_rows], [r[3] for r in t_rows], "^-", color="#F59F00", label="이 분해의 샘플링 비용 1.707^t")
tt = np.array([r[0] for r in t_rows]); ax.semilogy(tt, (1 / np.cos(np.pi / 8) ** 2) ** tt, "d-", color="#2F9E44", label=r"최적 범위 $\xi^t = 2^{0.228t}$")
ax.set_xlabel("T 게이트 수 t"); ax.set_ylabel("비용"); ax.set_title("(3) Clifford+T: t 에 지수"); ax.legend(fontsize=8)

ax = axs[1, 0]
ax.plot(depths, S_haar, "o-", color="#4C6EF5", label="무작위 2큐비트 게이트 (Haar)")
ax.plot(depths, S_cliff, "s-", color="#F03E3E", label="무작위 Clifford (고전으로 쉬움)")
ax.axhline(page, ls="--", color="#868E96", lw=0.8, label="Page 값")
ax.axhline(S_tfim, ls=":", color="#2F9E44", lw=1.2, label="TFIM 바닥상태 (임계점: log n 성장)")
ax.set_xlabel("벽돌 회로 깊이"); ax.set_ylabel("반쪽 사슬 엔트로피 (비트)"); ax.set_title("(4a) 얽힘 성장, n = 14"); ax.legend(fontsize=7)

ax = axs[1, 1]; chis = [2, 8, 32, 128]
for d, c in zip(sorted(mps_rows), ("#74C0FC", "#4C6EF5", "#364FC7", "#1C2D7A")):
    ax.semilogx(chis, mps_rows[d], "o-", color=c, label=f"무작위 회로 깊이 {d}")
ax.semilogx(chis, row_tfim, "s--", color="#2F9E44", label="TFIM 바닥상태")
ax.set_xticks(chis); ax.set_xticklabels([str(c) for c in chis])
ax.set_xlabel("결합 차원 χ (가운데 자름)"); ax.set_ylabel("절단 충실도"); ax.set_title("(4b) MPS 절단: 얽힘이 비용"); ax.legend(fontsize=7)

ax = axs[1, 2]
phis_r, phis_f = forr_hist
ax.hist(phis_r, bins=30, color="#868E96", alpha=0.8, label="무작위 (f, g)")
ax.hist(phis_f, bins=30, color="#4C6EF5", alpha=0.8, label="forrelated (f, g)")
ax.set_xlabel("Φ(f, g)"); ax.set_ylabel("빈도 (200쌍)"); ax.set_title("(5) Forrelation, n = 10 (N = 1024)"); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)), "ch15_results.png"), dpi=120)
print("6) saved ch15_results.png")
