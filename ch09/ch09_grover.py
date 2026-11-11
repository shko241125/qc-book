"""9장 Grover 탐색과 진폭 증폭 — 오라클·확산 연산자를 게이트로 조립하고 최적 반복 횟수를 numpy 로 검증한다.
실행: python ch09_grover.py
환경: Python 3.13, numpy 2.5 (SDK 없이 동작). 큐비트 순서: 빅 엔디안(3장 규약).
"""
import numpy as np

np.set_printoptions(precision=4, suppress=True, linewidth=120)
rng = np.random.default_rng(seed=7)

I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)

def kron(*ops):
    out = np.array([[1]], dtype=complex) if ops[0].ndim == 2 else np.array([1], dtype=complex)
    for o in ops: out = np.kron(out, o)
    return out
def single(g, q, n): ops = [I]*n; ops[q] = g; return kron(*ops)
def Hn(n): return kron(*([H]*n))
def Xn(n): return kron(*([X]*n))
def ket0n(n): v = np.zeros(2**n, dtype=complex); v[0] = 1; return v
def probs(psi): return abs(psi)**2

# ── 1. 오라클과 확산 연산자 ─────────────────────────────────────────
def oracle(marked, n):                  # O = I - 2 Σ|w><w|  : 표시 원소에 -1
    d = np.ones(2**n, dtype=complex); d[list(marked)] = -1; return np.diag(d)
def mcz(n):                             # 다중 제어 Z: |1...1> 에만 -1
    d = np.ones(2**n, dtype=complex); d[-1] = -1; return np.diag(d)
def diffusion_gates(n):                 # D = H^n X^n (MCZ) X^n H^n  = 2|s><s| - I  (전역 위상 제외)
    return Hn(n) @ Xn(n) @ mcz(n) @ Xn(n) @ Hn(n)
def diffusion_formula(n):
    s = Hn(n) @ ket0n(n); return 2*np.outer(s, s.conj()) - np.eye(2**n)
n = 3
print("1) 확산 연산자: 게이트 분해 H X MCZ X H == 2|s><s|-I (부호까지)?", np.allclose(diffusion_gates(n), -diffusion_formula(n)),
      " / 전역 위상 -1 차이 → 둘 다 Grover 반복에 동일")
print("   n=3 의 다중제어 Z 는 CCZ = H·Toffoli·H (4장: Toffoli = 6 CNOT, 7 T)")

# ── 2. Grover 반복: 성공 확률 vs 반복 횟수 ────────────────────────────
def grover_probs(n, marked, kmax):
    N = 2**n; O = oracle(marked, n); D = diffusion_formula(n); G = D @ O
    psi = Hn(n) @ ket0n(n); out = []
    for k in range(kmax + 1):
        out.append(probs(psi)[list(marked)].sum()); psi = G @ psi
    return np.array(out)
def theory(n, M, k):
    N = 2**n; th = np.arcsin(np.sqrt(M/N)); return np.sin((2*k+1)*th)**2
def k_opt(n, M):
    N = 2**n; th = np.arcsin(np.sqrt(M/N)); return int(np.floor(np.pi/(4*th)))
print("\n2) 단일 해: 성공 확률 vs 반복 k  (이론 sin²((2k+1)θ), θ = arcsin(1/√N))")
rows = []
for n in (3, 4, 6, 8, 10):
    w = int(rng.integers(0, 2**n)); kmax = k_opt(n, 1) + 3
    p = grover_probs(n, [w], kmax); ko = k_opt(n, 1)
    ok = np.allclose(p, [theory(n, 1, k) for k in range(kmax + 1)])
    rows.append((n, ko, p[ko], p[ko+1], p[0]))
    print(f"   n={n:>2d} N={2**n:>5d}: 최적 k = ⌊π/4·√N⌋ = {ko:>2d}, P(k_opt) = {p[ko]:.4f}, P(k_opt+1) = {p[ko+1]:.4f} (지나침), P(0) = 1/N = {p[0]:.4f}, 이론 일치 {ok}")
# n=3 전체 곡선 (그림용)
p3 = grover_probs(3, [5], 8)
print("   n=3, 표시=|101>: k=0..8 →", np.round(p3, 3), " (k=2 에서 0.945, k=6 에서 1.000 = sin²(13θ)≈1 — 주기적으로 오르내린다)")

# ── 3. 다중 해 M 개 ───────────────────────────────────────────────
print("\n3) 다중 해: n=8, N=256")
for M in (1, 2, 4, 16, 64):
    marked = rng.choice(256, M, replace=False); ko = k_opt(8, M)
    p = grover_probs(8, marked, ko + 1)
    print(f"   M={M:>2d}: 최적 k = {ko:>2d}, P = {p[ko]:.4f}  (이론 {theory(8, M, ko):.4f}),  √(N/M) = {np.sqrt(256/M):.1f}")
# 연습문제 2: M = N/4 → θ = π/6, 1번에 sin²(3·π/6) = 1
p = grover_probs(4, list(range(4)), 1)
print(f"   연습문제 2: N=16, M=4 → θ=π/6, k=1 에서 P = {p[1]:.4f} (정확히 1)")
p = grover_probs(2, [3], 1)
print(f"   연습문제 1: N=4, M=1 → θ=π/6, k=1 에서 P = {p[1]:.4f} (정확히 1)")

# ── 4. 진폭 증폭: A 가 H^n 이 아닌 일반 상태 준비 회로일 때 ────────────────
def amp_amp(A, good, kmax):
    """Q = -A S_0 A† S_f.  A: 상태 준비 유니터리, good: '좋은' 기저 상태 인덱스"""
    N = A.shape[0]; S_f = oracle(good, int(np.log2(N)))
    S_0 = np.eye(N, dtype=complex); S_0[0, 0] = -1
    Q = -A @ S_0 @ A.conj().T @ S_f
    psi = A @ ket0n(int(np.log2(N))); out = []
    for k in range(kmax + 1):
        out.append(probs(psi)[good].sum()); psi = Q @ psi
    return np.array(out)
n = 5
# 난수 유니터리 A (Haar 근사) → 초기 성공 확률 a = sin²θ
Zr = rng.normal(size=(2**n, 2**n)) + 1j*rng.normal(size=(2**n, 2**n)); Q_, R_ = np.linalg.qr(Zr); A = Q_ @ np.diag(np.diag(R_)/abs(np.diag(R_)))
good = [7]
a = probs(A @ ket0n(n))[7]; th = np.arcsin(np.sqrt(a)); ko = int(np.floor(np.pi/(4*th)))
p = amp_amp(A, good, ko + 2)
print(f"\n4) 진폭 증폭 (난수 A, n={n}): 초기 P = {a:.4f} → θ = {th:.4f}, k_opt = {ko}, P(k_opt) = {p[ko]:.4f}, 이론 sin²((2k+1)θ) = {np.sin((2*ko+1)*th)**2:.4f}")
print("   → 가속은 A 가 무엇이든 1/√a 번: '검색'이 아니라 '성공 확률 증폭'이 본질")

# ── 5. 해 수를 모를 때: 무작위 반복 횟수 (연습문제 3) ──────────────────
def unknown_M_search(n, marked, trials=200):
    """BBHT: m=1 부터 λ=6/5 배씩 키우며 k ∈ [0, m) 균등 난수로 반복, 측정 결과 검증"""
    N = 2**n; total_q = []
    for _ in range(trials):
        m = 1.0; queries = 0
        while True:
            k = int(rng.integers(0, int(np.ceil(m))))
            p = grover_probs(n, marked, k)[k]; queries += k + 1   # 오라클 k회 + 검증 1회
            if rng.random() < p: total_q.append(queries); break
            m = min(m * 1.2, np.sqrt(N))
    return np.mean(total_q)
print("\n5) 해 수를 모를 때 (BBHT 무작위 반복, n=8):")
for M in (1, 4, 16):
    marked = rng.choice(256, M, replace=False)
    print(f"   M={M:>2d}: 평균 총 오라클 호출 {unknown_M_search(8, marked):.1f}  (알 때 k_opt+1 = {k_opt(8, M)+1}, √(N/M) = {np.sqrt(256/M):.1f})")

# ── 6. 노이즈: 반복마다 위상 감쇠 → 성공률 vs k ─────────────────────────
def phase_damp_all(rho, lam, N_q):
    K0 = np.array([[1, 0], [0, np.sqrt(1-lam)]]); K1 = np.array([[0, 0], [0, np.sqrt(lam)]])
    for q in range(N_q):
        A_ = single(K0, q, N_q); B_ = single(K1, q, N_q)
        rho = A_ @ rho @ A_.conj().T + B_ @ rho @ B_.conj().T
    return rho
def grover_noisy(n, marked, kmax, lam):
    O = oracle(marked, n); D = diffusion_formula(n); G = D @ O
    psi = Hn(n) @ ket0n(n); rho = np.outer(psi, psi.conj()); out = []
    for k in range(kmax + 1):
        out.append(np.real(np.diag(rho))[list(marked)].sum())
        rho = G @ rho @ G.conj().T; rho = phase_damp_all(rho, lam, n)
    return np.array(out)
n = 6; w = [37]; kmax = 12
print(f"\n6) 노이즈: n={n}, 반복마다 큐비트당 위상 감쇠 λ → 성공 확률 (k_opt = {k_opt(6,1)})")
noisy = {}
for lam in (0.0, 0.01, 0.03, 0.1):
    p = grover_noisy(n, w, kmax, lam); noisy[lam] = p
    print(f"   λ={lam:.2f}: max P = {p.max():.4f} at k = {p.argmax()},  P(k_opt) = {p[k_opt(6,1)]:.4f}")

# ── 7. 그림 ───────────────────────────────────────────────────────
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Noto Sans CJK JP"; plt.rcParams["axes.unicode_minus"] = False
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 3.4))
ks = np.arange(9); ax1.plot(ks, p3, 'o-', color="#4C6EF5", label="numpy")
kk = np.linspace(0, 8, 200); ax1.plot(kk, np.sin((2*kk+1)*np.arcsin(1/np.sqrt(8)))**2, '--', color="#868E96", label="sin²((2k+1)θ)")
ax1.set_xlabel("반복 k"); ax1.set_ylabel("P(표시 원소)"); ax1.set_title("Grover n=3: 지나치면 줄어든다"); ax1.legend(frameon=False, fontsize=8)
ns = [3, 4, 6, 8, 10]; ax2.plot(ns, [r[1] for r in rows], 'o-', color="#4C6EF5", label="최적 k"); ax2.plot(ns, [np.pi/4*np.sqrt(2**m) for m in ns], '--', color="#868E96", label="π/4·√N")
ax2.set_xlabel("큐비트 수 n"); ax2.set_ylabel("반복 횟수"); ax2.set_yscale("log"); ax2.set_title("최적 반복 ∝ √N"); ax2.legend(frameon=False, fontsize=8)
for lam, p in noisy.items(): ax3.plot(range(kmax+1), p, 'o-', ms=3, label=f"λ={lam}")
ax3.set_xlabel("반복 k"); ax3.set_ylabel("P(표시 원소)"); ax3.set_title("위상 감쇠가 반복마다 깎는다 (n=6)"); ax3.legend(frameon=False, fontsize=8)
plt.tight_layout(); plt.savefig("ch09_results.png", dpi=150); print("7) saved ch09_results.png")

# ── 8. 게이트 비용 기록 (2큐비트 게이트 / 깊이 / T-count) ─────────────
# 4절의 n = 3, 표시 |101> 회로를 게이트 리스트로 다시 짜서 위의 행렬 G = D·O 와 같은지 확인한 뒤 센다(책 앞 「기호와 비용 단위」의 규약).
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tools"))
import gate_cost as gc
orc = [gc.G("x", 1), gc.G("ccz", 0, 1, 2), gc.G("x", 1)]                                   # |101> 에만 −1
dif = [gc.G("h", q) for q in range(3)] + [gc.G("x", q) for q in range(3)] + [gc.G("ccz", 0, 1, 2)] \
    + [gc.G("x", q) for q in range(3)] + [gc.G("h", q) for q in range(3)]                  # 2|s><s| − I (전역 위상 제외)
it = orc + dif
assert gc.same_up_to_phase(gc.unitary(it, 3), diffusion_formula(3) @ oracle([5], 3))
full = [gc.G("h", q) for q in range(3)] + it + it
psi_g = gc.apply(full, ket0n(3), 3)
assert abs(abs(psi_g[5]) ** 2 - p3[2]) < 1e-12
print("\n8) 게이트 비용 기록 (n = 3, 표시 |101>; 전결합, CCZ = CX 6 + T 7)")
print("   반복 1회:        ", gc.fmt(it))
print("   전체(H 층 + k=2):", gc.fmt(full))
