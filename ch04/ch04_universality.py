"""4장 범용성과 게이트 분해 — {H,T} 근사, Toffoli 분해, 기저 게이트 변환을 numpy로 검증한다.
실행: python ch04_universality.py
환경: Python 3.13, numpy 2.5 (SDK 없이 동작). 큐비트 순서: 빅 엔디안(3장 규약).
"""
import numpy as np
from itertools import product

np.set_printoptions(precision=4, suppress=True)
rng = np.random.default_rng(seed=7)

I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
S = np.array([[1, 0], [0, 1j]], dtype=complex)
T = np.array([[1, 0], [0, np.exp(1j*np.pi/4)]], dtype=complex)
Tdg = T.conj().T
SX = np.array([[1+1j, 1-1j], [1-1j, 1+1j]], dtype=complex) / 2
def Rx(t): return np.array([[np.cos(t/2), -1j*np.sin(t/2)], [-1j*np.sin(t/2), np.cos(t/2)]])
def Ry(t): return np.array([[np.cos(t/2), -np.sin(t/2)], [np.sin(t/2), np.cos(t/2)]], dtype=complex)
def Rz(t): return np.array([[np.exp(-1j*t/2), 0], [0, np.exp(1j*t/2)]])

def kron(*ops):
    out = np.array([[1]], dtype=complex)
    for o in ops: out = np.kron(out, o)
    return out

def dist(A, B):
    """전역 위상 무시 거리: 1 - |Tr(A†B)|/d"""
    d = A.shape[0]
    return 1 - abs(np.trace(A.conj().T @ B)) / d

def su2_angle(A, B):
    """두 SU(2) 원소 사이의 회전각 (전역 위상 무시). cos(θ/2) = |Tr(A†B)|/2"""
    return 2 * np.arccos(min(1.0, abs(np.trace(A.conj().T @ B)) / 2))

# ── 1. HTH T 가 무리수 각 회전임을 수치로 확인 ────────────────────
U0 = H @ T @ H @ T          # 범용성의 핵심 원소. 회전축·각 계산
# U = e^{iα} exp(-iθ n·σ/2) 에서 cos(θ/2) = |Tr U|/2
theta = 2*np.arccos(abs(np.trace(U0))/2)
print(f"1) HTHT 회전각 θ = {theta:.6f} rad = {theta/np.pi:.6f} π  (π의 무리수 배 → 거듭제곱이 구를 조밀하게 덮는다)")
print("   cos(θ/2) = cos²(π/8) =", round(np.cos(np.pi/8)**2, 6), " 실측:", round(abs(np.trace(U0))/2, 6))

# ── 2. {H, T} 단어 집합으로 임의 게이트 근사: 길이 vs 최소 오차 ─────
target = Rx(1.0) @ Ry(0.7)   # 임의의 목표 게이트
gates = {'H': H, 'T': T}
print("2) 목표 U = Rx(1.0)Ry(0.7) 를 {H,T} 단어로 근사 (전수 탐색)")
print("   길이≤L | 단어 수(길이 L) | 최소 거리(1-|Tr|/2) | 최소 회전각 오차(rad)")
best_word_overall = None
cum = (np.inf, None, None)
for L in range(1, 15):
    best = (np.inf, None)
    for word in product('HT', repeat=L):
        # HH = I 인 단어는 건너뛰어 탐색을 줄인다
        if any(word[i] == word[i+1] == 'H' for i in range(L-1)): continue
        M = I.copy()
        for g in word: M = M @ gates[g]
        d = dist(M, target)
        if d < best[0]: best = (d, ''.join(word), M)
    if best[0] < cum[0]: cum = best
    ang = su2_angle(cum[2], target)
    print(f"   {L:4d} | {2**L:7d} | {cum[0]:.6f} | {ang:.4f}")
    best_word_overall = cum
print("   최선의 단어:", best_word_overall[1])
print("   → 길이가 늘수록 오차가 줄지만 느리다. Solovay–Kitaev 는 이것을 O(log^c(1/ε)) 로 만든다.")

# ── 3. Toffoli 분해: 6 CNOT + 단일 큐비트 게이트 (N&C Fig. 4.9와 동치인 S 없는 15게이트 판) ────
def ccx_matrix():
    M = np.eye(8, dtype=complex); M[6, 6] = M[7, 7] = 0; M[6, 7] = M[7, 6] = 1
    return M
CCX = ccx_matrix()
def cnot(ctrl, tgt, n=3):   # 빅 엔디안 n큐비트 CNOT
    P0 = np.diag([1, 0]).astype(complex); P1 = np.diag([0, 1]).astype(complex)
    ops0 = [I]*n; ops0[ctrl] = P0
    ops1 = [I]*n; ops1[ctrl] = P1; ops1[tgt] = X
    return kron(*ops0) + kron(*ops1)
def single(g, q, n=3):
    ops = [I]*n; ops[q] = g; return kron(*ops)
# 회로 순서(왼쪽이 먼저): H(2) CX(1,2) Tdg(2) CX(0,2) T(2) CX(1,2) Tdg(2) CX(0,2) T(1) T(2) H(2) CX(0,1) T(0) Tdg(1) CX(0,1)
seq = [single(H,2), cnot(1,2), single(Tdg,2), cnot(0,2), single(T,2), cnot(1,2), single(Tdg,2), cnot(0,2),
       single(T,1), single(T,2), single(H,2), cnot(0,1), single(T,0), single(Tdg,1), cnot(0,1)]
M = np.eye(8, dtype=complex)
for g in seq: M = g @ M         # 시간 순서대로 왼쪽에 곱한다
print("3) Toffoli = 6 CNOT + 7 T/T† + 2 H (N&C Fig.4.9와 동치인 15게이트 판, Qiskit ccx 정의):", np.allclose(M, CCX))
print("   CNOT 6개, T/T† 7개, H 2개 — N&C 그림 4.9(S 하나를 더 쓴 16게이트)에서 S를 T 배치로 흡수한 판")

# ── 4. 임의 단일 큐비트 게이트 → IBM 기저 {rz, sx, x} ──────────────
def zyz(Uin):
    V = Uin / np.sqrt(np.linalg.det(Uin))
    gamma = 2*np.arctan2(abs(V[1,0]), abs(V[0,0]))
    if np.isclose(np.sin(gamma/2), 0): beta, delta = 2*np.angle(V[1,1]), 0.0
    elif np.isclose(np.cos(gamma/2), 0): beta, delta = 2*np.angle(V[1,0]), 0.0
    else: beta = np.angle(V[1,1]) + np.angle(V[1,0]); delta = np.angle(V[1,1]) - np.angle(V[1,0])
    return beta, gamma, delta
def random_unitary():
    A = rng.normal(size=(2,2)) + 1j*rng.normal(size=(2,2))
    Q, R = np.linalg.qr(A); return Q @ np.diag(np.diag(R)/abs(np.diag(R)))
U = random_unitary()
b, g, d = zyz(U)
# Ry(γ) = Rz(-π/2) SX Rz(π-γ)... 간단히: Ry(γ) ≅ Rz(-π/2) · Rx(γ) · Rz(π/2),  Rx(γ) ≅ SX · Rz(γ) · SX† ...
# 표준 결과: U ≅ Rz(β+π) · SX · Rz(γ+π) · SX · Rz(δ)  (Qiskit 'ZSXZSXZ' 분해)
M_isa = Rz(b + np.pi) @ SX @ Rz(g + np.pi) @ SX @ Rz(d)
print(f"4) ZSXZSXZ 분해: U ≅ Rz({b+np.pi:.3f})·SX·Rz({g+np.pi:.3f})·SX·Rz({d:.3f})  거리 = {dist(M_isa, U):.2e}")
print("   → 임의 단일 큐비트 게이트 = SX 2개 + Rz 3개. Rz 는 가상 게이트(2장)이므로 실제 비용은 SX 2개.")

# ── 5. 연결성 제약: 선형 사슬에서 CNOT(0,2) 는 SWAP 이 필요 ─────────
CX02 = cnot(0,2)
SWAP12 = cnot(1,2) @ cnot(2,1) @ cnot(1,2)
routed = SWAP12 @ cnot(0,1) @ SWAP12          # SWAP(1,2) · CX(0,1) · SWAP(1,2)
print("5) CX(0,2) = SWAP(1,2)·CX(0,1)·SWAP(1,2):", np.allclose(routed, CX02), "→ CNOT 1개가 7개가 된다")

# ── 6. 그림: 근사 오차 vs 단어 길이 ─────────────────────────────
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Noto Sans CJK JP"; plt.rcParams["axes.unicode_minus"] = False
Ls, errs = [], []
best = np.inf
for L in range(1, 15):
    for word in product('HT', repeat=L):
        if any(word[i] == word[i+1] == 'H' for i in range(L-1)): continue
        Mw = I.copy()
        for gg in word: Mw = Mw @ gates[gg]
        best = min(best, su2_angle(Mw, target))
    Ls.append(L); errs.append(best)
fig, ax = plt.subplots(figsize=(6, 3.4))
ax.semilogy(Ls, errs, 'o-', color="#4C6EF5")
ax.set_xlabel("{H,T} 단어 길이 ≤ L"); ax.set_ylabel("최소 회전각 오차 (rad)"); ax.set_title("전수 탐색: 길이가 늘면 오차는 줄지만 느리게")
plt.tight_layout(); plt.savefig("ch04_results.png", dpi=150); print("6) saved ch04_results.png")

# ── 연습문제 검산 (8절) ──────────────────────────────────────────
# 8절 해답이 인용하는 수치·행렬을 전부 다시 계산한다. 위쪽 코드·난수(rng)·그림은 건드리지 않는다.
print("── 연습문제 검산 (8절) ──")

# 문제 1: {H, S} 가 만드는 게이트는 전역 위상을 무시하면 24개 (너비 우선 탐색)
def ex_canon(Min):
    """위상 정규화: 첫 번째 0이 아닌 성분을 양의 실수로 맞춘 뒤 반올림해 사전 키로 쓴다"""
    f = Min.flatten(); k = int(np.flatnonzero(abs(f) > 1e-9)[0])
    return tuple((np.round(f * (abs(f[k]) / f[k]), 6) + 0.0).tolist())
def ex_bfs(gens, cap=5000):
    """I 에서 출발해 생성원을 계속 곱한다. 새 원소가 안 나오면 멈춤(닫힘=True), cap 을 넘으면 포기(닫힘=False)"""
    seen = {ex_canon(I): I}; frontier = [I]; depth = 0
    while frontier and len(seen) <= cap:
        nxt = []
        for Mb in frontier:
            for G in gens:
                N = G @ Mb; key = ex_canon(N)
                if key not in seen: seen[key] = N; nxt.append(N)
        frontier = nxt
        if nxt: depth += 1
    return seen, depth, not frontier
ex_cl, ex_depth, ex_closed = ex_bfs([H, S])
print(f"[연습 1] {{H, S}} 너비 우선 탐색: 서로 다른 게이트(전역 위상 무시) = {len(ex_cl)}개, 닫힘 = {ex_closed} (단어 길이 {ex_depth} 까지에서 다 나오고 그 뒤로 새 원소 없음)")
ex_ht, ex_depth_ht, ex_closed_ht = ex_bfs([H, T], cap=5000)
print(f"[연습 1] {{H, T}} 같은 탐색: 단어 길이 {ex_depth_ht} 까지에서 이미 {len(ex_ht)}개 (상한 5000 초과), 닫힘 = {ex_closed_ht} → 끝나지 않는다")
def ex_so3(Uin):
    """단일 큐비트 게이트가 Bloch 구에 일으키는 3×3 회전행렬 R_ij = Tr(σ_i U σ_j U†)/2"""
    return np.array([[np.trace(a @ Uin @ b @ Uin.conj().T).real / 2 for b in (X, Y, Z)] for a in (X, Y, Z)])
ex_rots = [ex_so3(Mb) for Mb in ex_cl.values()]
ex_signed_perm = all(np.allclose(np.sort(abs(R), axis=1), [[0, 0, 1]]*3) and np.allclose(np.sort(abs(R), axis=0).T, [[0, 0, 1]]*3)
                     and np.isclose(np.linalg.det(R), 1) for R in ex_rots)
ex_distinct = len({tuple((np.round(R).astype(int)).flatten().tolist()) for R in ex_rots})
print(f"[연습 1] 24개 각각의 Bloch 회전이 좌표축을 좌표축으로 보내는 고유 회전(det=+1)? {ex_signed_perm},  서로 다른 회전 = {ex_distinct}개"
      " = 정육면체 회전군의 크기 24")

# 문제 2: CNOT = (I⊗H)·CZ·(I⊗H), Toffoli 의 2큐비트 게이트 수는 그대로 6
CNOT = cnot(0, 1, n=2); CZ = np.diag([1, 1, 1, -1]).astype(complex)   # 2큐비트 행렬(3장 규약: 제어 q0, 표적 q1)
print("[연습 2] np.allclose(kron(I,H) @ CZ @ kron(I,H), CNOT) =", np.allclose(kron(I, H) @ CZ @ kron(I, H), CNOT))
def ex_cz(a, b_, n=3):      # 빅 엔디안 n큐비트 CZ: 두 큐비트가 모두 1 인 기저에만 -1
    return np.diag([-1 if (i >> (n-1-a)) & 1 and (i >> (n-1-b_)) & 1 else 1 for i in range(2**n)]).astype(complex)
# 3절 Toffoli 회로(시간 순서)를 기호로 다시 적고, CX(c,t) 마다 H(t)·CZ(c,t)·H(t) 로 바꾼다
ex_circ = [('h', 2), ('cx', 1, 2), ('tdg', 2), ('cx', 0, 2), ('t', 2), ('cx', 1, 2), ('tdg', 2), ('cx', 0, 2),
           ('t', 1), ('t', 2), ('h', 2), ('cx', 0, 1), ('t', 0), ('tdg', 1), ('cx', 0, 1)]
ex_1q = {'h': H, 't': T, 'tdg': Tdg}
ex_M_cx = np.eye(8, dtype=complex); ex_M_cz = np.eye(8, dtype=complex); ex_n_cz = 0
for ex_op in ex_circ:
    if ex_op[0] == 'cx':
        ex_M_cx = cnot(ex_op[1], ex_op[2]) @ ex_M_cx
        ex_M_cz = single(H, ex_op[2]) @ ex_cz(ex_op[1], ex_op[2]) @ single(H, ex_op[2]) @ ex_M_cz; ex_n_cz += 1
    else:
        ex_M_cx = single(ex_1q[ex_op[0]], ex_op[1]) @ ex_M_cx; ex_M_cz = single(ex_1q[ex_op[0]], ex_op[1]) @ ex_M_cz
print(f"[연습 2] Toffoli 의 CNOT 을 전부 H·CZ·H 로 바꾼 회로 = Toffoli ? {np.allclose(ex_M_cz, CCX)}  (CNOT 판도 {np.allclose(ex_M_cx, CCX)}),"
      f"  2큐비트 게이트 = CZ {ex_n_cz}개")
print("[연습 2] 연속된 CNOT 쌍: CX(0,2)·CX(1,2) = H(2)·CZ(0,2)·CZ(1,2)·H(2) (가운데 H·H = I 상쇄):",
      np.allclose(cnot(0, 2) @ cnot(1, 2), single(H, 2) @ ex_cz(0, 2) @ ex_cz(1, 2) @ single(H, 2)))
print(f"[연습 2] H 는 √X·Rz 로 녹는다: dist(Rz(π/2)·SX·Rz(π/2), H) = {abs(dist(Rz(np.pi/2) @ SX @ Rz(np.pi/2), H)):.1e}")

# 문제 3: 제어-Rz(t) = (I⊗A)·CNOT·(I⊗B)·CNOT,  A = Rz(t/2), B = Rz(-t/2), C = I
ex_t = 0.9
ex_A, ex_B, ex_C = Rz(ex_t/2), Rz(-ex_t/2), I
print("[연습 3] ABC = I :", np.allclose(ex_A @ ex_B @ ex_C, I), "  X Rz(-t/2) X = Rz(t/2) :", np.allclose(X @ ex_B @ X, ex_A),
      "  A X B X = Rz(t) :", np.allclose(ex_A @ X @ ex_B @ X, Rz(ex_t)), " (t = 0.9, α = 0)")
ex_CRz = np.diag([1, 1, np.exp(-1j*ex_t/2), np.exp(1j*ex_t/2)])
print("[연습 3] np.allclose(kron(I,A) @ CNOT @ kron(I,B) @ CNOT, diag(1, 1, e^{-it/2}, e^{it/2})) =",
      np.allclose(kron(I, ex_A) @ CNOT @ kron(I, ex_B) @ CNOT, ex_CRz), " (t = 0.9)")
print("[연습 3] 같은 행렬 = |0><0|⊗I + |1><1|⊗Rz(t) (제어-Rz 의 정의):",
      np.allclose(ex_CRz, kron(np.diag([1, 0]), I) + kron(np.diag([0, 1]), Rz(ex_t))))
ex_th = 0.37
ex_ZZ = kron(Z, Z)
ex_expZZ = np.cos(ex_th) * np.eye(4) - 1j * np.sin(ex_th) * ex_ZZ      # (Z⊗Z)² = I 이므로 exp(-iθ Z⊗Z) = cosθ·I - i sinθ·Z⊗Z
print("[연습 3] CNOT·(I⊗Rz(2θ))·CNOT = exp(-iθ Z⊗Z)  (θ = 0.37):", np.allclose(CNOT @ kron(I, Rz(2*ex_th)) @ CNOT, ex_expZZ))
ex_lam = 0.9
ex_CP = np.diag([1, 1, 1, np.exp(1j*ex_lam)])
ex_CRz_lam = np.diag([1, 1, np.exp(-1j*ex_lam/2), np.exp(1j*ex_lam/2)])
print("[연습 3] CP(λ) = CRz(λ)·(P(λ/2)⊗I)  (λ = 0.9):", np.allclose(ex_CP, ex_CRz_lam @ kron(np.diag([1, np.exp(1j*ex_lam/2)]), I)))
