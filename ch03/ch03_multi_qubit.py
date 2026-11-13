"""3장 다중 큐비트 게이트와 얽힘 — 텐서곱·CNOT·Bell·GHZ·부분추적·CHSH를 numpy로 검증한다.
실행: python ch03_multi_qubit.py
환경: Python 3.13, numpy 2.5 (SDK 없이 동작)
규약: 이 파일은 N&C 교과서 순서(빅 엔디안: 왼쪽 큐비트 = 상위 비트)를 쓴다.
      Qiskit 순서(리틀 엔디안)로 바꾸는 함수 to_qiskit_order() 를 끝에 둔다.
"""
import os
import numpy as np
from itertools import product

np.set_printoptions(precision=4, suppress=True)
rng = np.random.default_rng(seed=7)

I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
def Ry(t): return np.array([[np.cos(t/2), -np.sin(t/2)], [np.sin(t/2), np.cos(t/2)]], dtype=complex)

ket0 = np.array([1, 0], dtype=complex); ket1 = np.array([0, 1], dtype=complex)
def kron(*ops):
    out = np.array([[1]], dtype=complex)
    for o in ops: out = np.kron(out, o)
    return out
def ket(bits):  # '01' -> |0>⊗|1>
    v = np.array([1], dtype=complex)
    for b in bits: v = np.kron(v, ket1 if b == '1' else ket0)
    return v

# ── 1. 2큐비트 게이트 (빅 엔디안: 큐비트0 = 왼쪽 = 상위 비트) ─────────
CNOT = np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]], dtype=complex)   # control=q0, target=q1
CZ   = np.diag([1,1,1,-1]).astype(complex)
SWAP = np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]], dtype=complex)
P0 = np.outer(ket0, ket0.conj()); P1 = np.outer(ket1, ket1.conj())
CNOT_built = kron(P0, I) + kron(P1, X)        # |0><0|⊗I + |1><1|⊗X
print("1) CNOT = |0><0|⊗I + |1><1|⊗X :", np.allclose(CNOT, CNOT_built))
print("   CZ = (I⊗H) CNOT (I⊗H) :", np.allclose(CZ, kron(I,H) @ CNOT @ kron(I,H)))
print("   CZ 는 대칭(제어·표적 교환 불변):", np.allclose(CZ, SWAP @ CZ @ SWAP))
print("   CNOT 은 비대칭:", np.allclose(CNOT, SWAP @ CNOT @ SWAP))
print("   SWAP = CNOT01·CNOT10·CNOT01 :", np.allclose(SWAP, CNOT @ (SWAP@CNOT@SWAP) @ CNOT))
print("   (H⊗H) CNOT (H⊗H) = CNOT 뒤집힘 :", np.allclose(kron(H,H) @ CNOT @ kron(H,H), SWAP@CNOT@SWAP))
assert np.allclose(CNOT, CNOT_built), "5절 게이트 항등식: CNOT = |0><0|⊗I + |1><1|⊗X"
assert np.allclose(CZ, kron(I,H) @ CNOT @ kron(I,H)), "5절 게이트 항등식: CZ = (I⊗H)CNOT(I⊗H)"
assert np.allclose(CZ, SWAP @ CZ @ SWAP), "5절 게이트 항등식: CZ 는 대칭"
assert not np.allclose(CNOT, SWAP @ CNOT @ SWAP), "5절 게이트 항등식: CNOT 은 비대칭(False)"
assert np.allclose(SWAP, CNOT @ (SWAP@CNOT@SWAP) @ CNOT), "5절 게이트 항등식: SWAP = CNOT01·CNOT10·CNOT01"
assert np.allclose(kron(H,H) @ CNOT @ kron(H,H), SWAP@CNOT@SWAP), "5절 게이트 항등식: (H⊗H)CNOT(H⊗H) = 뒤집힌 CNOT"

# ── 2. Bell 상태 ─────────────────────────────────────────────────
phi_plus = CNOT @ kron(H, I) @ ket('00')
print("2) |Φ+> = CNOT (H⊗I)|00> =", phi_plus, " = (|00>+|11>)/√2")
def probs(psi): return np.abs(psi)**2
print("   P(00,01,10,11) =", probs(phi_plus))

def sample(psi, shots):
    n = int(np.log2(len(psi)))
    idx = rng.choice(len(psi), size=shots, p=probs(psi))
    labels = [format(i, f'0{n}b') for i in idx]
    return {k: labels.count(k) for k in sorted(set(labels))}
c = sample(phi_plus, 1024)
print("   1024샷:", c, "-> 두 비트가 항상 같다")
assert np.allclose(phi_plus, [0.7071, 0, 0, 0.7071], atol=5e-5), "5절 Bell 상태: 진폭 (0.7071, 0, 0, 0.7071)"
assert np.allclose(probs(phi_plus), [0.5, 0, 0, 0.5]), "5절 Bell 상태: 확률 (0.5, 0, 0, 0.5)"
assert c == {"00": 517, "11": 507}, "5절 Bell 상태: 1024샷 = 00이 517회, 11이 507회, 01·10은 0회"

# 분리 가능성 검사: 2큐비트 순수 상태가 곱상태인지 (진폭 행렬의 랭크)
def is_product(psi):
    M = psi.reshape(2, 2)
    return np.linalg.matrix_rank(M, tol=1e-9) == 1
print("   |Φ+> 곱상태?", is_product(phi_plus), " |+>|0> 곱상태?", is_product(kron(H,I)@ket('00')))
assert np.linalg.matrix_rank(phi_plus.reshape(2, 2), tol=1e-9) == 2, "5절 Bell 상태: 진폭 행렬의 랭크 2(얽힘)"
assert is_product(kron(H,I)@ket('00')), "5절 Bell 상태: |+>|0> 는 랭크 1(곱상태)"

# ── 3. 부분 추적과 순수도 ────────────────────────────────────────
def partial_trace_B(rho):   # 2큐비트 밀도행렬에서 B를 추적 → A의 밀도행렬
    r = rho.reshape(2,2,2,2)
    return np.einsum('ijkj->ik', r)
rho = np.outer(phi_plus, phi_plus.conj())
rhoA = partial_trace_B(rho)
print("3) ρ_A = Tr_B|Φ+><Φ+| =\n", rhoA)
print("   순수도 Tr(ρ_A²) =", np.trace(rhoA@rhoA).real, "(1이면 순수, 0.5면 최대 혼합)")
rho_prod = np.outer(kron(H,I)@ket('00'), (kron(H,I)@ket('00')).conj())
print("   곱상태 |+>|0> 의 ρ_A 순수도 =", np.trace(partial_trace_B(rho_prod)@partial_trace_B(rho_prod)).real)
# 얽힘 엔트로피
ev = np.linalg.eigvalsh(rhoA); ev = ev[ev > 1e-12]
print("   얽힘 엔트로피 S(ρ_A) =", float(-(ev*np.log2(ev)).sum()), "bit")
assert np.allclose(rhoA, np.diag([0.5, 0.5])), "5절 부분 추적 표: |Φ+> 의 ρ_A = diag(0.5, 0.5)"
assert abs(np.trace(rhoA@rhoA).real - 0.5000) < 5e-5, "5절 부분 추적 표: |Φ+> 의 순수도 0.5000"
assert abs(float(-(ev*np.log2(ev)).sum()) - 1.000) < 5e-4, "5절 부분 추적 표: |Φ+> 의 S(ρ_A) = 1.000 bit"
assert np.allclose(partial_trace_B(rho_prod), np.outer(H @ ket0, (H @ ket0).conj())), "5절 부분 추적 표: |+>|0> 의 ρ_A = |+><+|"
assert abs(np.trace(partial_trace_B(rho_prod)@partial_trace_B(rho_prod)).real - 1.0000) < 5e-5, "5절 부분 추적 표: |+>|0> 의 순수도 1.0000"

# ── 4. GHZ (3큐비트) ─────────────────────────────────────────────
CNOT01 = kron(CNOT, I)                       # control q0, target q1
CNOT12 = kron(I, CNOT)                       # control q1, target q2
ghz = CNOT12 @ CNOT01 @ kron(H, I, I) @ ket('000')
c_ghz = sample(ghz, 1024)
print("4) GHZ =", ghz, " 1024샷:", c_ghz)
assert np.allclose(ghz, [0.7071, 0, 0, 0, 0, 0, 0, 0.7071], atol=5e-5), "5절 GHZ: 진폭은 |000>, |111> 에만 0.7071"
assert c_ghz == {"000": 520, "111": 504}, "5절 GHZ: 1024샷 = 000이 520회, 111이 504회"

# ── 5. CHSH: 얽힘은 고전 상관으로 흉내낼 수 없다 ─────────────────
def expect(psi, op): return np.vdot(psi, op @ psi).real
def obs(theta): return np.cos(theta)*Z + np.sin(theta)*X      # 측정 축 (xz 평면)
a, a2, b, b2 = 0, np.pi/2, np.pi/4, -np.pi/4
E = lambda t1, t2: expect(phi_plus, kron(obs(t1), obs(t2)))
S = E(a,b) + E(a,b2) + E(a2,b) - E(a2,b2)     # CHSH: 한 항만 부호가 다르다
print(f"5) CHSH S = {S:.4f}  (고전 한계 2, 양자 최대 2√2 = {2*np.sqrt(2):.4f})")
print("   각 항 E(a,b), E(a,b'), E(a',b), E(a',b') =", [float(round(E(a,b),4)), float(round(E(a,b2),4)), float(round(E(a2,b),4)), float(round(E(a2,b2),4))])
rho_cl = 0.5*np.outer(ket('00'),ket('00')) + 0.5*np.outer(ket('11'),ket('11'))
Ecl = lambda t1, t2: np.trace(rho_cl @ kron(obs(t1), obs(t2))).real
S_cl = Ecl(a,b) + Ecl(a,b2) + Ecl(a2,b) - Ecl(a2,b2)
print(f"   같은 각도에서 고전 혼합 상태의 S = {S_cl:.4f}  (≤ 2)")
assert abs(S - 2.8284) < 5e-5 and abs(S - 2*np.sqrt(2)) < 1e-9, "5절 CHSH 표: Bell 의 S = 2.8284 = 2√2"
assert np.allclose([E(a,b), E(a,b2), E(a2,b), E(a2,b2)], [0.7071, 0.7071, 0.7071, -0.7071], atol=5e-5), "5절 CHSH 표: Bell 열의 네 항"
assert np.allclose([Ecl(a,b), Ecl(a,b2), Ecl(a2,b), Ecl(a2,b2)], [0.7071, 0.7071, 0, 0], atol=5e-5), "5절 CHSH 표: 고전 혼합 열의 네 항"
assert abs(S_cl - 1.4142) < 5e-5 and S_cl <= 2, "5절 CHSH 표: 고전 혼합의 S = 1.4142"

# ── 6. 고전 상관(혼합 상태)과의 비교: 같은 Z-상관, 다른 X-상관 ───
rho_cl = 0.5*np.outer(ket('00'),ket('00')) + 0.5*np.outer(ket('11'),ket('11'))
def expect_rho(rho, op): return np.trace(rho@op).real
print("6) ZZ 상관: Bell", expect(phi_plus, kron(Z,Z)), " 고전혼합", expect_rho(rho_cl, kron(Z,Z)))
print("   XX 상관: Bell", expect(phi_plus, kron(X,X)), " 고전혼합", expect_rho(rho_cl, kron(X,X)))
assert abs(expect(phi_plus, kron(Z,Z)) - 1.0) < 1e-9 and abs(expect_rho(rho_cl, kron(Z,Z)) - 1.0) < 1e-9, "5절·7절: ZZ 상관은 Bell·고전 혼합 둘 다 1.0"
assert abs(expect(phi_plus, kron(X,X)) - 1.0) < 1e-9 and abs(expect_rho(rho_cl, kron(X,X))) < 1e-9, "5절·7절: XX 상관은 Bell 1.0, 고전 혼합 0"

# ── 7. 큐비트 순서: Qiskit(리틀 엔디안)으로 변환 ──────────────────
def to_qiskit_order(psi):
    n = int(np.log2(len(psi)))
    return psi.reshape([2]*n).transpose(list(range(n))[::-1]).reshape(-1)
psi = kron(H, I) @ ket('00')               # 빅엔디안: q0=|+>, q1=|0> → 진폭 [.707, 0, .707, 0]
print("7) |+>|0> 빅엔디안 진폭:", psi)
print("   같은 상태를 Qiskit 순서로:", to_qiskit_order(psi), " (Statevector 출력이 이렇게 나온다)")
assert np.allclose(psi, [0.707, 0, 0.707, 0], atol=5e-4), "5절 큐비트 순서: |+>|0> 빅 엔디안 진폭 (0.707, 0, 0.707, 0)"
assert np.allclose(to_qiskit_order(psi), [0.707, 0.707, 0, 0], atol=5e-4), "5절 큐비트 순서: Qiskit 순서 (0.707, 0.707, 0, 0)"
print("   Qiskit에서 '01' 은 q1=0, q0=1 을 뜻한다 (오른쪽 끝이 큐비트 0)")

# ── 8. 그림 ──────────────────────────────────────────────────────
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = ["Noto Sans CJK JP", "DejaVu Sans"]; plt.rcParams["axes.unicode_minus"] = False
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 3.4))
keys = ['00','01','10','11']; vals = [c.get(k,0) for k in keys]
ax1.bar(keys, vals, color="#4C6EF5"); ax1.set_title("|Φ+⟩ 측정 1024샷"); ax1.set_ylabel("횟수")
for i,v in enumerate(vals): ax1.text(i, v+10, str(v), ha="center")
thetas = np.linspace(0, np.pi, 200)
ax2.plot(thetas, [expect(phi_plus, kron(obs(np.pi/2), obs(t))) for t in thetas], color="#4C6EF5", label="Bell |Φ+⟩")
ax2.plot(thetas, [expect_rho(rho_cl, kron(obs(np.pi/2), obs(t))) for t in thetas], color="#868E96", ls="--", label="고전 혼합")
ax2.set_xticks([0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi]); ax2.set_xticklabels(["0","π/4","π/2","3π/4","π"])
ax2.set_xlabel("B의 측정 각 θ (A는 X축)"); ax2.set_ylabel("상관 E(π/2, θ)"); ax2.set_title("상관의 각도 의존성"); ax2.legend(frameon=False, fontsize=9)
plt.tight_layout(); plt.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)), "ch03_results.png"), dpi=150); print("8) saved ch03_results.png")

# ── 연습문제 검산 (8절) ──────────────────────────────────────────
# 8절 해답이 인용하는 수치·행렬을 전부 다시 계산한다. 위쪽 코드·난수(rng)·그림은 건드리지 않는다.
print("── 연습문제 검산 (8절) ──")
ex_r2 = 1 / np.sqrt(2)

# 문제 1: 나머지 Bell 상태 세 개와 정규직교성
ex_phi_m = kron(Z, I) @ phi_plus; ex_psi_p = kron(X, I) @ phi_plus; ex_psi_m = kron(Z @ X, I) @ phi_plus
print("[연습 1] |Φ-> = (Z⊗I)|Φ+> = (|00>-|11>)/√2 :", np.allclose(ex_phi_m, (ket('00') - ket('11')) * ex_r2))
print("[연습 1] |Ψ+> = (X⊗I)|Φ+> = (|01>+|10>)/√2 :", np.allclose(ex_psi_p, (ket('01') + ket('10')) * ex_r2))
print("[연습 1] |Ψ-> = (ZX⊗I)|Φ+> = (|01>-|10>)/√2 :", np.allclose(ex_psi_m, (ket('01') - ket('10')) * ex_r2))
assert np.allclose(ex_phi_m, (ket('00') - ket('11')) * ex_r2), "8절 문제 1: |Φ-> = (Z⊗I)|Φ+>"
assert np.allclose(ex_psi_p, (ket('01') + ket('10')) * ex_r2), "8절 문제 1: |Ψ+> = (X⊗I)|Φ+>"
assert np.allclose(ex_psi_m, (ket('01') - ket('10')) * ex_r2), "8절 문제 1: |Ψ-> = (ZX⊗I)|Φ+>"
ex_named = {'Φ+': phi_plus, 'Φ-': ex_phi_m, 'Ψ+': ex_psi_p, 'Ψ-': ex_psi_m}
ex_map = {}
for ex_ab in ('00', '01', '10', '11'):
    ex_v = CNOT @ kron(H, I) @ ket(ex_ab)
    ex_map[ex_ab] = [name for name, w in ex_named.items() if np.allclose(ex_v, w)]   # 위상까지 정확히 같은 것만
print("[연습 1] 일반형 CNOT(H⊗I)|ab> :", ", ".join(f"|{k}>→{'/'.join(v)}" for k, v in ex_map.items()))
ex_B = np.column_stack([CNOT @ kron(H, I) @ ket(ex_ab) for ex_ab in ('00', '01', '10', '11')])
print("[연습 1] np.allclose(B.conj().T @ B, np.eye(4)) =", np.allclose(ex_B.conj().T @ ex_B, np.eye(4)))
print("[연습 1] Bell 측정: (H⊗I)·CNOT 이 네 Bell 상태를 |00>,|01>,|10>,|11> 로 되돌린다:",
      np.allclose(kron(H, I) @ CNOT @ ex_B, np.eye(4)))
assert ex_map == {'00': ['Φ+'], '01': ['Ψ+'], '10': ['Φ-'], '11': ['Ψ-']}, "8절 문제 1: CNOT(H⊗I)|ab> 가 네 Bell 상태를 고른다"
assert np.allclose(ex_B.conj().T @ ex_B, np.eye(4)), "8절 문제 1: 네 Bell 상태는 정규직교(B†B = I)"
assert np.allclose(kron(H, I) @ CNOT @ ex_B, np.eye(4)), "8절 문제 1: Bell 측정 = 생성 회로를 거꾸로"

# 문제 2: cosθ|00> + sinθ|11> 의 얽힘 엔트로피 = H₂(cos²θ)
def ex_state(theta): return np.cos(theta) * ket('00') + np.sin(theta) * ket('11')
def ex_entropy(psi):
    """3절과 같은 방법: partial_trace_B 뒤 고유값 엔트로피"""
    rA = partial_trace_B(np.outer(psi, psi.conj()))
    w = np.linalg.eigvalsh(rA); w = w[w > 1e-12]
    return float(-(w * np.log2(w)).sum()) + 0.0, rA
def ex_H2(p): return float(-p*np.log2(p) - (1-p)*np.log2(1-p))
ex_t = np.pi / 8
ex_S, ex_rA = ex_entropy(ex_state(ex_t))
print("[연습 2] θ=π/8: ρ_A = diag(cos²θ, sin²θ) :", np.allclose(ex_rA, np.diag([np.cos(ex_t)**2, np.sin(ex_t)**2])),
      f" 대각 = ({ex_rA[0, 0].real:.4f}, {ex_rA[1, 1].real:.4f})")
print(f"[연습 2] θ=π/8: cos²(π/8) = {np.cos(ex_t)**2:.4f},  S(ρ_A) = {ex_S:.4f} bit,  H₂(cos²(π/8)) = {ex_H2(np.cos(ex_t)**2):.4f} bit")
print(f"[연습 2] θ=π/4: S = {ex_entropy(ex_state(np.pi/4))[0]:.4f} bit,  상태 = |Φ+> ? {np.allclose(ex_state(np.pi/4), phi_plus)}")
print(f"[연습 2] θ=0: S = {ex_entropy(ex_state(0))[0]:.4f} (곱상태? {is_product(ex_state(0))}),  "
      f"θ=π/2: S = {ex_entropy(ex_state(np.pi/2))[0]:.4f} (곱상태? {is_product(ex_state(np.pi/2))})")
ex_grid = np.linspace(0, np.pi/2, 181)                                # 0.5° 간격
ex_Ss = [ex_entropy(ex_state(t))[0] for t in ex_grid]
print(f"[연습 2] 0 ≤ θ ≤ π/2 격자(181점)에서 S 최대: θ = {ex_grid[int(np.argmax(ex_Ss))]/np.pi:.4f} π, S = {max(ex_Ss):.4f} bit")
print("[연습 2] CNOT·(Ry(2θ)⊗I)|00> = cosθ|00> + sinθ|11>  (θ = π/8, 0.3, 1.1):",
      all(np.allclose(CNOT @ kron(Ry(2*t), I) @ ket('00'), ex_state(t)) for t in (np.pi/8, 0.3, 1.1)))
assert np.allclose(ex_rA, np.diag([np.cos(ex_t)**2, np.sin(ex_t)**2])), "8절 문제 2: ρ_A = diag(cos²θ, sin²θ)"
assert abs(np.cos(ex_t)**2 - 0.8536) < 5e-5, "8절 문제 2: cos²(π/8) = 0.8536"
assert abs(ex_S - 0.601) < 5e-4 and abs(ex_S - ex_H2(np.cos(ex_t)**2)) < 1e-9, "8절 문제 2: S = H₂(0.8536) ≈ 0.601 bit"
assert abs(ex_entropy(ex_state(np.pi/4))[0] - 1) < 1e-9 and np.allclose(ex_state(np.pi/4), phi_plus), "8절 문제 2: θ = π/4 에서 최대 1 bit, |Φ+>"
assert abs(ex_entropy(ex_state(0))[0]) < 1e-9 and abs(ex_entropy(ex_state(np.pi/2))[0]) < 1e-9, "8절 문제 2: θ = 0, π/2 에서 S = 0"
assert is_product(ex_state(0)) and is_product(ex_state(np.pi/2)), "8절 문제 2: θ = 0, π/2 는 곱상태"
assert abs(ex_grid[int(np.argmax(ex_Ss))] - np.pi/4) < 1e-9, "8절 문제 2: 격자에서 S 최대는 θ = π/4"
assert all(np.allclose(CNOT @ kron(Ry(2*t), I) @ ket('00'), ex_state(t)) for t in (np.pi/8, 0.3, 1.1)), "8절 문제 2: Ry(2θ) 뒤 CNOT 으로 만든다"

# 문제 3: CNOT|+>|-> = |->|->  (위상 되돌림)
plus = H @ ket0; minus = H @ ket1
# 주의: 이 파일의 kron() 은 행렬용이라 벡터를 넣으면 (1,4) 모양이 된다 → 벡터끼리는 np.kron 을 쓴다
print("[연습 3] X|-> = -|-> :", np.allclose(X @ minus, -minus))
ex_out = CNOT @ np.kron(plus, minus)
print("[연습 3] CNOT @ kron(plus, minus) = kron(minus, minus) :", np.allclose(ex_out, np.kron(minus, minus)))
print(f"[연습 3] 입력 |+>|-> 와 출력의 겹침 |<+-|출력>|² = {abs(np.vdot(np.kron(plus, minus), ex_out))**2:.4f}"
      " → 제어 큐비트가 |+> 에서 |-> 로 바뀌었다(표적 |-> 는 그대로)")
print("[연습 3] (H⊗H) CNOT (H⊗H) = 뒤집힌 CNOT :", np.allclose(kron(H, H) @ CNOT @ kron(H, H), SWAP @ CNOT @ SWAP))
assert np.allclose(X @ minus, -minus), "8절 문제 3: X|-> = -|->"
assert np.allclose(ex_out, np.kron(minus, minus)), "8절 문제 3: CNOT|+>|-> = |->|->"
assert abs(np.vdot(np.kron(plus, minus), ex_out))**2 < 1e-12, "8절 문제 3: 출력은 입력과 직교(제어 큐비트가 바뀌었다)"
