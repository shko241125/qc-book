"""2장 단일 큐비트 게이트 — 행렬·회전·항등식을 numpy로 전부 검증한다.
실행: python ch02_single_qubit_gates.py
환경: Python 3.13, numpy 2.5 (SDK 없이 동작)
"""
import numpy as np

np.set_printoptions(precision=4, suppress=True)
rng = np.random.default_rng(seed=7)

# ── 1. 고정 게이트 ─────────────────────────────────────────────────
I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
S = np.array([[1, 0], [0, 1j]], dtype=complex)
T = np.array([[1, 0], [0, np.exp(1j * np.pi / 4)]], dtype=complex)

# ── 2. 회전 게이트와 일반 U ────────────────────────────────────────
def Rx(t): return np.array([[np.cos(t/2), -1j*np.sin(t/2)], [-1j*np.sin(t/2), np.cos(t/2)]])
def Ry(t): return np.array([[np.cos(t/2), -np.sin(t/2)], [np.sin(t/2), np.cos(t/2)]], dtype=complex)
def Rz(t): return np.array([[np.exp(-1j*t/2), 0], [0, np.exp(1j*t/2)]])
def P(lam): return np.array([[1, 0], [0, np.exp(1j*lam)]])          # phase gate
def U(theta, phi, lam):                                              # Qiskit UGate 정의
    return np.array([[np.cos(theta/2), -np.exp(1j*lam)*np.sin(theta/2)],
                     [np.exp(1j*phi)*np.sin(theta/2), np.exp(1j*(phi+lam))*np.cos(theta/2)]])

def is_unitary(M): return np.allclose(M.conj().T @ M, I)
def eq_up_to_phase(A, B):
    """A = e^{iγ} B 인지: 전역 위상만 다른 두 게이트는 같은 게이트"""
    idx = np.unravel_index(np.argmax(np.abs(B)), B.shape)
    gamma = A[idx] / B[idx]
    return np.allclose(A, gamma * B), np.angle(gamma)

ket0 = np.array([1, 0], dtype=complex)
def bloch(psi): return np.array([np.vdot(psi, O @ psi).real for O in (X, Y, Z)])

print("1) 유니터리 검사:", {n: is_unitary(M) for n, M in
      [("X", X), ("H", H), ("S", S), ("T", T), ("Rx(1.1)", Rx(1.1)), ("U(1,2,3)", U(1, 2, 3))]})

# ── 3. 항등식 ──────────────────────────────────────────────────────
print("2) 항등식")
print("   HZH = X :", np.allclose(H @ Z @ H, X))
print("   HXH = Z :", np.allclose(H @ X @ H, Z))
print("   T^2 = S :", np.allclose(T @ T, S), "  S^2 = Z :", np.allclose(S @ S, Z))
print("   X^2=Y^2=Z^2=H^2=I :", all(np.allclose(M @ M, I) for M in (X, Y, Z, H)))
print("   XYZ = iI :", np.allclose(X @ Y @ Z, 1j * I))
print("   Rz(pi) = -i Z (전역 위상만 다름):", eq_up_to_phase(Rz(np.pi), Z))
print("   P(lam) = e^{i lam/2} Rz(lam):", eq_up_to_phase(P(0.7), Rz(0.7)))
print("   H = U(pi/2, 0, pi):", np.allclose(H, U(np.pi/2, 0, np.pi)))
print("   X = U(pi, 0, pi):", np.allclose(X, U(np.pi, 0, np.pi)))
print("   Rz(t) = U(0,0,t) up to phase:", eq_up_to_phase(Rz(0.9), U(0, 0, 0.9)))

# ── 4. 같은 축 회전은 각이 더해진다 ───────────────────────────────
a, b = 0.4, 1.3
print("3) Rx(a)Rx(b) = Rx(a+b):", np.allclose(Rx(a) @ Rx(b), Rx(a + b)))
print("   Rx(a)Ry(b) = Ry(b)Rx(a)? :", np.allclose(Rx(a) @ Ry(b), Ry(b) @ Rx(a)), "(축이 다르면 교환 안 됨)")

# ── 5. ZYZ 분해: 임의 유니터리 = e^{iα} Rz(β) Ry(γ) Rz(δ) ───────────
def random_unitary():
    A = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
    Q, R = np.linalg.qr(A)
    return Q @ np.diag(np.diag(R) / np.abs(np.diag(R)))

def zyz(Uin):
    # 전역 위상 제거 → SU(2)
    det = np.linalg.det(Uin)
    V = Uin / np.sqrt(det)
    gamma = 2 * np.arctan2(np.abs(V[1, 0]), np.abs(V[0, 0]))
    if np.isclose(np.sin(gamma / 2), 0):
        beta, delta = 2 * np.angle(V[1, 1]), 0.0
    elif np.isclose(np.cos(gamma / 2), 0):
        beta, delta = 2 * np.angle(V[1, 0]), 0.0
    else:
        beta = np.angle(V[1, 1]) + np.angle(V[1, 0])
        delta = np.angle(V[1, 1]) - np.angle(V[1, 0])
    return beta, gamma, delta

Ur = random_unitary()
beta, gamma, delta = zyz(Ur)
rebuilt = Rz(beta) @ Ry(gamma) @ Rz(delta)
ok, alpha = eq_up_to_phase(Ur, rebuilt)
print(f"4) ZYZ 분해: beta={beta:.4f} gamma={gamma:.4f} delta={delta:.4f}  재구성 일치(위상 제외)={ok}, alpha={alpha:.4f}")

# ── 6. Bloch 궤적: |0> → H → T → H ────────────────────────────────
print("5) Bloch 궤적 |0> -H-> -T-> -H->")
psi = ket0.copy()
for name, G in [("start", I), ("H", H), ("T", T), ("H", H)]:
    psi = G @ psi
    print(f"   after {name:5s}: bloch={bloch(psi)}  P(1)={abs(psi[1])**2:.4f}")

# ── 7. 반례: Rz는 계산 기저 측정에 안 보인다 ─────────────────────
print("6) Rz 만 걸면 P(1)은 불변:")
for t in (0.3, 1.0, 2.5):
    print(f"   Rz({t}) H|0>: P(1)={abs((Rz(t) @ H @ ket0)[1])**2:.4f}  bloch={bloch(Rz(t) @ H @ ket0)}")
print("   그러나 H 를 한 번 더 걸면 위상이 확률로 바뀐다:")
for t in (0.0, np.pi/2, np.pi):
    print(f"   H Rz({t:.4f}) H|0>: P(1)={abs((H @ Rz(t) @ H @ ket0)[1])**2:.4f}")

# ── 8. 그림: Bloch 구 투영과 H-Rz(t)-H 곡선 ─────────────────────
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = ["Noto Sans CJK JP", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 3.6))
# 왼쪽: H -> Rz(t) -> H 에서 t 에 따른 P(1)
ts = np.linspace(0, 2*np.pi, 200)
p1 = [abs((H @ Rz(t) @ H @ ket0)[1])**2 for t in ts]
p1_rz_only = [abs((Rz(t) @ H @ ket0)[1])**2 for t in ts]
ax1.plot(ts, p1, color="#4C6EF5", label="H·Rz(t)·H |0⟩")
ax1.plot(ts, p1_rz_only, color="#868E96", ls="--", label="Rz(t)·H |0⟩")
ax1.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi]); ax1.set_xticklabels(["0", "π/2", "π", "3π/2", "2π"])
ax1.set_xlabel("t"); ax1.set_ylabel("P(1)"); ax1.set_title("Rz는 H로 감싸야 보인다"); ax1.legend(frameon=False, fontsize=9)
# 오른쪽: 적도 평면(x-y) 위 궤적 |+> 에 Rz(t)
xs = [bloch(Rz(t) @ H @ ket0)[0] for t in ts]; ys = [bloch(Rz(t) @ H @ ket0)[1] for t in ts]
ax2.plot(np.cos(ts), np.sin(ts), color="#dee2e6"); ax2.plot(xs, ys, color="#4C6EF5")
for t, lab in [(0, "|+⟩"), (np.pi/2, "|+i⟩"), (np.pi, "|−⟩")]:
    v = bloch(Rz(t) @ H @ ket0); ax2.scatter(v[0], v[1], color="#4C6EF5", zorder=3); ax2.annotate(lab, (v[0], v[1]), xytext=(6, -14 if lab == "|+i⟩" else 6), textcoords="offset points")
ax2.set_aspect("equal"); ax2.set_xlabel("⟨X⟩"); ax2.set_ylabel("⟨Y⟩"); ax2.set_title("Rz(t)는 적도 위를 돈다")
plt.tight_layout(); plt.savefig("ch02_results.png", dpi=150)
print("7) saved ch02_results.png")

# ── 연습문제 검산 (8절) ──────────────────────────────────────────
# 8절 해답이 인용하는 수치·행렬을 전부 다시 계산한다. 위쪽 코드·난수(rng)·그림은 건드리지 않는다.
print("── 연습문제 검산 (8절) ──")
SX = np.array([[1+1j, 1-1j], [1-1j, 1+1j]], dtype=complex) / 2       # √X (문제 2에 주어진 행렬)

# 문제 1: HTH = e^{iπ/8} Rx(π/4)
ex_ok, ex_ph = eq_up_to_phase(T, Rz(np.pi/4))
print(f"[연습 1] T = e^{{iγ}} Rz(π/4): 일치={ex_ok}, γ={ex_ph:.4f}  (π/8 = {np.pi/8:.4f})")
print("[연습 1] H Rz(θ) H = Rx(θ)  (θ = 0.3, 1.0, 2.5):", all(np.allclose(H @ Rz(t) @ H, Rx(t)) for t in (0.3, 1.0, 2.5)))
ex_ok, ex_ph = eq_up_to_phase(H @ T @ H, Rx(np.pi/4))
print(f"[연습 1] eq_up_to_phase(H @ T @ H, Rx(np.pi/4)) = ({ex_ok}, {ex_ph:.4f})  → HTH = e^{{iπ/8}} Rx(π/4), 회전각 π/4 = 45°")
print("[연습 1] HTH = e^{iπ/8} Rx(π/4) (위상까지 직접 비교):", np.allclose(H @ T @ H, np.exp(1j*np.pi/8) * Rx(np.pi/4)))
ex_dev = max(abs(abs((H @ Rz(t) @ H @ ket0)[1])**2 - np.sin(t/2)**2) for t in np.linspace(0, 2*np.pi, 200))
print(f"[연습 1] H·Rz(t)·H|0> 의 P(1) = sin²(t/2): 200점 최대 편차 = {ex_dev:.1e}")

# 문제 2: H = e^{iπ/4} Rz(π/2) √X Rz(π/2)
print("[연습 2] √X·√X = X :", np.allclose(SX @ SX, X))
ex_ok, ex_ph = eq_up_to_phase(SX, Rx(np.pi/2))
print(f"[연습 2] √X = e^{{iγ}} Rx(π/2) (x축 90°): 일치={ex_ok}, γ={ex_ph:.4f}  (π/4 = {np.pi/4:.4f})")
ex_ok, ex_ph = eq_up_to_phase(Rz(np.pi/2) @ SX @ Rz(np.pi/2), H)
print(f"[연습 2] eq_up_to_phase(Rz(np.pi/2) @ SX @ Rz(np.pi/2), H) = ({ex_ok}, {ex_ph:.4f})  (-π/4 = {-np.pi/4:.4f})")
print("[연습 2] H = e^{iπ/4} Rz(π/2) √X Rz(π/2) (위상까지 직접 비교):",
      np.allclose(H, np.exp(1j*np.pi/4) * Rz(np.pi/2) @ SX @ Rz(np.pi/2)))
ex_axis = (X + Z) / np.sqrt(2)                                        # 축 n = (x+z)/√2 의 n·σ
ex_rot = np.cos(np.pi/2) * I - 1j * np.sin(np.pi/2) * ex_axis         # exp(-i·π/2·n·σ) = 그 축 180° 회전
ex_ok, ex_ph = eq_up_to_phase(Rz(np.pi/2) @ SX @ Rz(np.pi/2), ex_rot)
print(f"[연습 2] Rz(π/2)·√X·Rz(π/2) = (x+z)/√2 축 180° 회전 (위상 제외): 일치={ex_ok}, 위상={ex_ph:.4f}")

# 문제 3: S H S† = (Y+Z)/√2
ex_G = S @ H @ S.conj().T
print("[연습 3] S X S† = Y (S 는 x축을 y축으로):", np.allclose(S @ X @ S.conj().T, Y))
print("[연습 3] np.allclose(S @ H @ S.conj().T, (Y+Z)/np.sqrt(2)) =", np.allclose(ex_G, (Y + Z) / np.sqrt(2)))
print("[연습 3] S H S† = (√2/2)[[1, -i], [i, -1]] :", np.allclose(ex_G, np.sqrt(2)/2 * np.array([[1, -1j], [1j, -1]])))
ex_out = ex_G @ ket0
print("[연습 3] S H S†|0> = (1, i)/√2 = |+i> :", np.allclose(ex_out, np.array([1, 1j]) / np.sqrt(2)),
      " Bloch 좌표 =", np.round(bloch(ex_out), 4) + 0.0)
print("[연습 3] 비교: H|0> 의 Bloch 좌표 =", np.round(bloch(H @ ket0), 4) + 0.0, "(북극 → x축)")
