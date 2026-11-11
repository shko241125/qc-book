"""1장 큐비트와 상태벡터 — numpy만으로 큐비트 하나를 끝까지 계산한다.
실행: python ch01_qubit_statevector.py
환경: Python 3.13, numpy 2.5 (SDK 없이 동작)
"""
import numpy as np

np.set_printoptions(precision=4, suppress=True)
rng = np.random.default_rng(seed=7)

# ── 1. 계산 기저와 상태벡터 ─────────────────────────────────────────
ket0 = np.array([1, 0], dtype=complex)          # |0>
ket1 = np.array([0, 1], dtype=complex)          # |1>

def normalized(psi):
    return np.isclose(np.vdot(psi, psi).real, 1.0)

# 중첩 상태 |+> = (|0> + |1>)/sqrt(2)
plus = (ket0 + ket1) / np.sqrt(2)
print("1) |+> =", plus, "| normalized:", normalized(plus))

# ── 2. Born 규칙: 측정 확률 = |진폭|^2 ─────────────────────────────
def probs(psi):
    return np.abs(psi) ** 2

def sample(psi, shots):
    p = probs(psi)
    outcomes = rng.choice([0, 1], size=shots, p=p)
    return {"0": int((outcomes == 0).sum()), "1": int((outcomes == 1).sum())}

print("2) P(0), P(1) for |+> =", probs(plus))
counts_plus = sample(plus, shots=1024)
print("   1024 shots on |+> ->", counts_plus)

# ── 3. Bloch 구 좌표 (<X>, <Y>, <Z>) ──────────────────────────────
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)

def bloch(psi):
    return np.array([np.vdot(psi, P @ psi).real for P in (X, Y, Z)])

def state_from_angles(theta, phi):
    """|psi> = cos(theta/2)|0> + e^{i phi} sin(theta/2)|1>"""
    return np.array([np.cos(theta / 2), np.exp(1j * phi) * np.sin(theta / 2)])

for name, psi in [("|0>", ket0), ("|1>", ket1), ("|+>", plus),
                  ("|+i>", state_from_angles(np.pi / 2, np.pi / 2))]:
    print(f"3) Bloch({name}) = {bloch(psi)}")

# ── 4. 각도 θ를 바꾸며 예측 vs 측정 ───────────────────────────────
print("4) theta sweep: P(1) predicted = sin^2(theta/2) vs measured (2000 shots)")
rows = []
for k in range(5):
    theta = k * np.pi / 4
    psi = state_from_angles(theta, 0.0)
    pred = float(np.sin(theta / 2) ** 2)
    meas = sample(psi, shots=2000)["1"] / 2000
    rows.append((k, theta, pred, meas))
    print(f"   theta={k}pi/4  predicted={pred:.4f}  measured={meas:.4f}")

# ── 5. 반례: 큐비트는 확률 비트가 아니다 (H·H = I) ─────────────────
H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
psi_HH = H @ (H @ ket0)
print("5) H H |0> =", psi_HH, "-> 1024 shots:", sample(psi_HH, 1024))
print("   (고전 동전을 두 번 던지면 50:50이 유지되지만, 큐비트는 |0>로 되돌아온다)")

# ── 6. 전역 위상은 보이지 않는다 ───────────────────────────────────
phase = np.exp(1j * 0.9)
print("6) probs(|+>) =", probs(plus), " probs(e^{i0.9}|+>) =", probs(phase * plus))
print("   Bloch 좌표도 같다:", bloch(plus), bloch(phase * plus))

# ── 7. 그림 저장 ───────────────────────────────────────────────────
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Noto Sans CJK JP"
plt.rcParams["axes.unicode_minus"] = False

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 3.4))
ax1.bar(list(counts_plus.keys()), list(counts_plus.values()), color="#4C6EF5")
ax1.set_title("H|0> 측정 1024샷")
ax1.set_xlabel("결과"); ax1.set_ylabel("횟수")
for i, (k, v) in enumerate(counts_plus.items()):
    ax1.text(i, v + 10, str(v), ha="center")

thetas = np.linspace(0, np.pi, 200)
ax2.plot(thetas, np.sin(thetas / 2) ** 2, color="#868E96", label="예측 sin²(θ/2)")
ax2.scatter([r[1] for r in rows], [r[3] for r in rows], color="#4C6EF5", zorder=3, label="측정 (2000샷)")
ax2.set_xticks([0, np.pi / 4, np.pi / 2, 3 * np.pi / 4, np.pi])
ax2.set_xticklabels(["0", "π/4", "π/2", "3π/4", "π"])
ax2.set_xlabel("θ"); ax2.set_ylabel("P(1)"); ax2.set_title("θ에 따른 P(1)")
ax2.legend(frameon=False)
plt.tight_layout()
plt.savefig("ch01_results.png", dpi=150)
print("7) saved ch01_results.png")
