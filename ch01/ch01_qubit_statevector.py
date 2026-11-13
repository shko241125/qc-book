"""1장 큐비트와 상태벡터 — numpy만으로 큐비트 하나를 끝까지 계산한다.
실행: python ch01_qubit_statevector.py
환경: Python 3.13, numpy 2.5 (SDK 없이 동작)
"""
import os
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
assert normalized(plus), "5절: |+> 정규화 검사 True"
assert np.allclose(plus, [0.7071, 0.7071], atol=5e-5), "5절: |+> 상태벡터 (0.7071, 0.7071)"

# ── 2. Born 규칙: 측정 확률 = |진폭|^2 ─────────────────────────────
def probs(psi):
    return np.abs(psi) ** 2

def sample(psi, shots):
    p = probs(psi)
    outcomes = rng.choice([0, 1], size=shots, p=p)
    return {"0": int((outcomes == 0).sum()), "1": int((outcomes == 1).sum())}

print("2) P(0), P(1) for |+> =", probs(plus))
assert np.allclose(probs(plus), [0.5, 0.5]), "5절: |+> 의 P(0)=P(1)=0.5 (예측 512회)"
counts_plus = sample(plus, shots=1024)
print("   1024 shots on |+> ->", counts_plus)
assert counts_plus == {"0": 517, "1": 507}, "5절: |+> 1024샷 = 0이 517회, 1이 507회"

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
for psi, r_md in [(ket0, (0, 0, 1)), (ket1, (0, 0, -1)), (plus, (1, 0, 0)),
                  (state_from_angles(np.pi / 2, np.pi / 2), (0, 1, 0))]:
    assert np.allclose(bloch(psi), r_md), "5절: Bloch 좌표 표(|0>, |1>, |+>, |+i>)"

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
for r, pred_md, n1_md in zip(rows, (0.0000, 0.1464, 0.5000, 0.8536, 1.0000), (0, 290, 1054, 1693, 2000)):
    assert abs(r[2] - pred_md) < 5e-5, "5절: θ 표의 예측 P(1)"
    assert round(r[3] * 2000) == n1_md, "5절: θ 표의 측정 P(1) (0.0000, 0.1450, 0.5270, 0.8465, 1.0000 × 2000샷)"
for r, diff_md in zip(rows, (0.0000, -0.0014, +0.0270, -0.0071, 0.0000)):
    assert abs((r[3] - r[2]) - diff_md) < 5e-5, "5절: θ 표의 차이 열"
assert abs((rows[2][3] - rows[2][2]) / np.sqrt(0.25 / 2000) - 2.4) < 0.05, "5절·연습 2: π/2 의 +0.027 은 표준오차의 약 2.4배"

# ── 5. 반례: 큐비트는 확률 비트가 아니다 (H·H = I) ─────────────────
H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
psi_HH = H @ (H @ ket0)
counts_HH = sample(psi_HH, 1024)
print("5) H H |0> =", psi_HH, "-> 1024 shots:", counts_HH)
assert np.allclose(psi_HH, [1, 0]), "7절: HH|0> 의 상태벡터는 (1, 0)"
assert counts_HH == {"0": 1024, "1": 0}, "7절: HH|0> 1024샷 = 0이 1024회, 1이 0회"
print("   (고전 동전을 두 번 던지면 50:50이 유지되지만, 큐비트는 |0>로 되돌아온다)")

# ── 6. 전역 위상은 보이지 않는다 ───────────────────────────────────
phase = np.exp(1j * 0.9)
print("6) probs(|+>) =", probs(plus), " probs(e^{i0.9}|+>) =", probs(phase * plus))
print("   Bloch 좌표도 같다:", bloch(plus), bloch(phase * plus))
assert np.allclose(probs(plus), probs(phase * plus)), "7절: 전역 위상 e^{0.9i} 는 확률을 바꾸지 않는다"
assert np.allclose(bloch(plus), bloch(phase * plus)), "7절: 전역 위상 e^{0.9i} 는 Bloch 좌표를 바꾸지 않는다"

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
plt.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)), "ch01_results.png"), dpi=150)
print("7) saved ch01_results.png")
