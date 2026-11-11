"""19장 후보 문제 선정 — numpy 검증 스크립트
작업계획서의 후보 10개 각각에 대해 판정을 가르는 '결정적 숫자'를 작은 크기에서 직접 잰다.
 (1) 음성 피치 추정을 QFT로: 필요한 샷 수 × 상태 준비 게이트 수 vs 고전 FFT 한 번
 (2) 화자 임베딩 양자 커널: 커널 집중(큐비트 수에 따른 비대각 원소의 평균·표준편차), RBF 커널과의 분리 성능 비교
 (3) 분자 바닥상태 VQE: 큐비트 수와 FCI 차원(조합 계산)
 (4) Heisenberg 사슬 동역학: 자유 페르미온 여부, 얽힘 성장과 MPS 결합 차원, Trotter 오차와 잡음 생존 확률
 (5) HP 격자 단백질 접힘: 자기 회피 보행 전수 열거(정답·탐색 공간), 턴 인코딩 큐비트 수
 (7) Schwinger 모형: Majorana 차수(자유 페르미온 여부)
 (8) 희소 선형계: 조건수와 켤레 기울기법 반복 수의 크기 의존성
 (6)(9)(10)은 앞 장의 수치를 다시 쓰므로 계산하지 않는다.
규약: big-endian. 난수는 seed = 19 고정.
"""
import itertools, math, os, time
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import expm_multiply, cg
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
matplotlib.rcParams["font.family"] = ["Noto Sans CJK JP", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

T0 = time.time()
rng = np.random.default_rng(19)

# =====================================================================
# (1) 음성 피치 추정을 QFT/QPE로
# =====================================================================
print("=" * 78); print("(1) 음성 피치 추정을 QFT로 — 적재 비용과 샷 수"); print("=" * 78)
fs, N = 16000, 1024                                   # 16 kHz, 64 ms 프레임, n = 10 큐비트
f0 = 180.0
tt = np.arange(N) / fs
clean = sum(np.sin(2 * np.pi * k * f0 * tt + 0.3 * k) / k for k in range(1, 11))
noise = rng.normal(size=N); noise *= np.linalg.norm(clean) / np.linalg.norm(noise) / 10 ** (10 / 20)   # SNR 10 dB
x = (clean + noise) * np.hanning(N)
X = np.fft.fft(x); P = np.abs(X) ** 2; P /= P.sum()   # 진폭 인코딩 뒤 QFT 를 재면 구간 k 가 P_k 확률로 나온다
k_true = int(round(f0 * N / fs)); k_fft = int(np.argmax(P[: N // 2]))
t_fft = min((lambda s: (np.fft.fft(x), time.perf_counter() - s)[1])(time.perf_counter()) for _ in range(200))
print(f"프레임 N = {N} (n = 10 큐비트), f0 = {f0} Hz → 구간 {f0*N/fs:.2f}; FFT 최대 구간 {k_fft} ({k_fft*fs/N:.1f} Hz), numpy FFT 한 번 {t_fft*1e6:.1f} µs")
p_fund = P[k_true - 1:k_true + 2].sum() + P[N - k_true - 1:N - k_true + 2].sum()
print(f"기본 주파수 구간(±1)에 실린 확률 {p_fund:.3f}  (나머지는 배음과 잡음)")
N_pit, P_pit, k_pit, fs_pit = N, P.copy(), k_true, fs                      # (R) 절에서 다시 쓴다
shots_list = [1, 2, 4, 8, 16, 32, 64, 128]; pitch_rows = []
for s in shots_list:
    ok = 0
    for _ in range(2000):
        ks = rng.choice(N, size=s, p=P); ks = np.minimum(ks, N - ks)            # 음의 주파수 접기
        vals, cnts = np.unique(ks, return_counts=True); est = vals[np.argmax(cnts)]
        ok += abs(est - k_true) <= 1
    pitch_rows.append((s, ok / 2000))
print("샷 수별 '최빈 구간 = 기본 주파수(±1 구간)' 성공률 (2000회): " + ", ".join(f"{s}: {p:.3f}" for s, p in pitch_rows))
s99 = next((s for s, p in pitch_rows if p >= 0.99), None)
prep = 2 ** 10 - 2                                     # 임의 실수 진폭 상태 준비의 CNOT 수(균일 제어 회전 분해), 큐비트 10개
print(f"성공률 99%에 필요한 샷 수 ≈ {s99};  샷마다 상태를 다시 준비해야 하므로 CNOT ≈ {s99} × {prep} = {s99*prep:,} (+ QFT {10*9//2}개의 제어 위상)")
print(f"고전: FFT 한 번(덧셈·곱셈 약 N·log₂N = {N*10:,}회, {t_fft*1e6:.1f} µs), 표본 {N}개를 읽는 것 자체가 O(N)")

# =====================================================================
# (2) 화자 임베딩 양자 커널
# =====================================================================
print("\n" + "=" * 78); print("(2) 화자 임베딩 양자 커널 — 커널 집중과 RBF 대비 분리 성능"); print("=" * 78)
K_spk, per, D = 4, 25, 192
cent = rng.normal(size=(K_spk, D)); cent /= np.linalg.norm(cent, axis=1, keepdims=True)
emb = np.vstack([c + 0.085 * rng.normal(size=(per, D)) for c in cent]); emb /= np.linalg.norm(emb, axis=1, keepdims=True)
lab = np.repeat(np.arange(K_spk), per)
cosm = emb @ emb.T; same = lab[:, None] == lab[None, :]; off = ~np.eye(len(lab), dtype=bool)
print(f"합성 임베딩: 화자 {K_spk}명 × 발화 {per}개, {D}차원; 코사인 유사도 평균 — 같은 화자 {cosm[same & off].mean():.3f}, 다른 화자 {cosm[~same].mean():.3f}")


def hadamard_all(v, n):
    """H^{⊗n} 을 벡터(배치 × 2^n)에 적용 (빠른 Walsh–Hadamard 변환)."""
    v = v.reshape(v.shape[0], *([2] * n))
    for ax in range(1, n + 1):
        a, b = np.take(v, 0, axis=ax), np.take(v, 1, axis=ax)
        v = np.stack([a + b, a - b], axis=ax) / np.sqrt(2)
    return v.reshape(v.shape[0], -1)


def zz_feature_states(F, reps=2):
    """Qiskit ZZFeatureMap 정의(6장): 반복마다 H^{⊗n}, P(2x_i), 쌍마다 CX·P(2(π−x_i)(π−x_j))·CX.
    대각 위상 = exp(i[Σ_i 2x_i b_i + Σ_{i<j} 2(π−x_i)(π−x_j)(b_i⊕b_j)])."""
    m, n = F.shape
    bits = (np.arange(2 ** n)[:, None] >> np.arange(n - 1, -1, -1)) & 1          # big-endian
    phase = 2 * F @ bits.T
    for i in range(n):
        for j in range(i + 1, n):
            phase += 2 * np.outer((np.pi - F[:, i]) * (np.pi - F[:, j]), bits[:, i] ^ bits[:, j])
    Dg = np.exp(1j * phase)
    psi = np.zeros((m, 2 ** n), dtype=complex); psi[:, 0] = 1
    for _ in range(reps): psi = Dg * hadamard_all(psi, n)
    return psi


def loo_1nn(Kmat, lab):
    """자기 자신을 뺀 최근접 이웃(커널 값 최대)의 화자가 같은 비율 — 군집이 갈라지는 정도의 간단한 지표."""
    Km = Kmat.copy(); np.fill_diagonal(Km, -np.inf)
    return float((lab[np.argmax(Km, axis=1)] == lab).mean())


def alignment(Kmat, lab):
    """중심화 커널–목표 정렬 A(K, YYᵀ) ∈ [−1, 1]."""
    m = len(lab); Hc = np.eye(m) - 1 / m; Y = (lab[:, None] == lab[None, :]).astype(float)
    Kc, Yc = Hc @ Kmat @ Hc, Hc @ Y @ Hc
    return float((Kc * Yc).sum() / np.linalg.norm(Kc) / np.linalg.norm(Yc))


embc = emb - emb.mean(0); U_, S_, Vt_ = np.linalg.svd(embc, full_matrices=False)
print("    n   설명 분산   | ZZ 커널: 비대각 평균  표준편차   1-NN   정렬  | RBF 커널: 1-NN   정렬 | 차이 σ_K 를 2σ 로 구분하는 샷 수 8K(1−K)/σ_K²")
kern_rows = []
for n in (2, 3, 4, 6, 8, 10):
    F = embc @ Vt_[:n].T                                                        # PCA n 성분
    Fs = np.pi * (F - F.min(0)) / (F.max(0) - F.min(0))                         # 특성을 [0, π] 로
    psi = zz_feature_states(Fs); Kq = np.abs(psi @ psi.conj().T) ** 2
    d2 = ((Fs[:, None, :] - Fs[None, :, :]) ** 2).sum(-1); Kr = np.exp(-d2 / np.median(d2[off]))
    mq, sq_ = Kq[off].mean(), Kq[off].std()
    shots = 8 * mq * (1 - mq) / sq_ ** 2                                        # 3절 식 (2), Δ = σ_K
    kern_rows.append((n, (S_[:n] ** 2).sum() / (S_ ** 2).sum(), mq, sq_, loo_1nn(Kq, lab), alignment(Kq, lab), loo_1nn(Kr, lab), alignment(Kr, lab), shots))
    r = kern_rows[-1]
    print(f"   {n:2d}    {r[1]:.3f}      |      {r[2]:.4f}       {r[3]:.4f}    {r[4]:.2f}   {r[5]:+.3f} |       {r[6]:.2f}   {r[7]:+.3f} |   {r[8]:9.0f}")
print("   대역폭 조절(n = 8): 특성을 [0, c·π] 로 줄이면")
bw_rows = []
F8 = embc @ Vt_[:8].T; F8n = (F8 - F8.min(0)) / (F8.max(0) - F8.min(0))
G_RBF = np.concatenate([np.geomspace(0.01, 100, 60), np.geomspace(100, 1e4, 30)[1:]])      # 폭 격자(위쪽을 넓힘)
for c_bw in (1.0, 0.5, 0.25, 0.1, 0.05, 0.02):
    psi = zz_feature_states(np.pi * c_bw * F8n); Kq = np.abs(psi @ psi.conj().T) ** 2
    d2 = ((F8n[:, None, :] - F8n[None, :, :]) ** 2).sum(-1)
    # 같은 데이터에서 양자 커널을 가장 잘 흉내 내는 RBF 폭을 격자 탐색 → 두 커널의 상관
    best_corr = max(np.corrcoef(Kq[off], np.exp(-d2[off] / g))[0, 1] for g in G_RBF)
    bw_rows.append((c_bw, Kq[off].mean(), Kq[off].std(), loo_1nn(Kq, lab), alignment(Kq, lab), best_corr))
    print(f"     c = {c_bw:4.2f}: 비대각 평균 {bw_rows[-1][1]:.4f}, 표준편차 {bw_rows[-1][2]:.4f}, 1-NN {bw_rows[-1][3]:.2f}, 정렬 {bw_rows[-1][4]:+.3f}, 가장 닮은 RBF 커널과의 상관 {best_corr:.3f}")
print(f"   참고: 192차원 코사인 유사도의 1-NN = {loo_1nn(cosm, lab):.2f}, 정렬 = {alignment(cosm, lab):+.3f};  무작위 상태 쌍의 평균 충실도 1/2ⁿ: n=10 → {1/2**10:.4f}")

# =====================================================================
# (3) 분자 바닥상태 VQE — 큐비트 수와 FCI 차원
# =====================================================================
print("\n" + "=" * 78); print("(3) 분자 바닥상태 — 큐비트 수(스핀 궤도 수)와 FCI 차원 C(M, N↑)·C(M, N↓)"); print("=" * 78)
mols = [("H₂ / STO-3G", 2, 2), ("LiH / STO-3G", 6, 4), ("H₂O / STO-3G", 7, 10), ("N₂ / cc-pVDZ", 28, 14), ("FeMoco 활성 공간 (54e, 54o)", 54, 54)]
fci_rows = []
print("    분자 / 기저                     공간 궤도 M   전자 수   큐비트 2M    FCI 차원")
for name, M, ne in mols:
    dim = math.comb(M, ne // 2) ** 2; fci_rows.append((name, M, ne, 2 * M, dim))
    print(f"    {name:28s}   {M:3d}        {ne:3d}       {2*M:4d}      {dim:.3e}" + ("  ← 정확 대각화가 순식간" if dim < 1e6 else ""))

# =====================================================================
# (4) Heisenberg 사슬 실시간 동역학
# =====================================================================
print("\n" + "=" * 78); print("(4) Heisenberg 사슬 동역학 — 얽힘 성장, Trotter 오차, 잡음 생존 확률"); print("=" * 78)
PAULI = {"I": sp.identity(2, format="csr", dtype=complex), "X": sp.csr_matrix([[0, 1], [1, 0]], dtype=complex),
         "Y": sp.csr_matrix([[0, -1j], [1j, 0]]), "Z": sp.csr_matrix([[1, 0], [0, -1]], dtype=complex)}


def pauli(n, spec):
    out = sp.identity(1, format="csr", dtype=complex)
    for q in range(n): out = sp.kron(out, PAULI[spec.get(q, "I")], format="csr")
    return out


def heisenberg(n, J=1.0):
    H = sp.csr_matrix((2 ** n, 2 ** n), dtype=complex)
    for i in range(n - 1):
        for a in "XYZ": H = H + J * pauli(n, {i: a, i + 1: a})
    return H


def majorana_degree(spec, n, axis):
    """18장 5절 (D)와 같은 함수: Jordan–Wigner 로 옮긴 파울리 문자열의 Majorana 차수."""
    deg, odd_right = 0, 0
    for j in range(n - 1, -1, -1):
        pj = spec.get(j, "I"); is_axis, is_id = (pj == axis), (pj == "I")
        if odd_right: is_axis, is_id = is_id, is_axis
        deg += 0 if is_id else (2 if is_axis else 1)
        if pj not in ("I", axis): odd_right ^= 1
    return deg


def schmidt(psi, n, cut): return np.linalg.svd(psi.reshape(2 ** cut, -1), compute_uv=False) ** 2
def entropy_bits(p): p = p[p > 1e-15]; return float(-(p * np.log2(p)).sum())
def chi_needed(p, infid=1e-6): return int(np.searchsorted(np.cumsum(np.sort(p)[::-1]), 1 - infid) + 1)


hterms = [{i: a, i + 1: a} for i in range(5) for a in "XYZ"]
print("Majorana 최대 차수(끈 축 X/Y/Z): " + ", ".join(f"{ax}: {max(majorana_degree(s, 6, ax) for s in hterms)}" for ax in "XYZ") + "  → 2차가 아님(자유 페르미온 아님)")
n = 16; H16 = heisenberg(n).tocsc()
neel = np.zeros(2 ** n, dtype=complex); neel[int("01" * (n // 2), 2)] = 1
ts = np.arange(0, 3.01, 0.25); heis_rows = []
for t, ph in zip(ts, expm_multiply(-1j * H16, neel, start=0, stop=ts[-1], num=len(ts), endpoint=True)):
    p = schmidt(ph, n, n // 2); heis_rows.append((t, entropy_bits(p), chi_needed(p)))
print(f"Néel 상태에서 급랭, n = {n}, H = Σ(XX+YY+ZZ):")
print("   " + ", ".join(f"t={t:.2f}: S={abs(S):.2f} χ={c}" for t, S, c in heis_rows if abs(t * 2 - round(t * 2)) < 1e-9))
grow = [r for r in heis_rows if 0.25 <= r[0] <= 1.25]
rate_h = np.polyfit([r[0] for r in grow], [r[1] for r in grow], 1)[0]
S_sat = np.mean([r[1] for r in heis_rows if r[0] >= 2.0])
print(f"t = 0.25~1.25 구간 S 증가 속도 ≈ {rate_h:.2f} 비트/단위시간;  t ≥ 2 에서 S ≈ {S_sat:.2f} 비트로 포화(반쪽 8큐비트의 {S_sat/8:.0%}), χ ≈ {heis_rows[-1][2]}")
# Trotter 오차 (n = 10, t = 2): 결합마다 exp(−i dt (XX+YY+ZZ)) 는 CNOT 3개(18장: 일반 2큐비트 게이트)
n_t = 10; Ht = heisenberg(n_t).tocsc(); neel_t = np.zeros(2 ** n_t, dtype=complex); neel_t[int("01" * (n_t // 2), 2)] = 1
exact = expm_multiply(-2.0j * Ht, neel_t)
bond = [sum(pauli(n_t, {i: a, i + 1: a}) for a in "XYZ").tocsc() for i in range(n_t - 1)]
even, odd = sum(bond[0::2]), sum(bond[1::2])
trot_rows = []
for r in (5, 10, 20, 40):
    dt = 2.0 / r; psi = neel_t.copy()
    for _ in range(r):                                                           # 2차 Trotter: 짝/2 · 홀 · 짝/2
        psi = expm_multiply(-0.5j * dt * even, psi); psi = expm_multiply(-1j * dt * odd, psi); psi = expm_multiply(-0.5j * dt * even, psi)
    infid = 1 - abs(np.vdot(exact, psi)) ** 2
    ncx = 3 * ((n_t // 2) * (r + 1) + (n_t // 2 - 1) * r)                        # 이웃한 짝 층은 합쳐진다: 짝 층 r+1 번, 홀 층 r 번
    trot_rows.append((r, infid, ncx, (1 - 3e-3) ** ncx))
    print(f"   n = {n_t}, t = 2, 2차 Trotter {r:2d}단계: 상태 비충실도 {infid:.2e}, CNOT {ncx:4d}개, 2큐비트 오류 0.3%일 때 무오류 확률 {(1-3e-3)**ncx:.3f}")
S_target = 13.3                                        # χ ≈ 10⁴ = 2^13.3 을 "MPS 가 버거워지는" 기준으로 둔다
t_hard = S_target / rate_h; n_hard = 2 * math.ceil(S_target / (S_sat / 8))
print(f"MPS 기준 χ ≈ 10⁴ (S ≈ {S_target} 비트): 포화 비율 {S_sat/8:.0%} 를 그대로 쓰면 반쪽 길이 ≥ {n_hard//2} → n ≳ {n_hard}, 시간 t ≈ {t_hard:.1f}")
r_hard = int(math.ceil(t_hard / 0.1))
ncx_hard = 3 * ((n_hard // 2) * (r_hard + 1) + (n_hard // 2 - 1) * r_hard)
print(f"   그 크기(n = {n_hard}, dt = 0.1, {r_hard}단계)의 CNOT ≈ {ncx_hard:,}개 → 무오류 확률 (1−0.003)^{ncx_hard} = {(1-3e-3)**ncx_hard:.1e}")

# =====================================================================
# (5) 단백질 접힘 HP 격자 모형
# =====================================================================
print("\n" + "=" * 78); print("(5) HP 격자 모형 — 자기 회피 보행 전수 열거"); print("=" * 78)
MOVES = [(1, 0), (0, 1), (-1, 0), (0, -1)]


def enumerate_hp(seq):
    """2D 정사각 격자의 자기 회피 보행을 모두 돌며 HP 에너지(이웃하지만 사슬에서 잇닿지 않은 H–H 쌍마다 −1)를 센다.
    반환: 보행 수, 최저 에너지, 최저 에너지 보행 수(회전·반사 포함)."""
    L = len(seq); count = 0; best = 0; nbest = 0
    pos = [(0, 0)]; occ = {(0, 0): 0}
    def energy():
        e = 0
        for i, (x_, y_) in enumerate(pos):
            if seq[i] != "H": continue
            for dx, dy in ((1, 0), (0, 1)):
                j = occ.get((x_ + dx, y_ + dy))
                if j is not None and seq[j] == "H" and abs(i - j) > 1: e -= 1
        return e
    def rec(k):
        nonlocal count, best, nbest
        if k == L:
            count += 1; e = energy()
            if e < best: best, nbest = e, 1
            elif e == best: nbest += 1
            return
        x_, y_ = pos[-1]
        for dx, dy in MOVES:
            q = (x_ + dx, y_ + dy)
            if q in occ: continue
            pos.append(q); occ[q] = k; rec(k + 1); pos.pop(); del occ[q]
    rec(1)
    return count, best, nbest


SAW_KNOWN = [4, 12, 36, 100, 284, 780, 2172, 5916, 16268, 44100, 120292]        # OEIS A001411, 걸음 수 1~11
seq_full = "HPHPPHHPHPPH"; hp_rows = []
print("    잔기 수 N   서열            보행 수 c(N−1)   알려진 값   최저 에너지   최저 보행 수   턴 인코딩 비트 2(N−1)   4^(N−1)        열거 시간")
for L in range(4, 13):
    seq = seq_full[:L]; t1 = time.time(); cnt, e0, n0 = enumerate_hp(seq); dt_ = time.time() - t1
    hp_rows.append((L, cnt, e0, n0, 2 * (L - 1), 4 ** (L - 1), dt_))
    print(f"      {L:2d}       {seq:12s}    {cnt:8d}       {SAW_KNOWN[L-2]:8d}       {e0:3d}          {n0:5d}            {2*(L-1):2d}              {4**(L-1):.2e}      {dt_:.2f} s")
assert all(r[1] == SAW_KNOWN[r[0] - 2] for r in hp_rows)
growth = (hp_rows[-1][1] / hp_rows[-3][1]) ** 0.5
print(f"보행 수의 걸음당 증가율 ≈ {growth:.3f} (알려진 연결 상수 2.638);  비트열 4^(N−1) 가운데 유효한 접힘의 비율: N = 12 에서 {hp_rows[-1][1]/hp_rows[-1][5]:.3%}")
t12 = hp_rows[-1][6]; t30 = t12 * 2.638 ** 18
print(f"이 코드의 속도로 N = 30 을 전수 열거하면 약 {t12:.2f} s × 2.638¹⁸ ≈ {t30:.1e} s ≈ {t30/86400:.0f}일 (보행 수만 센 외삽)")

# =====================================================================
# (7) Schwinger 모형 — 자유 페르미온인가
# =====================================================================
print("\n" + "=" * 78); print("(7) 격자 Schwinger 모형 (게이지장 소거 뒤) — Majorana 차수"); print("=" * 78)
ns = 6
sch = [{i: a, i + 1: a} for i in range(ns - 1) for a in "XY"] + [{i: "Z"} for i in range(ns)] \
      + [{k: "Z", l: "Z"} for k in range(ns) for l in range(k + 1, ns - 1)]       # 전기장 에너지 Σ_n L_n² 에서 나오는 장거리 ZZ
print(f"항: 도약 XX+YY {2*(ns-1)}개, 질량·전기장 Z {ns}개, 장거리 ZZ {sum(1 for s in sch if list(s.values()) == ['Z', 'Z'])}개")
print("Majorana 최대 차수: " + ", ".join(f"축 {ax}: {max(majorana_degree(s, ns, ax) for s in sch)}" for ax in "XYZ") + "  → 최선(축 Z)에서도 4차: 상호작용 모형")

# =====================================================================
# (8) 희소 선형계 — 조건수와 켤레 기울기법
# =====================================================================
print("\n" + "=" * 78); print("(8) 희소 선형계 — 조건수 κ 와 켤레 기울기법(CG) 반복 수의 크기 의존성"); print("=" * 78)


def cg_iters(A, b, tol=1e-8):
    it = [0]
    def cb(_): it[0] += 1
    xsol, info = cg(A, b, rtol=tol, maxiter=200000, callback=cb)
    return it[0], np.linalg.norm(A @ xsol - b) / np.linalg.norm(b)


def poisson1d(N): return sp.diags([-1.0, 2.0, -1.0], [-1, 0, 1], shape=(N, N), format="csr")
def poisson2d(L): I = sp.identity(L, format="csr"); T = poisson1d(L); return (sp.kron(I, T) + sp.kron(T, I)).tocsr()
def well_conditioned(N, rng):
    """대각 우세 희소 행렬 I + 0.2·(무작위 대칭 희소, 행마다 성분 약 4개): κ 가 N 과 무관하게 유지된다."""
    R = sp.random(N, N, density=2.0 / N, random_state=rng, format="csr", data_rvs=lambda k: rng.choice([-1.0, 1.0], size=k))
    S = (R + R.T) * 0.5
    return (sp.identity(N, format="csr") + 0.2 * S).tocsr()


lin_rows = []
print("    행렬 족                N        nnz/N    κ            CG 반복(1e-8)   √κ")
for N in (256, 1024, 4096):
    L = int(round(np.sqrt(N)))
    for name, A in (("1D Poisson", poisson1d(N)), ("2D Poisson", poisson2d(L)), ("잘 조건화된 희소", well_conditioned(N, rng))):
        if name == "1D Poisson":
            kap = (np.sin(N * np.pi / (2 * (N + 1))) / np.sin(np.pi / (2 * (N + 1)))) ** 2
        elif name == "2D Poisson":
            kap = (np.sin(L * np.pi / (2 * (L + 1))) / np.sin(np.pi / (2 * (L + 1)))) ** 2
        else:
            ev = np.linalg.eigvalsh(A.toarray()); kap = ev[-1] / ev[0]
        b = rng.normal(size=A.shape[0]); it, res = cg_iters(A, b)
        lin_rows.append((name, A.shape[0], A.nnz / A.shape[0], kap, it))
        print(f"    {name:14s}    {A.shape[0]:6d}      {A.nnz/A.shape[0]:.1f}    {kap:10.3e}     {it:6d}        {np.sqrt(kap):8.1f}")
print("    읽는 법: 고전 CG 의 비용 ≈ nnz·(반복 수). 양자 선형계의 호출 수는 κ 에 선형(최적 알고리즘 κ·log(1/ε)), 크기에는 다항로그.")
for name in ("1D Poisson", "2D Poisson", "잘 조건화된 희소"):
    rows = [r for r in lin_rows if r[0] == name]
    sl_k = np.polyfit(np.log([r[1] for r in rows]), np.log([r[3] for r in rows]), 1)[0]
    sl_i = np.polyfit(np.log([r[1] for r in rows]), np.log([max(r[4], 1) for r in rows]), 1)[0]
    print(f"    {name:14s}: κ ∝ N^{sl_k:.2f},  CG 반복 ∝ N^{sl_i:.2f}  → 고전 비용 ∝ N^{1+sl_i:.2f},  양자 호출 ∝ N^{sl_k:.2f}·polylog N")

# =====================================================================
# (R) 검토 보강 — 본문 주장을 받치는 추가 확인 (자체 난수 seed = 1920)
# =====================================================================
print("\n" + "=" * 78); print("(R) 검토 보강"); print("=" * 78)
rngR = np.random.default_rng(1920)
# R1. 피치: 2의 거듭제곱 사이의 샷 수, 그리고 14장의 CZ 시간(68 ns)으로 차례로 돌릴 때의 시간
fine = {}
for s in (20, 24, 28):
    ok = 0
    for _ in range(20000):
        ks = rngR.choice(N_pit, size=s, p=P_pit); ks = np.minimum(ks, N_pit - ks); v, c = np.unique(ks, return_counts=True)
        ok += abs(v[np.argmax(c)] - k_pit) <= 1
    fine[s] = ok / 20000
s99f = next((s for s in (20, 24, 28) if fine[s] >= 0.99), 32)
print("R1) 피치 성공률(2만 회): " + ", ".join(f"{s}샷 {v:.3f}" for s, v in fine.items()) + f" → 99%를 넘는 첫 값 {s99f}샷, CNOT ≈ {s99f} × {prep} = {s99f*prep:,}")
print(f"    ±1 구간 = ±{fs_pit/N_pit:.1f} Hz;  32샷 회로를 CZ 68 ns 로 차례로 돌리면 {32*prep*68e-9*1e3:.1f} ms = FFT({t_fft*1e6:.0f} µs)의 약 {32*prep*68e-9/t_fft:.0f}배")
# R2. 커널: 주성분마다 화자 정보가 얼마나 있는가, 화자 수를 늘린 대조 실험
Fall = embc @ Vt_[:10].T
bc = [np.var([Fall[lab == k, j].mean() for k in range(K_spk)]) / Fall[:, j].var() for j in range(10)]
print("R2) 주성분별 화자 간 분산 비율(화자 4명): " + ", ".join(f"{b:.3f}" for b in bc))
for Kc_ in (8, 11):
    cent_c = rngR.normal(size=(Kc_, D)); cent_c /= np.linalg.norm(cent_c, axis=1, keepdims=True)
    emb_c = np.vstack([c + 0.085 * rngR.normal(size=(per, D)) for c in cent_c]); emb_c /= np.linalg.norm(emb_c, axis=1, keepdims=True)
    lab_c = np.repeat(np.arange(Kc_), per); off_c = ~np.eye(len(lab_c), dtype=bool)
    ec = emb_c - emb_c.mean(0); Vc = np.linalg.svd(ec, full_matrices=False)[2]
    Fc = ec @ Vc[:10].T; Fcs = np.pi * (Fc - Fc.min(0)) / (Fc.max(0) - Fc.min(0))
    psi_c = zz_feature_states(Fcs); Kq_c = np.abs(psi_c @ psi_c.conj().T) ** 2
    d2c = ((Fcs[:, None, :] - Fcs[None, :, :]) ** 2).sum(-1); Kr_c = np.exp(-d2c / np.median(d2c[off_c]))
    chance = (per - 1) / (Kc_ * per - 1)
    print(f"    대조: 화자 {Kc_}명, n = 10 → ZZ 커널 비대각 평균 {Kq_c[off_c].mean():.4f}, 표준편차 {Kq_c[off_c].std():.4f}, 1-NN {loo_1nn(Kq_c, lab_c):.2f} (찍기 {chance:.2f}), 정렬 {alignment(Kq_c, lab_c):+.3f};  RBF 1-NN {loo_1nn(Kr_c, lab_c):.2f}, 정렬 {alignment(Kr_c, lab_c):+.3f}")
# R3. 분자: N₂/cc-pVDZ 에서 1s 궤도 둘을 얼리면
print(f"R3) N₂/cc-pVDZ, 1s 동결(궤도 26, 전자 10): FCI 차원 C(26,5)² = {math.comb(26, 5)**2:.2e}")
# R4. Heisenberg: χ 와 2^S 의 차이, χ 를 직접 외삽한 경계, 오류율을 낮췄을 때, 적분 가능성을 깬 J1–J2 사슬
print("R4) χ / 2^S: " + ", ".join(f"t={t:.2f}: {c/2**S:.1f}" for t, S, c in heis_rows if 0.25 <= t <= 2.0 and abs(t * 4 - round(t * 4)) < 1e-9 and abs(t * 2 - round(t * 2)) < 1e-9))
fitx = [r for r in heis_rows if 0.25 <= r[0] <= 1.5]
b_chi, a_chi = np.polyfit([r[0] for r in fitx], [math.log2(r[2]) for r in fitx], 1)
bS, aS = np.polyfit([r[0] for r in grow], [r[1] for r in grow], 1)
t_hard2 = (math.log2(1e4) - a_chi) / b_chi; S_hard2 = aS + bS * t_hard2; n_hard2 = 2 * math.ceil(S_hard2 / (S_sat / 8))
r_hard2 = int(math.ceil(t_hard2 / 0.1)); ncx_hard2 = 3 * ((n_hard2 // 2) * (r_hard2 + 1) + (n_hard2 // 2 - 1) * r_hard2)
print(f"    표의 χ 를 외삽(t = 0.25~1.5): log₂χ ≈ {a_chi:.2f} + {b_chi:.2f}·t → χ = 10⁴ 은 t ≈ {t_hard2:.1f}, 그때 S ≈ {S_hard2:.1f} 비트 → n ≳ {n_hard2}")
print(f"    그 크기(n = {n_hard2}, dt = 0.1, {r_hard2}단계)의 CNOT ≈ {ncx_hard2:,}개 → 무오류 확률 {(1-3e-3)**ncx_hard2:.1e}")
print(f"    n = {n_hard} 경계 회로: rG = 0.003 × {ncx_hard:,} = {3e-3*ncx_hard:.1f};  오류율이 3×10⁻⁴ 이면 rG = {3e-4*ncx_hard:.1f}, 무오류 확률 {(1-3e-4)**ncx_hard:.2f}")
H12 = H16.copy()
nH = 16; tsH = np.arange(0, 3.01, 0.25)
for i in range(nH - 2):
    for a in "XYZ": H12 = H12 + 0.5 * pauli(nH, {i: a, i + 2: a})
j2_rows = []
for t, ph in zip(tsH, expm_multiply(-1j * H12.tocsc(), neel, start=0, stop=tsH[-1], num=len(tsH), endpoint=True)):
    p = schmidt(ph, nH, nH // 2); j2_rows.append((t, entropy_bits(p), chi_needed(p)))
g2 = [r for r in j2_rows if 0.25 <= r[0] <= 1.25]
print(f"    J₁–J₂ 사슬(J₂ = 0.5, 비적분), 같은 Néel 급랭 n = {nH}: " + ", ".join(f"t={t:.1f}: S={S:.2f} χ={c}" for t, S, c in j2_rows if t in (1.0, 2.0, 3.0)) + f";  S 증가 속도 {np.polyfit([r[0] for r in g2], [r[1] for r in g2], 1)[0]:.2f} 비트/단위시간")
# R5. HP: 균등 무작위 비트열의 TTS₉₉ (턴 인코딩 비트열 4^(N−1) 개 중 최저 접힘의 비율)
print("R5) 무작위 비트열 추출의 TTS₉₉: " + ", ".join(f"N={L}: {math.log(0.01)/math.log(1 - n0/tot):,.0f}" for L, cnt, e0, n0, bits_, tot, dt_ in hp_rows if L >= 6))
# R6. Schwinger: 얽힘이 자라는가 — H = Σ(XX+YY) + m Σ(−1)ⁿ Z_n + J Σ_n L_n²,  L_n = Σ_{k≤n} (Z_k + (−1)^k)/2, 맨 진공 |1010…⟩ 에서 급랭
nS = 12
hopS = sum(pauli(nS, {i: a, i + 1: a}) for i in range(nS - 1) for a in "XY")
bitsS = ((np.arange(2 ** nS)[:, None] >> np.arange(nS - 1, -1, -1)) & 1); zS = 1 - 2 * bitsS; sgn = (-1.0) ** np.arange(nS)
Lf = np.cumsum((zS + sgn) / 2, axis=1)[:, :-1]
vac = np.zeros(2 ** nS, dtype=complex); vac[int("10" * (nS // 2), 2)] = 1
tS = np.arange(0, 3.01, 0.25)
for Jg, mg in ((0.0, 0.0), (0.5, 0.0), (2.0, 1.0), (4.0, 1.0)):       # J = 0 은 자유 페르미온(견줄 대상 아님), J = 0.5 는 상호작용이 있는 빠른 경우
    HS = (hopS + sp.diags(mg * (zS * sgn).sum(1) + Jg * (Lf ** 2).sum(1))).tocsc()
    Smax = 0.0; cmax = 0
    for ph in expm_multiply(-1j * HS, vac, start=0, stop=tS[-1], num=len(tS), endpoint=True):
        p = schmidt(ph, nS, nS // 2); Smax = max(Smax, entropy_bits(p)); cmax = max(cmax, chi_needed(p))
    print(f"R6) Schwinger {nS}사이트, J = {Jg}, m = {mg}: t ≤ 3 에서 반쪽 엔트로피 최대 {Smax:.2f} 비트, χ 최대 {cmax}")

# =====================================================================
# (E) 연습문제 검산 (맨 앞의 R1b 는 5절 (1)의 보강 확인)
# =====================================================================
print("\n" + "=" * 78); print("(E) 연습문제 검산 (+ R1b: 5절 (1) 보강)"); print("=" * 78)
rngE = np.random.default_rng(1919)
# R1b: 프레임을 4배로 (N = 4096, n = 12) — 5절 (1) 끝 문장의 근거 (난수 흐름을 바꾸지 않으려고 이 자리에 둔다)
N4 = 4096; t4 = np.arange(N4) / fs
c4 = sum(np.sin(2 * np.pi * k * f0 * t4 + 0.3 * k) / k for k in range(1, 11)); z4 = rngE.normal(size=N4)
z4 *= np.linalg.norm(c4) / np.linalg.norm(z4) / 10 ** (10 / 20)
P4 = np.abs(np.fft.fft((c4 + z4) * np.hanning(N4))) ** 2; P4 /= P4.sum(); k4 = int(round(f0 * N4 / fs))
succ4 = {}
for s in (16, 32, 64):
    ok = 0
    for _ in range(2000):
        ks = rngE.choice(N4, size=s, p=P4); ks = np.minimum(ks, N4 - ks); v, c = np.unique(ks, return_counts=True)
        ok += abs(v[np.argmax(c)] - k4) <= 2                                   # 구간 폭이 1/4 이 됐으므로 ±2 구간(±7.8 Hz)으로 잡았다
    succ4[s] = ok / 2000
s99_4 = next((s for s in (16, 32, 64) if succ4[s] >= 0.99), None)
print(f"R1b) N = 4096: 성공률 " + ", ".join(f"{s}샷 {v:.3f}" for s, v in succ4.items()) + f" → 99%에 {s99_4}샷, CNOT ≈ {s99_4} × {2**12-2} = {s99_4*(2**12-2):,};  FFT 연산 ≈ N·log₂N = {N4*12:,}")
# E1: 고전 기술의 크기 — 2^{2M} 대 FCI 차원
for name, M, ne, nq, dim in fci_rows[2:4]:
    print(f"E1) {name}: 2^{nq} = {2**nq:.3e}, FCI 차원 {dim:.3e}, 비 {2**nq/dim:,.0f}배")
# E2: 커널 값 두 개의 차이를 2σ 로 구분하는 샷 수 s ≥ 8K(1−K)/Δ²
r10 = kern_rows[-1]; s_need = 8 * r10[2] * (1 - r10[2]) / r10[3] ** 2
pairs = 100 * 99 // 2
print(f"E2) n = 10: K̄ = {r10[2]:.4f}, Δ = σ_K = {r10[3]:.4f} → s ≥ 8K(1−K)/Δ² = {s_need:,.0f} 샷/쌍, 발화 100개의 쌍 {pairs}개 → 회로 실행 {s_need*pairs:.2e}회")
K20 = r10[2] * 2.0 ** -10
print(f"    같은 추세(K̄, σ_K ∝ 2⁻ⁿ)면 n = 20 에서 K̄ ≈ {K20:.1e}, s ≈ 8/K̄ = {8/K20:.1e} 샷/쌍")
# E3: 같은 예산 단위로 맞춘 TTS
def tts_runs(P): return math.log(0.01) / math.log(1 - P)
q_train, q_shots_per_eval, Pq = 2000, 1000, 0.02
sa_evals, Ps = 5000, 0.30
q_naive = q_train + tts_runs(Pq); q_real = q_train * q_shots_per_eval + tts_runs(Pq); sa = math.ceil(tts_runs(Ps)) * sa_evals
print(f"E3) QAOA: 표본 추출 {tts_runs(Pq):.0f}샷; '평가 수'로만 세면 {q_naive:,.0f}, 회로 실행 수로 세면 {q_real:,.0f};  담금질: {tts_runs(Ps):.1f} → {math.ceil(tts_runs(Ps))}회 × {sa_evals} = {sa:,} 평가;  비율 {q_real/sa:.0f}배")

# ---------- 그림 ----------
fig, axs = plt.subplots(3, 2, figsize=(12, 12.6))
ax = axs[0, 0]
ax.plot([r[0] for r in pitch_rows], [r[1] for r in pitch_rows], "o-", color="#4C6EF5", lw=1.8)
ax.axhline(0.99, color="#868E96", ls=":", lw=1); ax.set_xscale("log", base=2); ax.set_xticks(shots_list); ax.set_xticklabels([str(s) for s in shots_list])
ax.set_xlabel("샷 수 (샷마다 1,024개 표본을 다시 적재)"); ax.set_ylabel("기본 주파수 구간을 맞힌 비율")
ax.set_title("(1) 피치 추정: QFT는 샷마다 적재 O(N), FFT는 한 번에 끝")

ax = axs[0, 1]
ns_ = [r[0] for r in kern_rows]
ax.plot(ns_, [r[4] for r in kern_rows], "o-", color="#F03E3E", lw=1.8, label="ZZ 양자 커널: 1-NN 정확도")
ax.plot(ns_, [r[6] for r in kern_rows], "s-", color="#4C6EF5", lw=1.8, label="RBF 커널: 1-NN 정확도")
ax.plot(ns_, [r[5] for r in kern_rows], "o--", color="#F03E3E", lw=1.2, label="ZZ 양자 커널: 정렬")
ax.plot(ns_, [r[7] for r in kern_rows], "s--", color="#4C6EF5", lw=1.2, label="RBF 커널: 정렬")
ax.set_xlabel("특성(큐비트) 수 n"); ax.set_ylabel("분리 지표"); ax.set_ylim(-0.05, 1.05); ax.legend(fontsize=8)
ax.set_title("(2) 화자 임베딩: 양자 커널은 RBF를 넘지 못한다")

ax = axs[1, 0]
ax.semilogy(ns_, [r[2] for r in kern_rows], "o-", color="#343A40", lw=1.8, label="비대각 커널 값의 평균")
ax.semilogy(ns_, [r[3] for r in kern_rows], "s-", color="#F59F00", lw=1.8, label="비대각 커널 값의 표준편차")
ax.semilogy(ns_, [2.0 ** -k for k in ns_], "k:", lw=1, label="1/2ⁿ (무작위 상태 쌍)")
ax.set_xlabel("특성(큐비트) 수 n"); ax.set_ylabel("커널 값 (로그)"); ax.legend(fontsize=8)
ax.set_title("(3) 커널 집중: 큐비트가 늘면 값이 0 근처로 몰린다")

ax = axs[1, 1]
ax.plot([r[0] for r in heis_rows], [r[1] for r in heis_rows], "o-", color="#4C6EF5", lw=1.8, label="얽힘 엔트로피 S (비트)")
ax.set_xlabel("시간 t"); ax.set_ylabel("반쪽 사슬 S (비트)", color="#4C6EF5")
ax2 = ax.twinx(); ax2.semilogy([r[0] for r in heis_rows], [r[2] for r in heis_rows], "s--", color="#F03E3E", lw=1.5, label="MPS 결합 차원 χ")
ax2.set_ylabel("필요 결합 차원 χ (로그)", color="#F03E3E")
ax.set_title("(4) Heisenberg 급랭 n = 16: 얽힘이 빠르게 자라 포화한다")

ax = axs[2, 0]
ax.semilogy([r[0] for r in hp_rows], [r[5] for r in hp_rows], "s--", color="#868E96", lw=1.5, label="턴 비트열 4^(N−1)")
ax.semilogy([r[0] for r in hp_rows], [r[1] for r in hp_rows], "o-", color="#2F9E44", lw=1.8, label="자기 회피 보행 수(유효한 접힘)")
ax.set_xlabel("잔기 수 N"); ax.set_ylabel("개수 (로그)"); ax.legend(fontsize=8)
ax.set_title("(5) HP 격자: 탐색 공간과 유효한 접힘의 수")

ax = axs[2, 1]
for name, c_, mk in (("1D Poisson", "#F03E3E", "o"), ("2D Poisson", "#F59F00", "s"), ("잘 조건화된 희소", "#4C6EF5", "D")):
    rows = [r for r in lin_rows if r[0] == name]
    ax.loglog([r[1] for r in rows], [r[4] for r in rows], mk + "-", color=c_, lw=1.8, label=f"{name}: CG 반복")
    ax.loglog([r[1] for r in rows], [r[3] for r in rows], mk + ":", color=c_, lw=1.2, label=f"{name}: 조건수 κ")
ax.set_xticks([256, 1024, 4096]); ax.set_xticklabels(["256", "1024", "4096"]); ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_xlabel("행렬 크기 N"); ax.set_ylabel("CG 반복 수 / 조건수 (로그)"); ax.legend(fontsize=7)
ax.set_title("(6) 선형계: 양자 비용은 κ에, 고전 CG는 √κ에 비례")
fig.tight_layout(); fig.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)), "ch19_results.png"), dpi=120)
print(f"\n   (실행 시간: 총 {time.time()-T0:.0f} s)")
