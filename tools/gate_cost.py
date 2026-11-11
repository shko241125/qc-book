"""gate_cost.py — 게이트 리스트로 '세 숫자'(2큐비트 게이트 수 / 깊이 / T-count)를 세는 작은 도구.
실행(자체 검사): python3 tools/gate_cost.py
환경: Python 3.13, numpy 2.5 (SDK 없이 동작). 큐비트 순서: 빅 엔디안(3장 규약, 큐비트 q = 텐서 축 q).

규약(2부 '게이트 비용 기록'의 규칙):
  - 기저: CX(CNOT) + 임의 1큐비트 게이트, 전결합(all-to-all). CZ 는 2큐비트 게이트 1개로 센다.
  - 전개: SWAP = 3 CX; CP(θ) = 2 CX + 위상 회전 3개(각 큐비트에 θ/2, 두 CX 사이에 −θ/2), 단 CP(π) = CZ, CP(0) = 생략;
          제어 Rz/Ry/Rx = 2 CX + 회전 2개; exp(−iθ Z⊗Z) = 2 CX + Rz(2θ); Toffoli/CCZ = 6 CX + 7 T(4장·9장); Fredkin = Toffoli + 2 CX.
  - 2큐비트 게이트 수: 전개 후 CX·CZ 의 개수.
  - 깊이: 전개한 회로의 게이트를 리스트 순서대로 가능한 한 앞에(ASAP) 놓았을 때의 층 수. 1큐비트·2큐비트 게이트 모두 한 층, 측정 제외.
  - T-count: 전개 후 모든 게이트가 Clifford 이거나 π/4 의 배수 회전이면 정확(T = 1, Toffoli = 7).
          π/4 의 배수가 아닌 회전은 개수 R 로 세고 T ≈ 3·log2(1/ε)·R (Ross & Selinger 2016, arXiv:1403.2975; ε = 1e-6 이면 회전당 약 60).

게이트 = (이름, 큐비트 튜플, 각도 또는 None).
  1큐비트: h x y z s sdg t tdg | rx ry rz p(각도)        2큐비트 기본: cx cz
  합성(전개 대상): swap cp crz cry crx rzz ccx ccz cswap
"""
import numpy as np

_H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
_FIXED = {
    "h": _H,
    "x": np.array([[0, 1], [1, 0]], dtype=complex),
    "y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "z": np.diag([1, -1]).astype(complex),
    "s": np.diag([1, 1j]).astype(complex),
    "sdg": np.diag([1, -1j]).astype(complex),
    "t": np.diag([1, np.exp(1j * np.pi / 4)]),
    "tdg": np.diag([1, np.exp(-1j * np.pi / 4)]),
}
ROT_1Q = ("rx", "ry", "rz", "p")
PRIMITIVE_2Q = ("cx", "cz")
MERGEABLE = ("rx", "ry", "rz", "p", "rzz", "cp", "crz", "cry", "crx")


def G(name, *qubits, angle=None):
    """게이트 튜플을 만든다: G('cx', 0, 1), G('rz', 2, angle=0.3)."""
    return (name, tuple(int(q) for q in qubits), None if angle is None else float(angle))


def _mat1(name, angle):
    if name in _FIXED:
        return _FIXED[name]
    c, s = np.cos(angle / 2), np.sin(angle / 2)
    if name == "rx": return np.array([[c, -1j * s], [-1j * s, c]])
    if name == "ry": return np.array([[c, -s], [s, c]], dtype=complex)
    if name == "rz": return np.diag([np.exp(-1j * angle / 2), np.exp(1j * angle / 2)])
    if name == "p": return np.diag([1, np.exp(1j * angle)])
    raise ValueError(f"unknown 1-qubit gate {name}")


def _is_multiple(angle, unit, tol=1e-9):
    m = angle / unit
    return abs(m - round(m)) < tol


def _toffoli(a, b, t, with_h=True):
    """N&C 그림 4.9 의 Toffoli: CX 6, T/T† 7, H 2. with_h=False 이면 CCZ."""
    core = [G("cx", b, t), G("tdg", t), G("cx", a, t), G("t", t), G("cx", b, t), G("tdg", t), G("cx", a, t),
            G("t", b), G("t", t)]
    tail = [G("cx", a, b), G("t", a), G("tdg", b), G("cx", a, b)]
    return ([G("h", t)] + core + [G("h", t)] + tail) if with_h else (core + tail)


def expand(gates):
    """합성 게이트를 CX/CZ + 1큐비트 게이트로 푼다. 2π 의 배수인 회전(전역 위상뿐)은 버린다."""
    out = []
    for name, q, ang in gates:
        if name in _FIXED or name in PRIMITIVE_2Q:
            out.append((name, q, None))
        elif name in ROT_1Q:
            if not _is_multiple(ang, 2 * np.pi): out.append((name, q, ang))
        elif name == "swap":
            a, b = q; out += [G("cx", a, b), G("cx", b, a), G("cx", a, b)]
        elif name == "cp":
            c, t = q
            if _is_multiple(ang, 2 * np.pi): continue                      # CP(0) = I
            if _is_multiple(ang - np.pi, 2 * np.pi): out.append(G("cz", c, t)); continue   # CP(π) = CZ
            out += expand([G("p", c, angle=ang / 2), G("p", t, angle=ang / 2), G("cx", c, t), G("p", t, angle=-ang / 2), G("cx", c, t)])
        elif name in ("crz", "cry"):
            c, t = q; r = name[1:]
            if _is_multiple(ang, 4 * np.pi): continue
            out += expand([G(r, t, angle=ang / 2), G("cx", c, t), G(r, t, angle=-ang / 2), G("cx", c, t)])
        elif name == "crx":
            c, t = q
            if _is_multiple(ang, 4 * np.pi): continue
            out += [G("h", t)] + expand([G("crz", c, t, angle=ang)]) + [G("h", t)]
        elif name == "rzz":                                                   # exp(-i θ Z⊗Z)
            a, b = q
            if _is_multiple(ang, np.pi): continue                           # ±I
            out += expand([G("cx", a, b), G("rz", b, angle=2 * ang), G("cx", a, b)])
        elif name == "ccx":
            out += _toffoli(*q, with_h=True)
        elif name == "ccz":
            out += _toffoli(*q, with_h=False)
        elif name == "cswap":                                                 # Fredkin(c; a, b) = CX(b,a)·Toffoli(c,a→b)·CX(b,a)
            c, a, b = q; out += [G("cx", b, a)] + _toffoli(c, a, b) + [G("cx", b, a)]
        else:
            raise ValueError(f"unknown gate {name}")
    return out


_INV_NAME = {"t": "tdg", "tdg": "t", "s": "sdg", "sdg": "s"}


def inverse(gates):
    """게이트 리스트의 역(순서를 뒤집고 각도 부호를 바꾼다; h·x·y·z·cx·cz·swap·ccx·ccz·cswap 은 자기 역)."""
    return [(_INV_NAME.get(name, name), q, None if ang is None else -ang) for name, q, ang in reversed(gates)]


def two_qubit_count(gates):
    return sum(1 for name, _, _ in expand(gates) if name in PRIMITIVE_2Q)


def depth(gates):
    """ASAP 깊이: 리스트 순서대로, 각 게이트를 그 큐비트들이 비는 가장 이른 층에 놓는다."""
    level = {}
    for _, q, _ in expand(gates):
        l = max(level.get(x, 0) for x in q) + 1
        for x in q: level[x] = l
    return max(level.values()) if level else 0


def t_count(gates):
    """(정확한 T 수, π/4 의 배수가 아닌 회전 수 R). R = 0 이면 T 수가 정확한 T-count."""
    T = R = 0
    for name, _, ang in expand(gates):
        if name in ("t", "tdg"): T += 1
        elif name in ROT_1Q:
            if _is_multiple(ang, np.pi / 4): T += int(round(ang / (np.pi / 4))) % 2    # π/4 의 홀수 배 = Clifford·T
            else: R += 1
    return T, R


def t_estimate(R, eps=1e-6):
    """회전 R 개를 정밀도 ε 로 합성할 때의 T 수 추정 3·log2(1/ε)·R (Ross & Selinger 2016)."""
    return int(round(3 * np.log2(1 / eps) * R))


def count_ops(gates):
    ops = {}
    for name, _, _ in expand(gates): ops[name] = ops.get(name, 0) + 1
    return ops


def cost(gates):
    T, R = t_count(gates)
    return {"two_q": two_qubit_count(gates), "depth": depth(gates), "T": T, "R": R}


def fmt(gates, eps=1e-6):
    """한 줄 요약: '2큐비트 23, 깊이 31, T 9, 회전 R 12 (T ≈ 9 + 718)'."""
    c = cost(gates)
    s = f"2큐비트 {c['two_q']}, 깊이 {c['depth']}, T {c['T']}"
    if c["R"]:
        s += f", 회전 R {c['R']} (T ≈ {c['T']} + {t_estimate(c['R'], eps)})" if c["T"] else f", 회전 R {c['R']} (T ≈ {t_estimate(c['R'], eps)})"
    else:
        s += " (정확)"
    return s


def merge_adjacent(gates):
    """같은 큐비트들에 바로 이어 붙은 같은 종류의 회전(rx·ry·rz·p·rzz·cp·제어 회전)을 각도를 더해 하나로 합친다.
    '바로 이어 붙은' = 그 큐비트들의 직전 게이트가 모두 같은 게이트. 교환하는 게이트를 건너뛰는 재배열은 하지 않는다."""
    out, last = [], {}
    for name, q, ang in gates:
        prev = {last.get(x) for x in q}
        if name in MERGEABLE and len(prev) == 1 and None not in prev:
            i = prev.pop()
            if out[i][0] == name and out[i][1] == q:
                out[i] = (name, q, out[i][2] + ang); continue
        out.append((name, q, ang))
        for x in q: last[x] = len(out) - 1
    return out


def apply(gates, psi, n):
    """상태벡터(또는 (2^n, k) 배열의 각 열)에 게이트 리스트를 순서대로 적용. 합성 게이트는 전개해서 적용한다."""
    tail = psi.shape[1:]
    psi = np.asarray(psi, dtype=complex).reshape([2] * n + list(tail))
    for name, q, ang in expand(gates):
        if name == "cx":
            c, t = q; idx = [slice(None)] * psi.ndim; idx[c] = 1
            sub = psi[tuple(idx)]; ax = t if t < c else t - 1
            psi[tuple(idx)] = np.flip(sub, axis=ax)
        elif name == "cz":
            a, b = q; idx = [slice(None)] * psi.ndim; idx[a] = 1; idx[b] = 1
            psi[tuple(idx)] = -psi[tuple(idx)]
        else:
            psi = np.moveaxis(np.tensordot(_mat1(name, ang), psi, axes=([1], [q[0]])), 0, q[0])
    return psi.reshape([2 ** n] + list(tail))


def unitary(gates, n):
    """게이트 리스트의 2^n × 2^n 유니터리(검증용)."""
    return apply(gates, np.eye(2 ** n, dtype=complex), n)


def same_up_to_phase(A, B, tol=1e-9):
    """행렬(또는 벡터) A, B 가 전역 위상만 다르고 같은가."""
    A = np.asarray(A, dtype=complex); B = np.asarray(B, dtype=complex)
    if A.shape != B.shape: return False
    k = np.argmax(abs(B))
    if abs(B.flat[k]) < tol: return bool(np.allclose(A, B, rtol=0, atol=tol))
    ph = A.flat[k] / B.flat[k]
    return bool(abs(abs(ph) - 1) < tol and np.allclose(A, ph * B, rtol=0, atol=tol))


# ── 자체 검사 ───────────────────────────────────────────────────────
if __name__ == "__main__":
    def perm_matrix(f, n):
        U = np.zeros((2 ** n, 2 ** n), dtype=complex)
        for y in range(2 ** n): U[f(y), y] = 1
        return U
    def bits(y, n): return [(y >> (n - 1 - k)) & 1 for k in range(n)]
    def val(b): return int("".join(map(str, b)), 2)

    ok = True
    def check(label, cond, extra=""):
        global ok
        ok = ok and bool(cond)
        print(f"   {label}: {'OK' if cond else 'FAIL'}{extra}")

    print("gate_cost.py 자체 검사 (빅 엔디안, 전개 = CX/CZ + 1큐비트)")
    # 1) CP(θ): diag(1,1,1,e^{iθ}) — 2 CX + 위상 3
    th = 0.7; cp = [G("cp", 0, 1, angle=th)]
    check("CP(0.7) = diag(1,1,1,e^{iθ})", same_up_to_phase(unitary(cp, 2), np.diag([1, 1, 1, np.exp(1j * th)])), f"  → {fmt(cp)}")
    cs = [G("cp", 0, 1, angle=np.pi / 2)]
    check("CP(π/2) = 제어 S: T 3 (정확)", same_up_to_phase(unitary(cs, 2), np.diag([1, 1, 1, 1j])) and t_count(cs) == (3, 0), f"  → {fmt(cs)}")
    czg = [G("cp", 0, 1, angle=np.pi)]
    check("CP(π) = CZ: 2큐비트 1", same_up_to_phase(unitary(czg, 2), np.diag([1, 1, 1, -1])) and two_qubit_count(czg) == 1, f"  → {fmt(czg)}")
    # 2) Toffoli: 6 CX, 7 T, 유니터리 = 치환
    tof = [G("ccx", 0, 1, 2)]
    U_tof = perm_matrix(lambda y: y ^ 1 if (y >> 1) == 3 else y, 3)
    check("Toffoli = 치환 행렬, CX 6, T 7", same_up_to_phase(unitary(tof, 3), U_tof) and two_qubit_count(tof) == 6 and t_count(tof) == (7, 0), f"  → {fmt(tof)}")
    ccz = [G("ccz", 0, 1, 2)]
    check("CCZ = diag(1,…,1,−1), CX 6, T 7", same_up_to_phase(unitary(ccz, 3), np.diag([1] * 7 + [-1])) and cost(ccz)["two_q"] == 6 and t_count(ccz) == (7, 0), f"  → {fmt(ccz)}")
    # 3) SWAP = 3 CX, Fredkin = Toffoli + 2 CX
    sw = [G("swap", 0, 2)]
    check("SWAP(0,2) = 치환 행렬, CX 3", same_up_to_phase(unitary(sw, 3), perm_matrix(lambda y: val([bits(y, 3)[2], bits(y, 3)[1], bits(y, 3)[0]]), 3)) and two_qubit_count(sw) == 3, f"  → {fmt(sw)}")
    fr = [G("cswap", 0, 1, 2)]
    check("Fredkin = 제어 SWAP, CX 8, T 7", same_up_to_phase(unitary(fr, 3), perm_matrix(lambda y: val([1, bits(y, 3)[2], bits(y, 3)[1]]) if bits(y, 3)[0] else y, 3)) and cost(fr)["two_q"] == 8 and t_count(fr) == (7, 0), f"  → {fmt(fr)}")
    # 4) exp(-iθ ZZ), 제어 회전
    zz = [G("rzz", 0, 1, angle=th)]
    check("exp(−iθ Z⊗Z) = CX·Rz(2θ)·CX", same_up_to_phase(unitary(zz, 2), np.diag(np.exp(-1j * th * np.array([1, -1, -1, 1])))), f"  → {fmt(zz)}")
    for r in ("rz", "ry", "rx"):
        cr = [G("c" + r, 0, 1, angle=th)]; blk = np.eye(4, dtype=complex); blk[2:, 2:] = _mat1(r, th)
        check(f"제어 {r.capitalize()}(θ) = 2 CX + 회전 2", same_up_to_phase(unitary(cr, 2), blk) and two_qubit_count(cr) == 2, f"  → {fmt(cr)}")
    seqg = [G("h", 0), G("t", 0), G("cp", 0, 1, angle=0.3), G("s", 1), G("ccx", 0, 1, 2), G("ry", 2, angle=0.4), G("cswap", 2, 0, 1)]
    check("inverse: U·U† = I", same_up_to_phase(unitary(seqg + inverse(seqg), 3), np.eye(8)))
    # 5) T-count 분류: π/4 의 홀수 배 = 1, π/2 의 배수 = 0, 그 밖 = 회전 R
    check("t_count 분류", t_count([G("rz", 0, angle=np.pi / 4), G("p", 0, angle=3 * np.pi / 4), G("ry", 0, angle=np.pi / 2), G("rx", 0, angle=0.3), G("t", 0)]) == (3, 1))
    # 6) merge_adjacent: 이어 붙은 같은 회전만 합친다
    mg = merge_adjacent([G("rx", 0, angle=0.1), G("rx", 0, angle=0.2), G("rzz", 0, 1, angle=0.1), G("rzz", 0, 1, angle=0.1), G("rzz", 1, 2, angle=0.1), G("rzz", 0, 1, angle=0.1)])
    check("merge_adjacent", [g[0] for g in mg] == ["rx", "rzz", "rzz", "rzz"] and abs(mg[0][2] - 0.3) < 1e-12 and abs(mg[1][2] - 0.2) < 1e-12)
    # 7) BV (8장: n = 6, s = 101101): 2큐비트 4, 깊이 7, T 0
    n, s = 6, "101101"
    bv = [G("h", q) for q in range(n)] + [G("x", n), G("h", n)] + [G("cx", q, n) for q in range(n) if s[q] == "1"] + [G("h", q) for q in range(n)]
    psi = apply(bv, np.eye(2 ** (n + 1), dtype=complex)[0], n + 1).reshape(2 ** n, 2)
    p_in = (abs(psi) ** 2).sum(axis=1)
    check("BV(n = 6, s = 101101): 측정 = s (확률 1), 2큐비트 4, 깊이 7, T 0",
          abs(p_in[int(s, 2)] - 1) < 1e-12 and cost(bv) == {"two_q": 4, "depth": 7, "T": 0, "R": 0}, f"  → {fmt(bv)}")
    # 8) Grover n = 3 반복 1회 (9장: 표시 |101>): 오라클 X·CCZ·X + 확산 H X CCZ X H
    w = "101"
    orc = [G("x", q) for q in range(3) if w[q] == "0"] + [G("ccz", 0, 1, 2)] + [G("x", q) for q in range(3) if w[q] == "0"]
    dif = [G("h", q) for q in range(3)] + [G("x", q) for q in range(3)] + [G("ccz", 0, 1, 2)] + [G("x", q) for q in range(3)] + [G("h", q) for q in range(3)]
    it = orc + dif
    O = np.eye(8, dtype=complex); O[int(w, 2), int(w, 2)] = -1
    D = 2 * np.full((8, 8), 1 / 8) - np.eye(8)
    check("Grover n = 3 반복 1회 = D·O (전역 위상 제외), 2큐비트 12, T 14", same_up_to_phase(unitary(it, 3), D @ O) and cost(it)["two_q"] == 12 and t_count(it) == (14, 0), f"  → {fmt(it)}")
    print("   k = 2 전체(초기 H 3 포함):", fmt([G("h", q) for q in range(3)] + it + it))
    print("자체 검사:", "모두 통과" if ok else "실패 있음")
    raise SystemExit(0 if ok else 1)
