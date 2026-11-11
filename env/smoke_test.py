"""smoke_test.py — 설치된 SDK를 하나씩 불러 Bell 상태를 만들어 보고 표로 보고한다.
실행: python env/smoke_test.py            (설치 안 된 SDK는 "미설치"로 표시하고 넘어간다)
      python env/smoke_test.py --strict   (sdk.txt 의 SDK 중 하나라도 실패하면 종료 코드 1)

검사 내용
  1) Bell 상태 (|00> + |11>)/√2 의 확률이 {00: 0.5, 11: 0.5} 인가.
  2) 큐비트 순서: 큐비트 0에만 X를 건 2큐비트 상태를 각 SDK가 어떤 비트열로 보여 주는가(기대와 다르면 "값·순서 불일치").
     이 책은 빅 엔디안(q0가 왼쪽)이다. Qiskit은 리틀 엔디안이라 '01'로, Cirq·PennyLane·CUDA-Q·Braket은 '10'으로 보여야 한다(3장).
numpy 기준값은 SDK 없이 계산한다.
"""
import importlib.util
import sys

import numpy as np

H = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
X = np.array([[0, 1], [1, 0]])
I2 = np.eye(2)
CX = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])      # 빅 엔디안: 제어 q0, 표적 q1
REF_BELL = abs(CX @ np.kron(H, I2) @ np.eye(4)[0]) ** 2                        # [0.5, 0, 0, 0.5]
REF_X0 = abs(np.kron(X, I2) @ np.eye(4)[0]) ** 2                               # 빅 엔디안에서 '10'


def version(mod):
    try:
        from importlib.metadata import version as v
        return v(mod)
    except Exception:
        return "?"


def probs_from_counts(counts, shots):
    return {k.replace(" ", ""): v / shots for k, v in counts.items()}


def bell_ok(p, tol):
    """p: {'00': .., '11': ..} 꼴의 확률. 00 과 11 이 반반이고 나머지가 없는가."""
    return abs(p.get("00", 0) - 0.5) < tol and abs(p.get("11", 0) - 0.5) < tol and sum(p.values()) - p.get("00", 0) - p.get("11", 0) < tol


def t_numpy():
    assert np.allclose(REF_BELL, [0.5, 0, 0, 0.5]) and np.allclose(REF_X0, [0, 0, 1, 0])
    return np.__version__, True, "10 (기준: 빅 엔디안)"


def t_qiskit():
    from qiskit import QuantumCircuit
    from qiskit.quantum_info import Statevector
    qc = QuantumCircuit(2); qc.h(0); qc.cx(0, 1)
    p = {k: float(v) for k, v in Statevector(qc).probabilities_dict().items()}
    qx = QuantumCircuit(2); qx.x(0)
    order = max(Statevector(qx).probabilities_dict().items(), key=lambda kv: kv[1])[0]
    return version("qiskit"), bell_ok(p, 1e-9) and order == "01", f"{order} (리틀 엔디안: 01 이어야 한다)"


def t_aer():
    from qiskit import QuantumCircuit
    from qiskit_aer import AerSimulator
    qc = QuantumCircuit(2); qc.h(0); qc.cx(0, 1); qc.measure_all()
    counts = AerSimulator().run(qc, shots=4000, seed_simulator=7).result().get_counts()
    return version("qiskit-aer"), bell_ok(probs_from_counts(counts, 4000), 0.05), "샘플 4000샷"


def t_runtime():
    import qiskit_ibm_runtime  # noqa: F401
    from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2  # noqa: F401
    return version("qiskit-ibm-runtime"), True, "불러오기만 확인(계정 불필요)"


def t_pennylane():
    import pennylane as qp
    dev = qp.device("default.qubit", wires=2)

    @qp.qnode(dev)
    def bell():
        qp.Hadamard(wires=0); qp.CNOT(wires=[0, 1])
        return qp.probs(wires=[0, 1])

    @qp.qnode(dev)
    def x0():
        qp.PauliX(wires=0)
        return qp.probs(wires=[0, 1])

    p = np.asarray(bell()); order = format(int(np.argmax(np.asarray(x0()))), "02b")
    return version("pennylane"), bool(np.allclose(p, REF_BELL)) and order == "10", f"{order} (빅 엔디안: 10 이어야 한다)"


def t_cirq():
    import cirq
    q = cirq.LineQubit.range(2)
    sv = cirq.Simulator().simulate(cirq.Circuit(cirq.H(q[0]), cirq.CNOT(q[0], q[1]))).final_state_vector
    sx = cirq.Simulator().simulate(cirq.Circuit(cirq.X(q[0]), cirq.I(q[1]))).final_state_vector
    order = format(int(np.argmax(abs(sx) ** 2)), "02b")
    return version("cirq-core"), bool(np.allclose(abs(sv) ** 2, REF_BELL, atol=1e-6)) and order == "10", f"{order} (빅 엔디안: 10 이어야 한다)"


def t_qutip():
    import qutip
    from qutip import basis, tensor
    bell = (tensor(basis(2, 0), basis(2, 0)) + tensor(basis(2, 1), basis(2, 1))).unit()
    p = abs(bell.full().ravel()) ** 2
    x0 = tensor(basis(2, 1), basis(2, 0)).full().ravel()
    order = format(int(np.argmax(abs(x0) ** 2)), "02b")
    return qutip.__version__, bool(np.allclose(p, REF_BELL)) and order == "10", f"{order} (tensor 순서: 10 이어야 한다)"


def t_qutip_qip():
    import qutip_qip  # noqa: F401
    from qutip_qip.circuit import QubitCircuit
    qc = QubitCircuit(2); qc.add_gate("SNOT", targets=0); qc.add_gate("CNOT", controls=0, targets=1)
    return version("qutip-qip"), True, "회로 객체 생성만 확인"


def t_qualtran():
    import qualtran  # noqa: F401
    return version("qualtran"), True, "불러오기만 확인"


def t_cudaq():
    import cudaq
    k = cudaq.make_kernel(); q = k.qalloc(2); k.h(q[0]); k.cx(q[0], q[1]); k.mz(q)
    r = cudaq.sample(k, shots_count=4000)
    counts = {bits: n for bits, n in r.items()}
    kx = cudaq.make_kernel(); qx = kx.qalloc(2); kx.x(qx[0]); kx.mz(qx)
    order = max(cudaq.sample(kx, shots_count=100).items(), key=lambda kv: kv[1])[0]
    return version("cudaq"), bell_ok(probs_from_counts(counts, 4000), 0.05) and order == "10", f"{order} (q0 가 왼쪽: 10 이어야 한다); 대상 {cudaq.get_target().name}"


def t_braket():
    from braket.circuits import Circuit
    from braket.devices import LocalSimulator
    r = LocalSimulator().run(Circuit().h(0).cnot(0, 1), shots=4000).result()
    rx = LocalSimulator().run(Circuit().x(0).i(1), shots=100).result()
    order = max(rx.measurement_counts.items(), key=lambda kv: kv[1])[0]
    return version("amazon-braket-sdk"), bell_ok(probs_from_counts(dict(r.measurement_counts), 4000), 0.05) and order == "10", f"{order} (빅 엔디안: 10 이어야 한다)"


DIST = {"cirq": "cirq-core", "braket": "amazon-braket-sdk"}      # import 이름과 배포 이름이 다른 것
TESTS = [  # (표시 이름, import 이름, 함수, sdk.txt 에 들어 있는가)
    ("numpy (기준)", "numpy", t_numpy, True),
    ("Qiskit", "qiskit", t_qiskit, True),
    ("Qiskit Aer", "qiskit_aer", t_aer, True),
    ("qiskit-ibm-runtime", "qiskit_ibm_runtime", t_runtime, True),
    ("PennyLane", "pennylane", t_pennylane, True),
    ("Cirq", "cirq", t_cirq, True),
    ("QuTiP", "qutip", t_qutip, True),
    ("qutip-qip", "qutip_qip", t_qutip_qip, True),
    ("Qualtran", "qualtran", t_qualtran, True),
    ("CUDA-Q", "cudaq", t_cudaq, False),
    ("Amazon Braket SDK", "braket", t_braket, False),
]


def main():
    strict = "--strict" in sys.argv
    print(f"Python {sys.version.split()[0]}  ({sys.executable})")
    rows, bad = [], 0
    for name, mod, fn, required in TESTS:
        if importlib.util.find_spec(mod) is None:
            rows.append((name, "미설치", "-", "-")); bad += required
            continue
        try:
            ver, ok, note = fn()
            rows.append((name, "통과" if ok else "값·순서 불일치", ver, note)); bad += (not ok)
        except Exception as e:  # SDK 내부 오류도 표에 남긴다
            rows.append((name, "오류", version(DIST.get(mod, mod.replace("_", "-"))), f"{type(e).__name__}: {str(e)[:90]}")); bad += 1
    w = max(len(r[0]) for r in rows)
    print(f"\n{'SDK'.ljust(w)}  상태        버전        비고")
    for name, status, ver, note in rows:
        print(f"{name.ljust(w)}  {status.ljust(10)}  {str(ver).ljust(10)}  {note}")
    missing = [r[0] for r in rows if r[1] == "미설치"]
    if missing:
        print("\n미설치:", ", ".join(missing), "→ bash env/setup.sh (CUDA-Q는 --gpu, Braket은 --hardware)")
    sys.exit(1 if (strict and bad) else 0)


if __name__ == "__main__":
    main()
