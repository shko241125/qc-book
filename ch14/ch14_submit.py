"""14장 — IBM Open Plan 제출 스크립트. qiskit==2.6.2, qiskit-ibm-runtime==0.52.0, qiskit-aer==0.18.2 — 실행 완료(2026-11-11, fake·real)
  python ch14_submit.py fake           가짜 Heron(FakeMarrakesh) 보정값으로 로컬 실행 — QPU 를 쓰지 않는다
  python ch14_submit.py real           실제 장치를 고르고 트랜스파일해 제출 계획만 찍는다(제출하지 않는다)
  python ch14_submit.py real --yes     실제로 제출한다(Open 플랜 분량 약 10초)
Runtime 인터페이스: qiskit-ibm-runtime 0.52 의 최상위 Sampler·Estimator 는 클라이언트 쪽 프리미티브(executor_sampler·
executor_estimator)다. 0.50 부터 폐기 예정인 것은 서버 쪽 SamplerV2·EstimatorV2(qiskit_ibm_runtime.server_side)다.

순서: (0) 가짜 백엔드로 로컬 노이즈 에뮬레이션 → (1) 실제 백엔드 선택·트랜스파일 → (2) 제출 → (3) 결과·사용 시간 기록.
Open Plan 한도: 28일 롤링 윈도우 10분 QPU 시간. 사용량은 job.usage() (초)로 확인한다
(metrics()['usage'] 의 키는 runtime 0.51 에서 quantum_seconds → qpu_charge_time_seconds 로 바뀌었다).
"""
import json, pathlib, sys, time
import numpy as np
from qiskit import QuantumCircuit
from qiskit.transpiler import generate_preset_pass_manager
from qiskit.quantum_info import SparsePauliOp

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "env"))
from run_ibm import two_qubit_gates, gate_depth   # noqa: E402  깊이·2큐비트 수를 run_ibm·gate_cost 와 같은 규약으로(barrier·측정 제외)


# ---------- 회로 5개 (ch14_hardware.py 의 numpy 버전과 같은 논리 회로) ----------
def ghz(n=5):
    qc = QuantumCircuit(n, n); qc.h(0)
    for i in range(n - 1): qc.cx(i, i + 1)
    qc.measure(range(n), range(n)[::-1]); return qc     # 다른 회로처럼 책 순서(왼쪽이 q0)로 읽는다. 2026-11-11 실행분은 measure_all(Qiskit 순서)


def bv(s="1011"):
    n = len(s); qc = QuantumCircuit(n + 1, n)
    qc.h(range(n)); qc.x(n); qc.h(n)
    for q, bit in enumerate(s):
        if bit == "1": qc.cx(q, n)          # 빅 엔디안 s: q0 ↔ s[0]
    qc.h(range(n)); qc.measure(range(n), range(n)[::-1])   # 비트열을 s 그대로 읽기 위해 역순 매핑
    return qc


def grover(target="101", iters=2):
    n = 3; qc = QuantumCircuit(n, n); qc.h(range(n))
    flips = [q for q, b in enumerate(target) if b == "0"]
    for _ in range(iters):
        qc.x(flips); qc.ccz(0, 1, 2); qc.x(flips)          # 오라클
        qc.h(range(n)); qc.x(range(n)); qc.ccz(0, 1, 2); qc.x(range(n)); qc.h(range(n))   # 확산
    qc.measure(range(n), range(n)[::-1]); return qc


def qpe_t(t=3):
    qc = QuantumCircuit(t + 1, t); qc.x(t); qc.h(range(t))
    for j in range(t):                                   # q[j] 가 U^{2^j} (리틀 엔디안 계수)
        qc.cp(2 * np.pi * (2 ** j) / 8, j, t)
    from qiskit.circuit.library import QFTGate
    qc.append(QFTGate(t).inverse(), range(t))            # QFTGate 는 리틀 엔디안 정수 x = Σ q_j 2^j 를 내놓는다
    # → 고전 비트 c_j = q_j 로 그대로 매핑하면 Qiskit 문자열(c2 c1 c0)이 곧 x 의 이진 표기. (BV·Grover 와 달리 역순 매핑 없음)
    qc.measure(range(t), range(t)); return qc            # 결과 '001' = φ·8 = 1


def hva_tfim(gamma=0.2064, beta=0.3927, n=4):
    qc = QuantumCircuit(n); qc.h(range(n))
    for q in range(n - 1): qc.rzz(-2 * gamma, q, q + 1)  # e^{-iγ(-ZZ)} = RZZ(-2γ)
    qc.rx(-2 * beta, range(n))                           # e^{-iβ(-X)} = RX(-2β)
    return qc                                            # 측정은 Estimator 가 한다


H_TFIM = SparsePauliOp.from_sparse_list(
    [("ZZ", [i, i + 1], -1.0) for i in range(3)] + [("X", [i], -1.0) for i in range(4)], num_qubits=4)

CIRCUITS = {"ghz5": ghz(), "bv1011": bv(), "grover101": grover(), "qpe_t": qpe_t()}


def run(mode="fake", shots=4000, yes=False):
    if mode == "fake":
        # (0) 로컬 노이즈 에뮬레이션: Heron 가짜 백엔드의 보정 데이터로 Aer 노이즈 모델을 만든다
        from qiskit_ibm_runtime.fake_provider import FakeMarrakesh
        from qiskit_aer.primitives import SamplerV2, EstimatorV2
        backend = FakeMarrakesh()
        sampler, estimator = SamplerV2.from_backend(backend, seed=7), EstimatorV2.from_backend(backend)
    else:
        # (1) 실제 백엔드: Open Plan 인스턴스는 us-east 만; least_busy 로 대기열이 짧은 Heron 을 고른다
        from qiskit_ibm_runtime import Sampler, Estimator
        from run_ibm import get_service                  # 환경 변수 → .env → 저장 계정 순서(run_ibm.py 와 같다)
        service = get_service()
        backend = service.least_busy(operational=True, simulator=False, min_num_qubits=20,
                                     filters=lambda b: (getattr(b, "processor_type", None) or {}).get("family") == "Heron")   # Heron 으로 한정
        sampler, estimator = Sampler(mode=backend), Estimator(mode=backend)
        # 판독 오류 완화(TREX)는 Estimator 의 기본값(resilience_level=1). 이 장의 '노이즈 예측'과 같은 조건으로
        # 비교하려면 끄고, 완화된 값은 '판독 오류 제외 예측'(−4.555)과 비교한다.
        estimator.options.resilience_level = 0
        # 세션(Session) 은 Open Plan 에서 불가 → 작업 모드 또는 배치 모드로 제출
    print("backend:", backend.name)

    pm = generate_preset_pass_manager(optimization_level=3, backend=backend, seed_transpiler=7)
    isa = {k: pm.run(qc) for k, qc in CIRCUITS.items()}
    hva_isa = pm.run(hva_tfim()); H_isa = H_TFIM.apply_layout(hva_isa.layout)

    def stats(qc):                                       # 14장 표·부록 A 에 쓰는 트랜스파일 지표와 물리 큐비트
        used = sorted({qc.find_bit(q).index for i in qc.data if i.operation.name != "barrier" for q in i.qubits})
        return {"two_qubit_gates": two_qubit_gates(qc), "depth": gate_depth(qc), "ops": {k: int(v) for k, v in qc.count_ops().items()},
                "physical_qubits": used}
    st = {k: stats(qc) for k, qc in isa.items()}; st["hva"] = stats(hva_isa)
    for k, v in st.items():
        print(f"  {k:10s}: 2q={v['two_qubit_gates']}, depth={v['depth']}(게이트만), qubits={v['physical_qubits']}, ops={v['ops']}")

    if mode != "fake" and not yes:
        print(f"\n제출 계획: {backend.name} · Sampler PUB {len(isa)}개 × {shots}샷 + Estimator 1개(precision 0.02). "
              "2026-11-11 같은 구성의 실행은 QPU 6+3초를 썼다. 제출하려면 --yes 를 붙인다.")
        return None

    record = {"backend": backend.name, "shots": shots, "t": time.strftime("%Y-%m-%d %H:%M"), "bit_order": "book", "jobs": {}}

    def dump():
        with open(f"ch14_record_{mode}.json", "w") as fh:
            json.dump(record, fh, indent=1, ensure_ascii=False)

    job = sampler.run(list(isa.values()), shots=shots)   # 회로 4개를 한 작업(PUB 4개)으로 → 오버헤드 1회
    jid = job.job_id() if hasattr(job, "job_id") else "fake"
    print("sampler job id:", jid, flush=True)            # 끊겨도 service.job(id) 로 되찾는다
    # HVA 에너지: Estimator 가 Z·X 기저를 나눠 잰다 (관측량 그룹 2개). 두 작업을 먼저 다 제출하고 id 를 남긴 뒤 결과를 기다린다
    ejob = estimator.run([(hva_isa, H_isa)], precision=0.02)   # 목표 표준오차. 실제 표준오차는 결과의 stds(2026-11-11: 0.039)
    ejid = ejob.job_id() if hasattr(ejob, "job_id") else "fake"
    print("estimator job id:", ejid, flush=True)
    record["job_ids"] = {"sampler": jid, "estimator": ejid}
    dump()                                               # 결과를 받기 전에 끊겨도 두 id 가 파일에 남는다

    res = job.result()
    for k, pub in zip(isa.keys(), res):
        counts = pub.data.c.get_counts()
        record["jobs"][k] = {"job_id": jid, "counts": counts, **st[k]}
        print(k, sorted(counts.items(), key=lambda kv: -kv[1])[:3])
    d = ejob.result()[0].data
    ev, sd = float(d.evs), float(np.asarray(getattr(d, "stds", np.nan)))
    # fake 모드의 Aer Estimator 는 판독 오류 없이 기댓값을 계산한다 → 비교값은 −4.555, real(resilience 0) 은 −4.433
    print(f"HVA <H> = {ev:.4f} ± {sd:.4f} (이상 -4.6462 | 노이즈 예측: 판독 오류 포함 -4.433, 제외 -4.555)")
    record["jobs"]["hva"] = {"job_id": ejid, "evs": ev, "stds": sd, **st["hva"]}
    dump()                                               # 사용 시간 조회가 실패해도 측정 결과는 남는다
    if mode != "fake":                                   # (3) 사용 시간 기록 — 14장 표의 '사용 분' 열
        for j in (job, ejob):
            try:
                u = j.metrics().get("usage", {})
                secs = u.get("qpu_charge_time_seconds", u.get("quantum_seconds"))   # runtime 0.51 에서 키 이름이 바뀌었다
                if secs is None:
                    secs = j.usage()                     # RuntimeJobV2.usage() → 초(float). 집계 전에는 0 일 수 있다
            except Exception as e:                       # 사용량 때문에 기록을 잃지 않는다
                secs = f"조회 실패: {type(e).__name__}: {e}"
            record.setdefault("qpu_seconds", []).append(secs)
        print("qpu_seconds:", record["qpu_seconds"])
        dump()
        keep = pathlib.Path(__file__).resolve().parent.parent / "runs" / f"{time.strftime('%Y%m%d-%H%M%SZ', time.gmtime())}_{backend.name}_ch14_record.json"
        keep.write_text(json.dumps(record, indent=1, ensure_ascii=False), encoding="utf-8")   # 다음 실행이 덮어쓰지 않는 고정 사본
        print("고정 사본:", keep.name)
    return record


if __name__ == "__main__":
    args = sys.argv[1:]
    mode = next((a for a in args if a in ("fake", "real")), "fake")
    run(mode=mode, yes="--yes" in args)
