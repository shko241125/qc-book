"""run_ibm.py — IBM Quantum 장치(또는 로컬 시뮬레이터)에서 작은 회로를 돌리고 결과를 runs/ 에 남긴다.
작성됨, 미실행 — 이 스크립트를 쓴 환경에는 qiskit 이 없어 문법 검사만 했다. 첫 실행은 --local 로 한다.

  python env/run_ibm.py --local                      계정 없이 Aer 시뮬레이터로 같은 흐름을 확인한다
  python env/run_ibm.py --save                       API 키와 인스턴스 CRN 을 물어 이 컴퓨터에 저장한다(한 번만)
  python env/run_ibm.py --list                       내 인스턴스에서 쓸 수 있는 장치와 대기열
  python env/run_ibm.py --circuit bell               제출 계획만 보여 준다(--yes 가 없으면 제출하지 않는다)
  python env/run_ibm.py --circuit bell --yes         가장 한가한 장치에 1,024샷 제출
  python env/run_ibm.py --circuit ghz --n 5 --backend <이름> --shots 2048 --yes
  python env/run_ibm.py --fetch <job id>             전에 제출한 작업의 결과를 다시 받는다(제출 때 남긴 runs/pending_<id>.json 을 읽는다)

인증 순서: 환경 변수 QISKIT_IBM_TOKEN(+ QISKIT_IBM_INSTANCE) → 저장소 루트의 .env → --save 로 저장한 계정.
비트 순서: Qiskit 은 리틀 엔디안(오른쪽이 q0)이다. 결과 파일에는 Qiskit 순서와 이 책의 순서(왼쪽이 q0)를 둘 다 적는다.
과금: Open 플랜은 28일 이동 구간마다 QPU 10분. 대기 시간은 세지 않고 QPU 가 실제로 돈 시간만 센다. 작업이 끝나면 사용량을 찍는다.
주의: qiskit-ibm-runtime 0.50 부터 SamplerV2 는 폐기 예정(deprecated)이다(대체: qiskit_ibm_runtime.executor_sampler.Sampler).
      이 스크립트는 SamplerV2 를 쓰므로 requirements/sdk.txt 가 상한(<0.53)을 둔다. 폐기 경고가 떠도 동작은 한다(실행으로 확인하지는 못했다).
"""
import argparse
import datetime
import getpass
import json
import os
import pathlib
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parent.parent
CHANNEL = "ibm_quantum_platform"


def load_dotenv():
    """저장소 루트의 .env 를 읽어 비어 있는 환경 변수만 채운다(python-dotenv 없이)."""
    p = ROOT / ".env"
    if not p.exists():
        return
    for line in p.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        if line.startswith("export "):
            line = line[len("export "):].lstrip()
        k, v = line.split("=", 1)
        k, v = k.strip(), v.strip()
        if len(v) >= 2 and v[0] == v[-1] and v[0] in "\"'":
            v = v[1:-1]                      # 따옴표로 싼 값은 그대로(안의 # 도 값)
        else:
            v = v.split(" #", 1)[0].strip()  # 따옴표 없는 값 뒤의 주석은 버린다
        if v and not os.environ.get(k):
            os.environ[k] = v


def get_service():
    from qiskit_ibm_runtime import QiskitRuntimeService
    load_dotenv()
    token = os.environ.get("QISKIT_IBM_TOKEN")
    if token:
        kw = {"channel": CHANNEL, "token": token}
        if os.environ.get("QISKIT_IBM_INSTANCE"):
            kw["instance"] = os.environ["QISKIT_IBM_INSTANCE"]
        return QiskitRuntimeService(**kw)
    try:
        return QiskitRuntimeService()      # --save 로 저장한 기본 계정
    except Exception as e:
        sys.exit(f"IBM 계정에 연결하지 못했다({type(e).__name__}: {str(e)[:200]}).\n먼저 python env/run_ibm.py --save 를 하거나 .env 에 QISKIT_IBM_TOKEN 을 넣는다.")


def save_account():
    from qiskit_ibm_runtime import QiskitRuntimeService
    print("IBM Quantum Platform 대시보드에서 만든 API 키를 붙여 넣는다(화면에 보이지 않는다).")
    token = getpass.getpass("API 키: ").strip()
    if not token:
        sys.exit("키가 비어 있다.")
    crn = input("인스턴스 CRN (crn:v1:... 로 시작; 비우면 자동 선택): ").strip()
    kw = {"channel": CHANNEL, "token": token, "set_as_default": True, "overwrite": True}
    if crn:
        kw["instance"] = crn
    QiskitRuntimeService.save_account(**kw)
    print("저장했다(~/.qiskit/qiskit-ibm.json). 본인만 쓰는 컴퓨터에서만 저장한다. 지우려면 QiskitRuntimeService.delete_account().")


def build(name, n):
    from qiskit import QuantumCircuit
    if name == "bell":
        n = 2
    if n < 2:
        sys.exit("--n 은 2 이상이어야 한다.")
    qc = QuantumCircuit(n, name=f"{name}{n}")
    qc.h(0)
    for i in range(n - 1):
        qc.cx(i, i + 1)
    qc.measure_all()
    ideal = {"0" * n: 0.5, "1" * n: 0.5}       # Bell·GHZ 는 뒤집어도 같은 비트열
    return qc, ideal


def counts_of(pub):
    """SamplerV2 결과 하나에서 counts 를 꺼낸다(레지스터 이름은 measure_all 의 'meas')."""
    try:
        return pub.join_data().get_counts()
    except Exception:
        return pub.data.meas.get_counts()


def two_qubit_gates(circ):
    """2큐비트 게이트 수. measure_all() 이 넣는 barrier(지시문)는 게이트가 아니므로 세지 않는다."""
    return sum(1 for inst in circ.data if inst.operation.num_qubits == 2 and inst.operation.name != "barrier")


def gate_depth(circ):
    """게이트만 센 깊이(barrier 와 측정 제외) — tools/gate_cost.py 의 규약과 맞춘다."""
    return circ.depth(lambda inst: inst.operation.name not in ("barrier", "measure"))


def circuit_stats(isa):
    return {"two_qubit_gates": two_qubit_gates(isa), "depth": gate_depth(isa), "ops": {k: int(v) for k, v in isa.count_ops().items()}}


def final_usage(job, tries=6, wait=5):
    """job.usage() 는 집계가 끝나기 전에는 0 을 준다. 잠깐 기다렸다가 그래도 0 이면 None(미확정)으로 적는다."""
    for _ in range(tries):
        try:
            u = job.usage()
        except Exception:
            return None
        if u:
            return u
        time.sleep(wait)
    return None


def tvd(counts, ideal):
    shots = sum(counts.values())
    keys = set(counts) | set(ideal)
    return 0.5 * sum(abs(counts.get(k, 0) / shots - ideal.get(k, 0)) for k in keys)


def report(job_id, backend_name, qc, stats, counts, ideal, usage, shots, local, note=None):
    import qiskit
    book = {k[::-1]: v for k, v in counts.items()}
    top = sorted(counts.items(), key=lambda kv: -kv[1])[:6]
    print(f"\n장치 {backend_name} · 작업 {job_id}")
    print(f"변환 뒤 회로: 2큐비트 게이트 {stats['two_qubit_gates']}, 깊이 {stats['depth']}(측정 제외), 연산 {stats['ops']}" + (f"  [{note}]" if note else ""))
    print("가장 많이 나온 비트열(Qiskit 순서):", ", ".join(f"{k} {v}" for k, v in top))
    print(f"이상 분포와의 전변동 거리(TVD): {tvd(counts, ideal):.4f}   (0 이면 완벽. 샷 잡음만 있을 때의 평균은 약 {0.4 / shots ** 0.5:.3f})")
    if not local:
        print(f"QPU 사용량: {usage} 초   (Open 플랜 한도 600초 / 28일)" if usage is not None else "QPU 사용량: 미확정 — 집계가 늦다. 대시보드에서 확인한다.")
    now = datetime.datetime.now(datetime.timezone.utc)
    out = {
        "when_utc": now.isoformat(timespec="seconds"),
        "backend": backend_name, "local": local, "job_id": job_id, "circuit": qc.name, "shots": shots,
        **stats, "transpile_note": note,
        "counts_qiskit_order": counts, "counts_book_order": book, "tvd_vs_ideal": tvd(counts, ideal), "usage_seconds": usage,
        "versions": {"python": sys.version.split()[0], "qiskit": qiskit.__version__},
    }
    try:
        import qiskit_ibm_runtime
        out["versions"]["qiskit_ibm_runtime"] = qiskit_ibm_runtime.__version__
    except Exception:
        pass
    path = ROOT / "runs" / f"{now:%Y%m%d-%H%M%S}Z_{backend_name}_{qc.name}.json"
    path.parent.mkdir(exist_ok=True)
    path.write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
    print("저장:", path.relative_to(ROOT))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--save", action="store_true"); ap.add_argument("--list", action="store_true")
    ap.add_argument("--local", action="store_true"); ap.add_argument("--fetch", metavar="JOB_ID")
    ap.add_argument("--circuit", choices=["bell", "ghz"], default="bell"); ap.add_argument("--n", type=int, default=3)
    ap.add_argument("--backend"); ap.add_argument("--shots", type=int, default=1024)
    ap.add_argument("--yes", action="store_true", help="실제로 제출한다(무료 분량을 쓴다)")
    a = ap.parse_args()

    if a.save:
        return save_account()
    from qiskit.transpiler import generate_preset_pass_manager
    from qiskit_ibm_runtime import SamplerV2 as Sampler
    qc, ideal = build(a.circuit, a.n)

    if a.local:
        from qiskit_aer import AerSimulator
        backend = AerSimulator()
        isa = generate_preset_pass_manager(optimization_level=1, backend=backend).run(qc)
        job = Sampler(mode=backend).run([isa], shots=a.shots)
        counts = counts_of(job.result()[0])
        return report("local", "aer_simulator", qc, circuit_stats(isa), counts, ideal, None, a.shots, True)

    service = get_service()
    if a.list:
        try:
            crn = service.active_instance()          # CRN 은 로그에 그대로 남기지 않는다
            print("활성 인스턴스:", (crn[:24] + "…" + crn[-6:]) if crn and len(crn) > 34 else crn)
        except Exception:
            pass
        for b in service.backends(simulator=False, operational=True):
            st = b.status()
            print(f"{b.name:<22} {b.num_qubits:>4}큐비트  대기 {st.pending_jobs:>5}  {st.status_msg}")
        try:
            print("사용량:", service.usage())
        except Exception:
            pass
        return
    if a.fetch:
        pend = ROOT / "runs" / f"pending_{a.fetch}.json"
        meta = json.loads(pend.read_text(encoding="utf-8")) if pend.exists() else None
        if meta:
            qc, ideal = build(meta["circuit"], meta["n"])
        else:
            print("경고: runs/pending_<id>.json 이 없다. --circuit/--n 을 제출 때와 똑같이 줘야 TVD 가 맞는다.")
        job = service.job(a.fetch)
        print("상태:", job.status())
        counts = counts_of(job.result()[0])
        if len(next(iter(counts))) != qc.num_qubits:
            sys.exit(f"결과 비트 수({len(next(iter(counts)))})가 회로({qc.num_qubits}큐비트)와 다르다. 제출 때의 --circuit/--n 을 준다.")
        be = job.backend()
        if meta:
            stats, note = {k: meta[k] for k in ("two_qubit_gates", "depth", "ops")}, None
        else:                                   # 제출 때의 변환 결과를 모르면 다시 변환한 참고값을 적는다
            stats, note = circuit_stats(generate_preset_pass_manager(optimization_level=1, backend=be).run(qc)), "받을 때 다시 변환한 참고값"
        report(job.job_id(), be.name if be else "unknown", qc, stats, counts, ideal, final_usage(job), sum(counts.values()), False, note)
        if pend.exists():
            pend.unlink()
        return

    backend = service.backend(a.backend) if a.backend else service.least_busy(operational=True, simulator=False, min_num_qubits=qc.num_qubits)
    isa = generate_preset_pass_manager(optimization_level=1, backend=backend).run(qc)
    stats = circuit_stats(isa)
    print(f"제출 계획: {qc.name} · {a.shots}샷 · 장치 {backend.name}({backend.num_qubits}큐비트, 대기 {backend.status().pending_jobs}) · 2큐비트 게이트 {stats['two_qubit_gates']} · 깊이 {stats['depth']}")
    if not a.yes:
        print("--yes 가 없어 제출하지 않았다. 이 정도 회로는 보통 QPU 시간 수 초를 쓴다(작업 뒤에 실제 값을 찍는다).")
        return
    job = Sampler(mode=backend).run([isa], shots=a.shots)
    pend = ROOT / "runs" / f"pending_{job.job_id()}.json"       # 끊겨도 --fetch 가 같은 회로·같은 지표로 마무리하도록
    pend.parent.mkdir(exist_ok=True)
    pend.write_text(json.dumps({"job_id": job.job_id(), "backend": backend.name, "circuit": a.circuit, "n": qc.num_qubits,
                                "shots": a.shots, **stats}, indent=1), encoding="utf-8")
    print(f"제출했다. 작업 ID {job.job_id()} — 대기열에 따라 몇 분에서 몇 시간 걸린다. 끊겨도 --fetch {job.job_id()} 로 다시 받는다.")
    counts = counts_of(job.result()[0])
    report(job.job_id(), backend.name, qc, stats, counts, ideal, final_usage(job), a.shots, False)
    pend.unlink()


if __name__ == "__main__":
    main()
