"""부록 A 5절 — BV '0000' 확정 실험의 첫 작업 제출 (qiskit · qiskit-ibm-runtime). 이 판은 작성됨, 미실행이다(SDK 가 없는 환경에서 썼다).

  python appA/appA_bv_submit.py                       (= local) 계정 없이 Aer(잡음 없음)로 회로의 출력을 확인한다. 장치에 닿지 않는다. 먼저 이것부터.
  python appA/appA_bv_submit.py real                  계획만: 장치 조회, 지금의 보정값(init_error 포함), reset 지원 여부, 예상 QPU 시간. 제출하지 않는다
  python appA/appA_bv_submit.py real --yes            실제로 제출한다 — 사용자가 이 실행을 승인했을 때만(PUB 2개 × 4000샷, QPU 약 4초)
  python appA/appA_bv_submit.py real --resume <기록>  제출 뒤 결과를 받기 전에 끊겼을 때: 기록 파일의 작업 id 로 결과만 다시 받는다(QPU 사용 없음)
인자 규칙은 ch14/ch14_submit.py 와 같다(인자 없이 돌면 장치에 닿지 않는다 — env/test_all.sh --sdk 가 이 모드로 돌린다).

회로와 판정 규칙은 appA_bv_confirm.py 가 만든 appA/appA_bv_plan.json 에서 읽는다(먼저 `python appA/appA_bv_confirm.py`).
첫 작업은 회로 둘(P0 = 2026-11-11 의 BV 그대로, P1 = 보조가 |0⟩ 으로 끝나는 변형)뿐이다. 한 작업의 회로들이 어떤 순서로 도는지
(같은 회로가 연달아인지, 샷마다 번갈아인지) 모르기 때문이다 — 회로 둘이면 어느 쪽이든 결과가 갈린다(부록 A 5절).
회로는 그날 실제로 돈 ISA 게이트 목록을 같은 물리 큐비트에 그대로 올린 것이라 트랜스파일하지 않는다. Sampler 옵션은 기본값 그대로 둔다(그날과 같게).
결과는 runs/<UTC 시각>_<장치>_bv_confirm_record.json 에 제출 직후부터 쓴다(이미 있는 이름이면 멈춘다). counts 와 샷별 비트열을 함께 남긴다.
분석: python appA/appA_bv_confirm.py runs/<그 파일>.json
계정: 환경 변수 → .env → 저장 계정(env/run_ibm.py 의 get_service). 키·CRN 은 인자나 출력에 쓰지 않는다.
"""
import json, pathlib, sys, time
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "env"))
err = lambda e: f"{type(e).__name__}: {str(e)[:200]}"          # 기록에 남기는 오류 문구는 줄인다


def build(gates, n_qubits, mp=None):
    """게이트 목록 [[이름, [큐비트], [매개변수], [고전 비트]], …] → QuantumCircuit. mp: 큐비트 번호 바꾸기."""
    mp = mp or {}
    nc = 1 + max(g[3][0] for g in gates if g[0] == "measure")
    qc = QuantumCircuit(QuantumRegister(n_qubits, "q"), ClassicalRegister(nc, "c"))
    for name, qs, params, cb in gates:
        q = [mp.get(x, x) for x in qs]
        if name == "rz": qc.rz(params[0], q[0])
        elif name == "sx": qc.sx(q[0])
        elif name == "x": qc.x(q[0])
        elif name == "cz": qc.cz(q[0], q[1])
        elif name == "reset": qc.reset(q[0])
        elif name == "measure": qc.measure(q[0], cb[0])
        else: raise ValueError(f"모르는 게이트 {name}")
    return qc


def calib_now(backend, qubits):
    """지금 조회되는 보정값(값과 그 항목의 측정 시각). init_error 가 있으면 함께."""
    out = {}
    try:
        props = backend.properties(); out["last_update_date"] = str(props.last_update_date)
        for q in qubits:
            row = {}
            for name in ("T1", "T2", "readout_error", "prob_meas1_prep0", "prob_meas0_prep1", "init_error"):
                try:
                    v, t = props.qubit_property(q, name); row[name] = [float(v), str(t)]
                except Exception:
                    row[name] = None
            out[str(q)] = row
    except Exception as e:
        out["error"] = err(e)
    return out


def collect(job, record, names, dump):
    """결과를 받아 counts 와 샷별 비트열을 기록에 넣는다."""
    res = job.result()
    for k, pub in zip(names, res):
        c = pub.data.c.get_counts(); record["pubs"][k]["counts"] = c
        try:                                           # 샷별 자료: 실행 순서와 몰림을 보려면 이것이 있어야 한다(배열 순서 = 실행 순서인지는 미확인)
            record["pubs"][k]["bitstrings"] = list(pub.data.c.get_bitstrings())
        except Exception as e:
            record["pubs"][k]["bitstrings"] = None; record["pubs"][k]["bitstrings_error"] = err(e)
        print(k, sorted(c.items(), key=lambda kv: -kv[1])[:4])
    dump()
    try:
        u = job.metrics().get("usage", {})
        secs = u.get("qpu_charge_time_seconds", u.get("quantum_seconds"))
        record["qpu_seconds"] = secs if secs is not None else job.usage()
    except Exception as e:
        record["qpu_seconds"] = f"조회 실패: {err(e)}"
    dump()
    print("qpu_seconds:", record["qpu_seconds"])


def main():
    args = sys.argv[1:]
    plan_path = HERE / "appA_bv_plan.json"
    if not plan_path.exists():
        sys.exit("appA/appA_bv_plan.json 이 없다. 먼저 python appA/appA_bv_confirm.py")
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    order, shots, roles = plan["order"], plan["shots"], plan["roles"]
    anc, data, idle = roles["ancilla"], roles["data_cz"], roles["idle_measured"]
    used = [anc] + data

    if "real" not in args:                             # local(기본값). 계정 없이: 큐비트 번호를 0…k−1 로 줄여 Aer 로 출력만 본다
        from qiskit_aer import AerSimulator
        allq = sorted(set(used + idle)); small = {q: i for i, q in enumerate(allq)}
        p3 = {q: len(allq) + i for i, q in enumerate(used)}                    # 번호만 바꾼 BV(둘째 작업의 P3)도 같은 출력인지 본다
        sim = AerSimulator(); want = plan["good_data_bits"]
        circs = {k: (plan["variants"][k], small) for k in order}
        circs.update({k: (g, small) for k, g in plan["job2_candidates"].items() if isinstance(g, list)})
        circs["P3(번호만 바꿈)"] = (plan["variants"]["P0"], {**small, **p3})
        for k, (g, mp) in circs.items():
            qc = build(g, len(allq) + len(used), mp)
            c = sim.run(qc, shots=shots).result().get_counts()
            tgt = "1111" if k == "PM" else want                                # PM 은 114번에 표지 X 가 있어 그 비트가 1 이다
            good = sum(v for s, v in c.items() if s[-4:] == tgt) / shots
            print(f"  {k:14s} ops={dict(qc.count_ops())} → {dict(sorted(c.items(), key=lambda kv: -kv[1])[:3])} | 데이터 비트 {tgt} 의 비율 {good:.4f}")
            assert good == 1.0, f"{k}: 잡음 없는 출력이 {tgt} 가 아니다"
        print("  모든 회로의 잡음 없는 출력이 기대와 같다.")
        return

    from qiskit_ibm_runtime import Sampler
    from run_ibm import get_service
    service = get_service()

    if "--resume" in args:                             # 이미 낸 작업의 결과만 다시 받는다
        rp = pathlib.Path(args[args.index("--resume") + 1]); record = json.loads(rp.read_text(encoding="utf-8"))
        if all(record["pubs"][k].get("counts") for k in record["order"]):
            sys.exit("이 기록에는 결과가 이미 있다. runs/ 의 기록은 고치지 않는다")
        job = service.job(record["job_id"])
        collect(job, record, record["order"], lambda: rp.write_text(json.dumps(record, indent=1, ensure_ascii=False) + "\n", encoding="utf-8"))
        print(f"다음: python appA/appA_bv_confirm.py {rp}")
        return

    backend = service.backend(plan["backend"])         # 같은 장치(ibm_fez)의 같은 물리 큐비트여야 뜻이 있다
    print("backend:", backend.name, "| 큐비트", backend.num_qubits)
    tgt = backend.target
    flips = []
    for k in order:                                    # CZ 는 대칭이다. target 에 반대 방향만 있으면 큐비트 순서를 바꿔 쓴다
        for g in plan["variants"][k]:
            if g[0] == "cz" and not tgt.instruction_supported("cz", tuple(g[1])):
                if tgt.instruction_supported("cz", tuple(g[1][::-1])): g[1] = g[1][::-1]; flips.append(tuple(g[1]))
                else: sys.exit(f"cz{g[1]} 가 이 장치의 target 에 어느 방향으로도 없다 — 계획을 다시 봐야 한다")
    try:
        has_reset = bool(tgt.instruction_supported("reset", (anc,)))
    except Exception:
        has_reset = "reset" in tgt.operation_names
    print(f"cz 방향을 바꾼 것 {flips or '없음'} | reset 지원(둘째 작업의 P2 용): {has_reset} | 지원 연산 {sorted(tgt.operation_names)}")
    circs = [build(plan["variants"][k], backend.num_qubits) for k in order]
    cal = calib_now(backend, sorted(set(used + idle)))
    row = cal.get(str(anc)) or {}
    print(f"지금의 보정값(보조 {anc}): " + ", ".join(f"{k} {v[0]:.4g} ({v[1]})" if v else f"{k} 없음" for k, v in row.items()))
    print(f"\n제출 계획: {backend.name} · Sampler PUB {len(circs)}개 {order} × {shots}샷. 예상 QPU 약 {plan['qpu_seconds_estimate']:.0f}초(2026-11-11 의 기록에 맞춘 어림).")
    if "--yes" not in args:
        print("제출하지 않았다. 제출하려면 real --yes (사용자 승인 후).")
        return

    stamp = time.strftime("%Y%m%d-%H%M%SZ", time.gmtime())
    keep = ROOT / "runs" / f"{stamp}_{backend.name}_bv_confirm_record.json"
    if keep.exists(): sys.exit(f"{keep.name} 이 이미 있다")
    sampler = Sampler(mode=backend)                    # 옵션은 기본값 그대로(2026-11-11 과 같다: 트월링·DD 꺼짐, init_qubits·rep_delay 기본값)
    record = {"experiment": "appA 5절 BV 0000 확정 실험 — 첫 작업 [P0, P1]", "backend": backend.name, "shots": shots,
              "t_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "order": order,
              "good_data_bits": plan["good_data_bits"], "bit_order": "qiskit(c_{m-1}…c0)", "rules": plan["rules"],
              "primitive": f"{type(sampler).__module__}.{type(sampler).__name__}", "calibration_at_submit": cal, "cz_flipped": flips,
              "pubs": {k: {"ops": {a: int(b) for a, b in c.count_ops().items()}} for k, c in zip(order, circs)}}
    for key, fn in (("rep_delay_default", lambda: float(backend.default_rep_delay)), ("sampler_options", lambda: repr(sampler.options)[:2000])):
        try:
            record[key] = fn()
        except Exception as e:
            record[key] = f"조회 실패: {err(e)}"
    dump = lambda: keep.write_text(json.dumps(record, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    job = sampler.run(circs, shots=shots)              # 회로 둘을 한 작업으로 — 같은 시각, 같은 보정 상태
    record["job_id"] = job.job_id(); dump()
    print("job id:", record["job_id"], "| 기록:", keep.name, flush=True)
    print(f"(결과를 받기 전에 끊기면: python appA/appA_bv_submit.py real --resume runs/{keep.name})", flush=True)
    collect(job, record, order, dump)
    print(f"다음: python appA/appA_bv_confirm.py runs/{keep.name}")
    print(f"      python env/fetch_ibm_job.py {record['job_id']} --out runs/{stamp}_{backend.name}_bv_confirm_isa_calib.json   # 그 시각의 보정값(QPU 사용 없음)")


if __name__ == "__main__":
    main()
