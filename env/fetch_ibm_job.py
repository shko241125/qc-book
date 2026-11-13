"""fetch_ibm_job.py — 끝난 IBM 작업에서 '실제로 돈 회로'와 그 시각의 보정값을 받아 JSON 으로 남긴다(QPU 를 쓰지 않는다).

  python env/fetch_ibm_job.py <job id> [<job id> ...] --out runs/<이름>.json

남기는 것: 작업마다 장치·생성 시각·PUB 별 게이트 목록(이름, 물리 큐비트, 각도) — 트랜스파일러가 만든 ISA 회로 그대로.
그리고 작업 생성 시각 직전의 보정값(backend.properties(datetime=...)): 쓰인 큐비트의 T1·T2·판독 오류,
쓰인 sx·x·cz 의 오류율. numpy 스크립트(부록 A)는 이 파일만 읽으므로 네트워크 없이 다시 돌릴 수 있다.
작업 입력의 구조(inputs['quantum_program'].items[i].circuit)는 qiskit-ibm-runtime 0.52 의 클라이언트 쪽 프리미티브 기준이다.
"""
import argparse
import json
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "env"))
from run_ibm import get_service  # noqa: E402  (.env → 저장 계정 순서의 인증을 그대로 쓴다)


def num(p):
    try:
        return float(p)
    except TypeError:                     # Estimator 작업의 회로는 기저 회전 각을 매개변수로 남긴다(값은 작업 안에 따로 있다)
        return str(p)


def gate_list(circ):
    idx = lambda q: circ.find_bit(q).index
    return [[i.operation.name, [idx(q) for q in i.qubits], [num(p) for p in i.operation.params],
             [circ.find_bit(c).index for c in i.clbits]]
            for i in circ.data if i.operation.name != "barrier"]


def calib(props, qubits, gates):
    q = {}
    for k in sorted(qubits):
        def qp(name):
            try:
                return props.qubit_property(k, name)[0]
            except Exception:
                return None
        q[k] = {"T1_us": props.t1(k) * 1e6, "T2_us": props.t2(k) * 1e6, "readout_error": props.readout_error(k),
                "p_meas1_prep0": qp("prob_meas1_prep0"), "p_meas0_prep1": qp("prob_meas0_prep1")}
    g = {}
    for name, qs in sorted({(n, tuple(qs)) for n, qs, _, _ in gates if n in ("sx", "x", "cz")}):
        try:
            g[f"{name}{list(qs)}"] = props.gate_error(name, list(qs))
        except Exception:
            g[f"{name}{list(qs)}"] = None
    return {"last_update_date": str(props.last_update_date), "qubits": q, "gates": g}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("job_ids", nargs="+"); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    service = get_service()
    out = {}
    for jid in a.job_ids:
        job = service.job(jid)
        be = job.backend()
        items = job.inputs["quantum_program"].items
        pubs = [gate_list(it.circuit) for it in items]
        used = {q for p in pubs for _, qs, _, _ in p for q in qs}
        props = be.properties(datetime=job.creation_date)
        out[jid] = {"backend": be.name, "created": str(job.creation_date), "status": str(job.status()),
                    "pubs": pubs, "calibration": calib(props, used, [g for p in pubs for g in p])}
        data = job.result()[0].data
        if hasattr(data, "evs"):                 # Estimator: 기댓값과 표준오차(목표 precision 이 아니라 실제 값)
            out[jid]["estimator"] = {"evs": float(data.evs), "stds": float(data.stds)}
        print(f"{jid}: {be.name}, PUB {len(pubs)}개, 큐비트 {sorted(used)}, 보정 시각 {props.last_update_date}")
    p = ROOT / a.out
    p.write_text(json.dumps(out, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    print("저장:", p.relative_to(ROOT))


if __name__ == "__main__":
    main()
