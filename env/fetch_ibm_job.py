"""fetch_ibm_job.py — 끝난 IBM 작업에서 '실제로 돈 회로'와 그 시각의 보정값을 받아 JSON 으로 남긴다(QPU 를 쓰지 않는다).

  python env/fetch_ibm_job.py <job id> [<job id> ...] --out runs/<이름>.json [--bitstrings]

남기는 것: 작업마다 장치·생성 시각·PUB 별 게이트 목록(이름, 물리 큐비트, 각도) — 트랜스파일러가 만든 ISA 회로 그대로.
그리고 작업 생성 시각 직전의 보정값(backend.properties(datetime=...)): 쓰인 큐비트의 T1·T2·판독 오류,
쓰인 sx·x·cz 의 오류율. numpy 스크립트(부록 A)는 이 파일만 읽으므로 네트워크 없이 다시 돌릴 수 있다.
작업 입력의 구조(inputs['quantum_program'].items[i].circuit)는 qiskit-ibm-runtime 0.52 의 클라이언트 쪽 프리미티브 기준이다.
2026-11-11 에 실행한 것은 이전 판이다. 그 뒤 더한 것(항목별 측정 시각, init_error, 게이트 길이, Sampler counts, 끝나지 않은 작업 건너뛰기,
이미 있는 파일 거절 — 검토 R11)은 작성됨, 미실행이다. 추가 항목은 받지 못하면 None 으로 남고 기존 항목은 그대로다.
--bitstrings(검토 R14, 작성됨·미실행): Sampler 작업의 샷별 비트열도 남긴다(PUB 마다 샷 수만큼의 문자열). 오답이 한때에 몰렸는지 볼 때 쓴다.
배열의 순서가 실행 순서라는 보장은 확인하지 못했다.
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
    """쓰인 큐비트·게이트의 보정값. 항목마다 측정 시각이 따로 있으므로(last_update_date 는 마지막 보정 갱신 시각일 뿐이다) 함께 남긴다.
    init_error 는 직전 실험이 |1⟩ 이었을 때 초기화 뒤에 남는 |1⟩ 의 비율 — prob_meas1_prep0 과 다른 양이다(부록 A 5절)."""
    q = {}
    for k in sorted(qubits):
        def qp(name):
            try:
                v, t = props.qubit_property(k, name)       # (값, 측정 시각)
                return v, str(t)
            except Exception:
                return None, None
        vals = {n: qp(n) for n in ("T1", "T2", "readout_error", "prob_meas1_prep0", "prob_meas0_prep1", "init_error", "readout_length")}
        q[k] = {"T1_us": props.t1(k) * 1e6, "T2_us": props.t2(k) * 1e6, "readout_error": props.readout_error(k),
                "p_meas1_prep0": vals["prob_meas1_prep0"][0], "p_meas0_prep1": vals["prob_meas0_prep1"][0],
                "init_error": vals["init_error"][0], "readout_length_s": vals["readout_length"][0],
                "measured_at": {n: t for n, (_, t) in vals.items()}}
    g, glen = {}, {}
    for name, qs in sorted({(n, tuple(qs)) for n, qs, _, _ in gates if n in ("sx", "x", "cz")}):
        key = f"{name}{list(qs)}"
        try:
            g[key] = props.gate_error(name, list(qs))
        except Exception:
            g[key] = None
        try:
            glen[key] = props.gate_length(name, list(qs))  # 초. 부록 A 의 재현은 지금 14장의 가정(32·68 ns)을 쓴다
        except Exception:
            glen[key] = None
    return {"last_update_date": str(props.last_update_date), "qubits": q, "gates": g, "gate_length_s": glen}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("job_ids", nargs="+"); ap.add_argument("--out", required=True); ap.add_argument("--bitstrings", action="store_true")
    a = ap.parse_args()
    p = (ROOT / a.out).resolve()
    if p.exists():
        sys.exit(f"{a.out} 이 이미 있다. runs/ 의 기록은 고치지 않는다 — 다른 이름을 준다.")
    service = get_service()
    out = {}
    for jid in a.job_ids:
        job = service.job(jid)
        be = job.backend()
        try:
            items = job.inputs["quantum_program"].items
        except (KeyError, AttributeError):
            sys.exit(f"{jid}: 작업 입력에 quantum_program 이 없다 — 클라이언트 쪽 프리미티브(runtime 0.50 이상)로 낸 작업만 읽는다.")
        pubs = [gate_list(it.circuit) for it in items]
        used = {q for pb in pubs for _, qs, _, _ in pb for q in qs}
        props = be.properties(datetime=job.creation_date)
        status = str(job.status())
        out[jid] = {"backend": be.name, "created": str(job.creation_date), "status": status,
                    "pubs": pubs, "calibration": calib(props, used, [g for pb in pubs for g in pb])}
        if "DONE" not in status.upper():         # 끝나지 않은 작업의 result() 는 끝날 때까지 막고, 실패한 작업은 예외를 낸다
            print(f"{jid}: 상태 {status} — 회로와 보정값만 받고 결과는 받지 않았다")
            continue
        res = job.result()
        data = res[0].data
        if hasattr(data, "evs"):                 # Estimator: 기댓값과 표준오차(목표 precision 이 아니라 결과와 함께 온 값)
            out[jid]["estimator"] = {"evs": float(data.evs), "stds": float(data.stds)}
        else:                                    # Sampler: PUB 별 counts 도 남긴다(실행 기록과 대조용)
            try:
                out[jid]["counts"] = [pub.join_data().get_counts() for pub in res]
            except Exception as e:
                out[jid]["counts"] = f"받지 못함: {type(e).__name__}: {e}"
            if a.bitstrings:
                try:
                    out[jid]["bitstrings"] = [list(pub.join_data().get_bitstrings()) for pub in res]
                except Exception as e:
                    out[jid]["bitstrings"] = f"받지 못함: {type(e).__name__}: {str(e)[:200]}"
        print(f"{jid}: {be.name}, PUB {len(pubs)}개, 큐비트 {sorted(used)}, 마지막 보정 갱신 {props.last_update_date}")
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(out, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    print("저장:", p)


if __name__ == "__main__":
    main()
