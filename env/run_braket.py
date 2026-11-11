"""run_braket.py — AWS Braket 장치(또는 로컬 시뮬레이터)에서 Bell·GHZ 회로를 돌린다.
작성됨, 미실행 — 이 스크립트를 쓴 환경에는 Braket SDK 가 없어 문법 검사만 했다. 첫 실행은 --local 로 한다.

  python env/run_braket.py --local                         계정 없이 LocalSimulator
  python env/run_braket.py --list                          지금 온라인인 QPU 와 ARN
  python env/run_braket.py --device <ARN> --shots 100      비용만 계산해 보여 준다
  python env/run_braket.py --device <ARN> --shots 100 --yes --max-usd 2
  python env/run_braket.py --fetch <task ARN>              전에 제출한 작업의 결과

인증: aws configure 로 넣은 자격 증명(또는 AWS_ACCESS_KEY_ID 등 환경 변수). Braket 을 처음 쓰면 콘솔에서 서비스 사용 동의를 한 번 해야 한다.
비용: 작업당 $0.30 + 샷당 요금(아래 표는 2026-11-11 에 aws.amazon.com/braket/pricing 에서 읽은 값). QPU 무료 분량은 없다.
비트 순서: Braket 의 비트열은 왼쪽이 q0 다(이 책과 같다).
"""
import argparse
import datetime
import json
import os
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
PER_TASK = 0.30
PER_SHOT = {  # ARN 에 들어 있는 이름 조각 → 샷당 달러. 게이트 회로를 받는 장치만(QuEra Aquila 는 아날로그 장치라 뺐다)
    "Cepheus": 0.000425, "Garnet": 0.00145, "Emerald": 0.0016, "IBEX": 0.0235, "Forte": 0.08,
}


def load_region():
    """저장소 루트 .env 의 AWS_DEFAULT_REGION 을 읽는다(없으면 us-east-1). 자격 증명은 aws configure 의 것을 쓴다."""
    p = ROOT / ".env"
    if p.exists() and not os.environ.get("AWS_DEFAULT_REGION"):
        for line in p.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if line.startswith("AWS_DEFAULT_REGION="):
                v = line.split("=", 1)[1].split(" #", 1)[0].strip().strip("\"'")
                if v:
                    os.environ["AWS_DEFAULT_REGION"] = v
    os.environ.setdefault("AWS_DEFAULT_REGION", "us-east-1")


def build(name, n):
    from braket.circuits import Circuit
    if name == "bell":
        n = 2
    c = Circuit().h(0)
    for i in range(n - 1):
        c.cnot(i, i + 1)
    return c, {"0" * n: 0.5, "1" * n: 0.5}, f"{name}{n}"


def tvd(counts, ideal):
    shots = sum(counts.values())
    return 0.5 * sum(abs(counts.get(k, 0) / shots - ideal.get(k, 0)) for k in set(counts) | set(ideal))


def save(tag, device, label, shots, counts, ideal, cost):
    out = {"when": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"), "device": device, "task": tag,
           "circuit": label, "shots": shots, "counts_book_order": counts, "tvd_vs_ideal": tvd(counts, ideal), "cost_usd": cost,
           "versions": {"python": sys.version.split()[0]}}
    try:
        from importlib.metadata import version
        out["versions"]["amazon-braket-sdk"] = version("amazon-braket-sdk")
    except Exception:
        pass
    path = ROOT / "runs" / f"{datetime.datetime.now(datetime.timezone.utc):%Y%m%d-%H%M%S}Z_braket_{label}.json"
    path.parent.mkdir(exist_ok=True)
    path.write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"가장 많이 나온 비트열: {sorted(counts.items(), key=lambda kv: -kv[1])[:6]}")
    print(f"이상 분포와의 TVD: {out['tvd_vs_ideal']:.4f} · 저장: {path.relative_to(ROOT)}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--local", action="store_true"); ap.add_argument("--list", action="store_true")
    ap.add_argument("--device", metavar="ARN"); ap.add_argument("--fetch", metavar="TASK_ARN")
    ap.add_argument("--circuit", choices=["bell", "ghz"], default="bell"); ap.add_argument("--n", type=int, default=3)
    ap.add_argument("--shots", type=int, default=100); ap.add_argument("--max-usd", type=float, default=2.0)
    ap.add_argument("--yes", action="store_true")
    a = ap.parse_args()
    if a.shots < 1:
        sys.exit("--shots 는 1 이상이어야 한다.")
    if a.circuit == "ghz" and a.n < 2:
        sys.exit("--n 은 2 이상이어야 한다.")
    circ, ideal, label = build(a.circuit, a.n)

    if a.local:
        from braket.devices import LocalSimulator
        r = LocalSimulator().run(circ, shots=a.shots).result()
        return save("local", "LocalSimulator", label, a.shots, dict(r.measurement_counts), ideal, 0.0)
    load_region()
    from braket.aws import AwsDevice, AwsQuantumTask
    if a.list:
        from braket.aws.aws_device import AwsDeviceType
        for d in AwsDevice.get_devices(types=[AwsDeviceType.QPU], statuses=["ONLINE"]):
            print(f"{d.provider_name:<10} {d.name:<22} {d.arn}")
        return
    if a.fetch:
        t = AwsQuantumTask(arn=a.fetch)
        state = t.state()
        print("상태:", state)
        if state in ("FAILED", "CANCELLED"):
            sys.exit(f"작업이 {state} 로 끝났다. 콘솔에서 사유를 확인한다.")
        r = t.result()
        if r is None:
            sys.exit("결과를 받지 못했다. 나중에 다시 --fetch 한다.")
        counts = dict(r.measurement_counts)
        n_bits = len(next(iter(counts)))
        ideal = {"0" * n_bits: 0.5, "1" * n_bits: 0.5}            # Bell·GHZ 전용 스크립트이므로 비트 수로 정한다
        label = "bell2" if n_bits == 2 else f"ghz{n_bits}"
        try:
            device = t.metadata().get("deviceArn", "?")
        except Exception:
            device = "?"
        return save(a.fetch, device, label, sum(counts.values()), counts, ideal, None)
    if not a.device:
        sys.exit("--device <ARN> 이 필요하다(--list 로 확인).")
    rate = next((v for k, v in PER_SHOT.items() if k.lower() in a.device.lower()), None)
    if rate is None:
        sys.exit("이 장치의 샷당 요금을 모른다. aws.amazon.com/braket/pricing 에서 확인해 PER_SHOT 표에 더한다.")
    cost = PER_TASK + a.shots * rate
    print(f"제출 계획: {label} · {a.shots}샷 · {a.device}\n예상 비용: $0.30 + {a.shots} × ${rate} = ${cost:.2f}  (한도 --max-usd {a.max_usd})")
    if cost > a.max_usd:
        sys.exit("예상 비용이 한도를 넘는다. 샷 수를 줄이거나 --max-usd 를 올린다.")
    if not a.yes:
        print("--yes 가 없어 제출하지 않았다.")
        return
    dev = AwsDevice(a.device)
    print("장치 상태:", dev.status)                       # OFFLINE·RETIRED 면 제출해도 돌지 않는다
    task = dev.run(circ, shots=a.shots)
    print(f"제출했다. 작업 ARN {task.id}\n장치가 열려 있는 시간대에만 실행되므로 오래 걸릴 수 있다. 끊기면 --fetch 로 다시 받는다.")
    r = task.result()
    save(task.id, a.device, label, a.shots, dict(r.measurement_counts), ideal, cost)


if __name__ == "__main__":
    main()
