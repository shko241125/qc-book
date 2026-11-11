"""extract_snippets.py — 각 장 md 의 SDK 코드 블록을 실행 가능한 파일로 뽑아낸다.
실행: python env/extract_snippets.py            → snippets/chNN/chNN_KK_<sdk>.py 와 snippets/index.json
책의 SDK 스니펫은 "작성됨, 미실행"으로 표시돼 있다. 이 도구로 뽑아서 env/run_snippets.py 로 돌리면
어느 스니펫이 실제 SDK에서 통과하는지 한 번에 알 수 있다. SDK 없이도 추출은 된다.

분류: 블록 안의 import 로 SDK를 판정한다. SDK import 가 없는 블록(numpy 발췌 등)은 건너뛴다.
"""
import json
import pathlib
import re
import shutil

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = ROOT / "snippets"
SDK_PATTERNS = [  # 먼저 맞는 것이 파일 이름에 쓰인다(여러 개면 전부 기록)
    ("cudaq", r"^\s*(import|from)\s+cudaq\b"),
    ("braket", r"^\s*(import|from)\s+braket\b"),
    ("qualtran", r"^\s*(import|from)\s+qualtran\b"),
    ("pennylane", r"^\s*(import|from)\s+pennylane\b"),
    ("cirq", r"^\s*(import|from)\s+cirq\b"),
    ("qutip", r"^\s*(import|from)\s+qutip(_qip)?\b"),
    ("qiskit", r"^\s*(import|from)\s+qiskit(_aer|_ibm_runtime)?\b"),
]
HARDWARE = re.compile(r"QiskitRuntimeService\s*\(|AwsDevice\s*\(|AwsQuantumTask\s*\(|least_busy\s*\(|IonQProvider\s*\(|qiskit\.remote|braket\.aws\.qubit")
FENCE = re.compile(r"^```(\w*)\s*$")


def blocks(md_text):
    """(시작 줄 번호, 언어, 코드) 를 차례로 낸다."""
    lines = md_text.split("\n"); i = 0
    while i < len(lines):
        m = FENCE.match(lines[i])
        if m and m.group(1):
            lang, start, j = m.group(1), i + 1, i + 1
            while j < len(lines) and not lines[j].startswith("```"): j += 1
            yield start + 1, lang, "\n".join(lines[start:j])
            i = j + 1
        else:
            i += 1


def main():
    shutil.rmtree(OUT, ignore_errors=True)        # 지난번 추출물은 지우고 다시 뽑는다(보고서는 runs/ 에 따로 남는다)
    index = []
    for md in sorted(ROOT.glob("ch*/ch*.md")):
        ch = md.parent.name; k = 0
        for line_no, lang, code in blocks(md.read_text(encoding="utf-8")):
            if lang != "python": continue
            sdks = [name for name, pat in SDK_PATTERNS if re.search(pat, code, re.M)]
            if not sdks: continue
            k += 1
            dest = OUT / ch / f"{ch}_{k:02d}_{sdks[0]}.py"
            dest.parent.mkdir(parents=True, exist_ok=True)
            header = (f"# 출처: {md.relative_to(ROOT)} {line_no}행의 코드 블록 (env/extract_snippets.py 가 뽑음)\n"
                      f"# SDK: {', '.join(sdks)}\n")
            dest.write_text(header + code.rstrip() + "\n", encoding="utf-8")
            index.append({"file": str(dest.relative_to(ROOT)), "chapter": ch, "md": str(md.relative_to(ROOT)), "line": line_no,
                          "sdks": sdks, "hardware": bool(HARDWARE.search(code)),
                          "marked_unexecuted": "미실행" in code.split("\n", 1)[0], "lines": code.count("\n") + 1})
    OUT.mkdir(exist_ok=True)
    (OUT / "index.json").write_text(json.dumps(index, ensure_ascii=False, indent=1), encoding="utf-8")
    by = {}
    for e in index: by[e["sdks"][0]] = by.get(e["sdks"][0], 0) + 1
    print(f"스니펫 {len(index)}개를 {OUT.relative_to(ROOT)}/ 에 뽑았다:", ", ".join(f"{s} {n}" for s, n in sorted(by.items())))
    print(f"실기 호출이 들어 있는 것: {sum(e['hardware'] for e in index)}개 (run_snippets.py 는 --hardware 없이는 건너뛴다)")


if __name__ == "__main__":
    main()
