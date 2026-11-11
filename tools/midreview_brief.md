# Mid-project review brief (after ch17 of 27) — common context for all reviewers

You are an independent reviewer of a Korean quantum-computing textbook written chapter by chapter.
You did NOT write it. Your job is to find things that are actually wrong, not to praise or restyle.

## Where things are
- Book root: `/home/claude/qc-book/`. One directory per chapter: `ch01` … `ch17`.
  Each has `chNN_<slug>.md` (the chapter text, source of truth), `chNN_<slug>.py` (numpy verification
  script; run with `MPLBACKEND=Agg python3 chNN_<slug>.py` inside the chapter directory), a results png, and an html build.
- Fresh stdout of every script (re-run today) is in
  `/tmp/claude-0/-home-claude/fd63146d-511d-5ae9-8dab-3100fb85f0db/scratchpad/mid/chNN.out`
  (if a file is missing or empty, run the script yourself; ch16 takes ~2 min, ch15 ~1 min, others seconds).
- Every chapter has the same 10 sections: 1 요약 / 2 직관(Quirk) / 3 정의·수식+교과서 비교 / 4 회로 구성 /
  5 시뮬레이션 결과 / 6 다른 SDK / 7 한계와 반례 / 8 연습문제와 해답 / 9 학습 노트 / 10 출처와 환경.

## Conventions you must know (do not report these as errors)
- Qubit order is big-endian in all formulas and numpy code (matches Nielsen & Chuang and Cirq); Qiskit snippets convert.
- SDK snippets (Qiskit, Cirq, PennyLane as `import pennylane as qp`, CUDA-Q, QuTiP) are marked "작성됨, 미실행"
  because PyPI is unreachable in this sandbox. You cannot run them either (`pip install` fails) — review them by reading
  and by checking official docs on the web if you have web tools.
- Text is Korean; code, API names and technical terms stay in English.
- Chapters were already reviewed 1–3 times for per-chapter numeric agreement with the scripts. Do not repeat a
  line-by-line number check unless your task says so; spend your effort on your assigned angle.

## Rules
- Do NOT edit anything under `/home/claude/qc-book/`. Write scratch code only under
  `/tmp/claude-0/-home-claude/fd63146d-511d-5ae9-8dab-3100fb85f0db/scratchpad/mid/<your-id>/`.
- Verify before you report: recompute with numpy/scipy, or cite the primary source you checked (URL).
  A finding without evidence must be labelled "unverified suspicion".
- Work and write your report in English. When you propose replacement text for the book, write that text in Korean,
  literally (no unicode escapes), ready to paste.
- Output: one markdown file at the path given in your task. For each finding give:
  `[severity] chNN §section — file:line` / quote of the current text (short, exact) / what is wrong and the evidence /
  proposed replacement text. Severity: **A** = factually or numerically wrong, **B** = misleading/incomplete or
  contradicts another chapter, **C** = minor (typo, wording, style rule). Sort by severity. End with a 3–6 line summary and
  an explicit list of what you checked and found correct, and what you could not check.
- Your final message back should be a short summary (counts by severity, the 5 most important findings, the report path).
