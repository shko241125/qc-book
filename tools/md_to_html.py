"""Build chNN_slug.html from a hand-written chNN_slug.md (figure embedded as data URI).
usage: python md_to_html.py <slug> <title> <png or -> [extra png ...]   (run inside the chapter directory)
"""
import base64, re, sys, pathlib, markdown

slug, title, png = sys.argv[1:4]
md = pathlib.Path(f"{slug}.md").read_text(encoding="utf-8")
md_html = re.sub(r"```latex\n(.*?)\n```", lambda m: "\n$$\n" + m.group(1) + "\n$$\n", md, flags=re.S)
for _png in [png] + sys.argv[4:]:          # 추가 그림(예: 도식 스크린샷)도 같은 방식으로 인라인
    if _png == "-": continue
    b = base64.b64encode(pathlib.Path(_png).read_bytes()).decode()
    md_html = md_html.replace(f"({_png})", f"(data:image/png;base64,{b})")
body = markdown.markdown(md_html, extensions=["tables", "fenced_code"])
# Python-Markdown 은 \~ 와 (표 밖의) \| 를 이스케이프로 보지 않아 역슬래시가 그대로 남는다. 코드·수식 밖에서만 지운다.
_parts = re.split(r"(<pre>.*?</pre>|<code>.*?</code>|\$\$.*?\$\$)", body, flags=re.S)
body = "".join(x if i % 2 else x.replace("\\~", "~").replace("\\|", "|") for i, x in enumerate(_parts))
html = f"""<!doctype html><html lang="ko"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{title}</title>
<script>window.MathJax={{tex:{{displayMath:[['$$','$$']]}}}};</script>
<script defer src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>
<style>
body{{max-width:760px;margin:40px auto;padding:0 16px;font-family:-apple-system,"Noto Sans KR","Apple SD Gothic Neo",sans-serif;line-height:1.65;color:#1f2328}}
pre{{background:#f6f8fa;padding:12px;border-radius:8px;overflow-x:auto;font-size:13px}}
code{{font-family:ui-monospace,Menlo,monospace}} table{{border-collapse:collapse;width:100%;margin:12px 0}}
th,td{{border:1px solid #d0d7de;padding:6px 10px;text-align:left}} th{{background:#f6f8fa}}
img{{max-width:100%}} h2{{margin-top:2em;border-bottom:1px solid #d0d7de;padding-bottom:4px}}
</style></head><body>{body}</body></html>"""
pathlib.Path(f"{slug}.html").write_text(html, encoding="utf-8")
print(f"{slug}.md {len(md)} chars; {slug}.html {len(html)} chars")
