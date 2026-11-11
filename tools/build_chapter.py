"""Build chNN_slug.md and chNN_slug.html from an exported doc (base64 markdown).
usage: python build_chapter.py <b64file> <slug> <title> <png or ->  [byline]
run inside the chapter directory.
"""
import base64, re, sys, pathlib, markdown

b64file, slug, title, png = sys.argv[1:5]
byline = sys.argv[5] if len(sys.argv) > 5 else "Gosh"
raw = pathlib.Path(b64file).read_text().strip()
raw += "=" * (-len(raw) % 4)
md = base64.b64decode(raw).decode("utf-8", errors="replace")
if "�" in md:
    i = md.find("�")
    sys.exit(f"base64 transcription corrupted near char {i}: ...{md[max(0,i-80):i]}")
MONTHS = {m: i for i, m in enumerate(["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"], 1)}
def _byline(m):
    return f"{m.group(3)}-{MONTHS[m.group(1)]:02d}-{int(m.group(2)):02d} · {byline}"
md = re.sub(r"^(\w{3}) (\d{1,2}), (\d{4}) · @\w+$", _byline, md, count=1, flags=re.M)

if png != "-":
    fig_md = f"\n![결과 그림]({png})\n\n"
    md = md.replace("\n## 6. ", fig_md + "## 6. ", 1)
pathlib.Path(f"{slug}.md").write_text(md, encoding="utf-8")

md_html = re.sub(r"```latex\n(.*?)\n```", lambda m: "\n$$\n" + m.group(1) + "\n$$\n", md, flags=re.S)
if png != "-":
    b = base64.b64encode(pathlib.Path(png).read_bytes()).decode()
    md_html = md_html.replace(f"({png})", f"(data:image/png;base64,{b})")
body = markdown.markdown(md_html, extensions=["tables", "fenced_code"])
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
