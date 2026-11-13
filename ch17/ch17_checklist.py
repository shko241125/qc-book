"""17장 발견 절차 — 3절의 도식 「발견 절차 · 5단계와 탈락 신호」(ch17_checklist.png)를 그린다.
실행: python ch17_checklist.py [출력 파일]   (인자가 없으면 이 파일 옆의 ch17_checklist.png)
계산은 없다. 상자의 글은 본문 3절의 체크리스트와 5절 후보 표에서 온 것이고, 고치려면 아래 STEPS 를 고친다.
좌표는 가로 1344 × 세로 996 의 격자(위가 0)이고 dpi 150 으로 2016 × 1494 픽셀을 만든다.
"""
import os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

matplotlib.rcParams["font.family"] = ["Noto Sans CJK JP", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False
matplotlib.rcParams["mathtext.fontset"] = "dejavusans"            # √N 의 근호만 수식 글꼴로 그린다(CJK 글꼴의 √ 는 전각이라 N 과 벌어진다)

TITLE = "다섯 단계 중 하나라도 걸리면 탈락한다"
# (단계 이름, 묻는 것) → 아니오일 때 (탈락 신호, 예·설명)
STEPS = [
    (("1. 문제 정의", "입력·출력·정확도가 정해졌나"), ("출력이 벡터 전체, 입력 적재 O(N)", "예: 피치 추정을 QFT로")),
    (("2. 구조 찾기", "주기성·국소성·희소성이 있나"), (r"구조 없음 — 최선이 $\sqrt{N}$", "예: 구조 없는 SAT 탐색")),
    (("3. 원시 요소 대응", "16장 부품으로 조립되나"), ("대응하는 부품 없음", "새 부품이 필요(18장 탐색 도구)")),
    (("4. 자원 추정", "Toffoli·논리 큐비트·초기 상태 겹침"), ("자원이 기술 전망 밖", "큐비트·실행 시간·겹침 γ")),
    (("5. 고전 베이스라인", "같은 지표에서 손익분기를 넘나"), ("고전이 같은 지표에서 더 빠름", "예: 1D 바닥상태, 몬테카를로 가격")),
]
PASS = ("가능 또는 조건부", "가정과 함께 19장 후보로")

W, H = 1344, 996
LX0, LX1 = 70, 601                    # 왼쪽(단계) 상자의 가로 범위
RX0, RX1 = 707, 1273                  # 오른쪽(탈락 신호) 상자의 가로 범위
TOPS = [113, 254, 396, 537, 679]      # 단계 상자의 위쪽 y
BOX_H = 100
PASS_TOP = 834
INK, SUB, LINE = "#1F1F1F", "#5F5F5F", "#C4C4B8"
RED, RED_FILL = "#D9453F", "#F9EDEC"
BLUE, BLUE_FILL = "#2B7BD9", "#E5EEF8"


def box(ax, x0, x1, top, edge, fill, lw):
    ax.add_patch(FancyBboxPatch((x0, top), x1 - x0, BOX_H, boxstyle="round,pad=0,rounding_size=12",
                                linewidth=lw, edgecolor=edge, facecolor=fill, zorder=2))


def two_lines(ax, xc, top, head, sub, bold):
    ax.text(xc, top + 35, head, ha="center", va="center", fontsize=16, color=INK, fontweight="bold" if bold else "normal", zorder=3)
    ax.text(xc, top + 71, sub, ha="center", va="center", fontsize=13.5, color=SUB, zorder=3)


def arrow(ax, x0, y0, x1, y1):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0), zorder=1,
                arrowprops=dict(arrowstyle="-|>", color=LINE, lw=1.6, mutation_scale=16, shrinkA=0, shrinkB=0))


fig = plt.figure(figsize=(W / 100, H / 100), dpi=150)
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(H, 0); ax.axis("off")
ax.text(42, 50, TITLE, ha="left", va="center", fontsize=19, color="#111111", fontweight="bold")

lxc, rxc = (LX0 + LX1) / 2, (RX0 + RX1) / 2
for i, ((head, sub), (rhead, rsub)) in enumerate(STEPS):
    top = TOPS[i]; yc = top + BOX_H / 2
    box(ax, LX0, LX1, top, LINE, "white", 1.6); two_lines(ax, lxc, top, head, sub, bold=True)
    box(ax, RX0, RX1, top, RED, RED_FILL, 1.6); two_lines(ax, rxc, top, rhead, rsub, bold=False)
    arrow(ax, LX1, yc, RX0, yc)                                              # 아니오 → 탈락 신호
    nxt = TOPS[i + 1] if i + 1 < len(TOPS) else PASS_TOP
    arrow(ax, lxc, top + BOX_H, lxc, nxt)                                    # 예 → 다음 단계
ax.text((LX1 + RX0) / 2, TOPS[0] + 28, "아니오", ha="center", va="center", fontsize=13.5, color=SUB)
ax.text(lxc + 27, (TOPS[0] + BOX_H + TOPS[1]) / 2, "예", ha="center", va="center", fontsize=13.5, color=SUB)
box(ax, LX0, LX1, PASS_TOP, BLUE, BLUE_FILL, 2.4); two_lines(ax, lxc, PASS_TOP, PASS[0], PASS[1], bold=True)

out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "ch17_checklist.png")
fig.savefig(out, dpi=150, facecolor="white")
print(f"saved {os.path.basename(out)} ({len(STEPS)}단계, 탈락 신호 {len(STEPS)}개, 통과 상자 1개)")
