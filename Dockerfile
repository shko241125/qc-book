# 책 실행 환경 이미지. 빌드: docker build -t qc-book .
# 실행: docker run --rm -it -v "$PWD":/work $([ -f .env ] && echo --env-file .env) qc-book   (.env 가 없으면 --env-file 을 빼야 한다)
# 컨테이너는 root 로 돈다 — Linux 에서는 runs/ 등에 쓴 파일이 root 소유가 된다(필요하면 --user "$(id -u):$(id -g)").
# 인터넷이 없는 곳으로 옮길 때: docker save qc-book | gzip > qc-book-image.tar.gz  →  docker load < qc-book-image.tar.gz
FROM python:3.12-slim
ENV PIP_NO_CACHE_DIR=1 MPLBACKEND=Agg PYTHONUNBUFFERED=1
WORKDIR /work
COPY requirements/ /tmp/requirements/
RUN pip install --upgrade pip && pip install -r /tmp/requirements/sdk.txt
# AWS Braket 경로까지 쓰려면 아래 줄의 주석을 푼다
# RUN pip install -r /tmp/requirements/hardware.txt
CMD ["bash"]
