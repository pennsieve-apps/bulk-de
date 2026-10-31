# Gene sets: downloaded and checked against pinned sha256 at build time (no network at run time).
FROM python:3.12.16-slim-bookworm@sha256:e697c81897302a4c4333bade000d98038fbb21eb064f8ee1fd42d02ebcd9d337 AS genesets
COPY genesets/build.py /build.py
RUN python /build.py /opt/genesets

FROM python:3.12.16-slim-bookworm@sha256:e697c81897302a4c4333bade000d98038fbb21eb064f8ee1fd42d02ebcd9d337

ENV PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

WORKDIR /processor

# Pinned Python dependencies (wheels only), plus the Lambda runtime interface client.
COPY requirements.txt /processor/requirements.txt
RUN pip install --timeout 120 --retries 10 --only-binary=:all: -r /processor/requirements.txt \
    && find /usr/local/lib/python3.12/site-packages -name '*.pyc' -delete \
    && find /usr/local/lib/python3.12/site-packages -type d -name tests -prune -exec rm -rf {} +

COPY --from=genesets /opt/genesets /opt/genesets
COPY processor/ /processor/processor/
COPY VERSION /processor/VERSION

ENV PYTHONPATH="/processor" \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    MPLBACKEND=Agg \
    MPLCONFIGDIR=/tmp/matplotlib \
    GENESETS_DIR=/opt/genesets

# Build matplotlib's font cache now, writable by the run-time user (Lambda's file system is
# read-only outside /tmp), and check that the stack imports and openpyxl uses defusedxml.
RUN python -c "import matplotlib.pyplot, pydeseq2.dds, gseapy, awslambdaric; from openpyxl.xml import DEFUSEDXML; assert DEFUSEDXML" \
    && chmod -R a+rwX /tmp/matplotlib

# Runs as uid 1000 on the platform; nothing here needs root.
USER 1000

# Runs main() once as a container, or serves processor.handler as a Lambda function (see processor/entry.py).
ENTRYPOINT ["python3", "-m", "processor.entry"]
