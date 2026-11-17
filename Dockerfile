FROM --platform=linux/amd64 python:3.12-slim

WORKDIR /app

# NCBI command-line tools (asn2gb, asn2fsa) for egapx2mss, in /app/bin (its default bin dir).
# These binaries expire periodically; rebuild with --no-cache to refresh.
RUN apt-get update && apt-get install -y --no-install-recommends wget \
 && wget -q http://archive.debian.org/debian/pool/main/libi/libidn/libidn11_1.33-3_amd64.deb \
 && dpkg -i libidn11_1.33-3_amd64.deb \
 && rm libidn11_1.33-3_amd64.deb \
 && mkdir -p bin \
 && wget -qO- https://ftp.ncbi.nih.gov/toolbox/ncbi_tools/cmdline/asn2gb.linux64.gz  | gunzip > bin/asn2gb \
 && wget -qO- https://ftp.ncbi.nih.gov/toolbox/ncbi_tools/cmdline/asn2fsa.linux64.gz | gunzip > bin/asn2fsa \
 && chmod +x bin/asn2gb bin/asn2fsa \
 && apt-get purge -y wget && apt-get autoremove -y \
 && rm -rf /var/lib/apt/lists/*

# Third-party dependencies first (cached on source-only changes).
COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

# Source only (no examples/, tests/). Editable install keeps /app/bin as egapx2mss's default bin dir.
COPY pyproject.toml README.md LICENSE ./
COPY src/ src/
RUN pip install --no-cache-dir -e .

CMD ["bash"]
