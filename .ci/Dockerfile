FROM python:3.11.0-slim

# Supress warning emitted by parallel
ENV LC_ALL C.UTF-8

# System packages
# RUN apt-get update && apt-get install -y ffmpeg fonts-humor-sans --no-install-recommends && rm -rf /var/lib/apt/lists/*

# Julia
ENV JULIA_PKG_SERVER ""
ENV JULIA_PATH /usr/local/julia/
ENV PATH $JULIA_PATH/bin:$PATH
COPY --from=julia:1.8.2 ${JULIA_PATH} /usr/local/

# Python dependencies e.g. matplotlib
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Julia environment
COPY Project.toml Manifest.toml ./
# Uncomment if you have src folder
COPY src/ src
RUN julia --threads=auto --color=yes --project="" -e 'import Pkg; Pkg.add("IJulia"); Pkg.build("IJulia")' && \
    julia --threads=auto --color=yes --project=@. -e 'import Pkg; Pkg.instantiate()'

CMD ["jupyter-book"]
