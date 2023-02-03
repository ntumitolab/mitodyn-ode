FROM ghcr.io/sosiristseng/docker-jupyterbook:0.13.1.7

# Julia
ENV JULIA_CI true
ENV JULIA_NUM_THREADS "auto"
ENV JULIA_PATH /usr/local/julia/
ENV JULIA_DEPOT_PATH /srv/juliapkg/
ENV PATH ${JULIA_PATH}/bin:${PATH}
COPY --from=julia:1.8.5 ${JULIA_PATH} ${JULIA_PATH}

# Python dependencies. e.g. matplotlib
RUN pip install --no-cache-dir matplotlib

# Julia environment
COPY Project.toml Manifest.toml .github/install_kernel.jl ./
COPY src/ src
RUN julia --color=yes --project="" install_kernel.jl && \
    julia --color=yes --project=@. -e 'import Pkg; Pkg.instantiate(); Pkg.resolve(); Pkg.precompile()'

CMD ["jupyter-book"]
