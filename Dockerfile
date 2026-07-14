FROM continuumio/miniconda3:25.1.1-2
LABEL version="1.3.4" \
      description="A pipeline for viral metagenomics analysis."

WORKDIR /app

COPY environment.yml /app/environment.yml
# Install the base environment (ViralUnity + core runtime deps).
RUN conda env create --quiet -f environment.yml && conda clean -a -y

ENV PATH=/opt/conda/envs/viralunity/bin:$PATH

COPY . /app/viralunity

WORKDIR /app/viralunity
RUN pip install . && rm -rf /root/.cache/pip
RUN viralunity --version > /app/viralunity-version.txt

# Run as an unprivileged user rather than root.
RUN useradd --create-home --uid 1000 viralunity \
    && chown -R viralunity:viralunity /app
USER viralunity

WORKDIR /tmp/
ENTRYPOINT ["viralunity"]
CMD ["--help"]

# NOTE (follow-up, see docs/embedding.md and CONTAINER_IMAGE_PLAN.md):
# per-rule conda envs are still created on first run. To pre-build them into
# the image, run `viralunity setup --pipelines all --conda-prefix <fixed path>`
# here and make the pipeline reuse that same --conda-prefix at runtime. This
# needs a fixed, world-readable prefix (not $HOME-derived) so the non-root
# runtime user can reuse the root-built envs; it requires build-time network
# and iterative testing, so it is intentionally left as a separate change.
