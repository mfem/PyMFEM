# ==============================================================================
# For alternative base images, see https://github.com/mfem/containers/tree/main 
# cpu-tpls includes OpenMPI, hypre, SuperLU, and PETSc
FROM ghcr.io/mfem/containers/cpu-tpls:latest

# ==============================================================================
# PyMFEM (will separate this into another layer later)
RUN git clone https://github.com/mfem/PyMFEM.git
ENV PYTHONPATH "${PYTHONPATH}:/home/euler/PyMFEM"
RUN python3 -m pip install jupyterlab

# ==============================================================================
# Install GLVis dependencies
USER root
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y \
    freeglut3-dev libglew-dev libsdl2-dev libfreetype-dev \
    libfontconfig1-dev libglm-dev libpng-dev xxd

# Git clone all relevant glvis repos
RUN git clone https://github.com/GLVis/glvis.git && \
    git clone https://github.com/GLVis/glvis-js.git && \
    git clone https://github.com/GLVis/pyglvis.git

# GLVis js / jupyter extension
#RUN git lfs install && git clone https://github.com/GLVis/glvis-js.git

# Build GLVis
WORKDIR glvis
RUN make

# Switch back to user
USER euler
WORKDIR /home/euler
# ==============================================================================
