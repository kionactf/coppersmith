FROM sagemath/sagemath:latest

RUN sudo apt-get update && sudo apt-get install -y tzdata  # avoid select timezone
RUN sudo apt-get update && sudo apt-get upgrade -y

RUN sudo apt-get install -y \
    vim \
    less \
    git \
    tmux \
    netcat

# for fplll, flatter
RUN sudo apt-get update \
    && sudo apt-get install -y \
    cmake \
    libtool \
    fplll-tools \
    libfplll-dev \
    libgmp-dev \
    libmpfr-dev \
    libeigen3-dev \
    libblas-dev \
    liblapack-dev

USER sage

RUN sage --pip install pycryptodome pwntools tqdm


RUN mkdir /home/sage/coppersmith
COPY --chown=sage:sage *.py /home/sage/coppersmith/
COPY --chown=sage:sage fplll /home/sage/coppersmith/fplll/
COPY --chown=sage:sage flatter /home/sage/coppersmith/flatter/

#fplll
WORKDIR /home/sage/coppersmith/fplll
RUN ./autogen.sh
RUN ./configure
RUN make -j4

#flatter
WORKDIR /home/sage/coppersmith/flatter
RUN mkdir build
WORKDIR /home/sage/coppersmith/flatter/build
RUN cmake ..
RUN make -j4

WORKDIR /home/sage/


## other lattice library download
RUN mkdir lattice_tools
WORKDIR /home/sage/lattice_tools

# defund/coppersmith
RUN git clone https://github.com/defund/coppersmith
RUN ln -s /home/sage/coppersmith/coppersmith.sage /home/sage/defund_coppersmith.sage
## load("/home/sage/lattice_tools/defund_coppersmith.sage")

# josephsurin/lattice-based-cryptanalysis
RUN git clone https://github.com/josephsurin/lattice-based-cryptanalysis
RUN sed -i "s/algorithm='msolve', //g" /home/sage/lattice_tools/lattice-based-cryptanalysis/lbc_toolkit/common/systems_solvers.sage
ENV PYTHONPATH=/home/sage/lattice_tools/lattice-based-cryptanalysis/:$PYTHONPATH

# jvdsn/crypto-attacks
RUN git clone https://github.com/jvdsn/crypto-attacks/
RUN mv crypto-attacks crypto_attacks
ENV PYTHONPATH=/home/sage/lattice_tools/crypto_attacks/:$PYTHONPATH

# rkm0959/Inequality_Solving_with_CVP
RUN git clone https://github.com/rkm0959/Inequality_Solving_with_CVP
RUN ln -s /home/sage/lattice_tools/Inequality_Solving_with_CVP/solver.sage /home/sage/lattice_tools/inequ_cvp_solve.sage
## load("/home/sage/lattice_tools/inequ_cvp_solve.sage

WORKDIR /home/sage


CMD ["/home/sage/sage/sage", "--python", "coppersmith/example.py"]
