FROM sagemath/sagemath:latest

RUN sudo apt-get update && sudo apt-get install -y tzdata  # avoid select timezone
RUN sudo apt-get update && sudo apt-get upgrade -y

RUN sudo apt-get install -y \
    vim \
    less \
    git \
    tmux \
    netcat

### for building fplll, flatter
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

RUN sage --pip install --no-cache-dir \
    pwntools \
    pycryptodome \
    z3-solver \
    tqdm

### prepare build fplll, flatter
RUN mkdir /home/sage/coppersmith
COPY --chown=sage:sage *.py /home/sage/coppersmith/
COPY --chown=sage:sage fplll /home/sage/coppersmith/fplll/
COPY --chown=sage:sage flatter /home/sage/coppersmith/flatter/

### build fplll
WORKDIR /home/sage/coppersmith/fplll
RUN ./autogen.sh
RUN ./configure
RUN make -j4

### build flatter
WORKDIR /home/sage/coppersmith/flatter
RUN mkdir build
WORKDIR /home/sage/coppersmith/flatter/build
RUN cmake ..
RUN make -j4


#### other lattice library download
RUN mkdir /home/sage/collection_lattice_tools
WORKDIR /home/sage/collection_lattice_tools
ENV PYTHONPATH=/home/sage/collection_lattice_tools/:$PYTHONPATH

### defund/coppersmith
RUN git clone https://github.com/defund/coppersmith
RUN ln -s /home/sage/collection_lattice_tools/coppersmith/coppersmith.sage /home/sage/collection_lattice_tools/defund_coppersmith.sage
## load("/home/sage/collection_lattice_tools/defund_coppersmith.sage")

### josephsurin/lattice-based-cryptanalysis
RUN git clone https://github.com/josephsurin/lattice-based-cryptanalysis
RUN sed -i "s/algorithm='msolve', //g" /home/sage/collection_lattice_tools/lattice-based-cryptanalysis/lbc_toolkit/common/systems_solvers.sage
ENV PYTHONPATH=/home/sage/collection_lattice_tools/lattice-based-cryptanalysis/:$PYTHONPATH

### jvdsn/crypto-attacks
RUN git clone https://github.com/jvdsn/crypto-attacks/
RUN mv crypto-attacks crypto_attacks
ENV PYTHONPATH=/home/sage/collection_lattice_tools/crypto_attacks/:$PYTHONPATH

### rkm0959/Inequality_Solving_with_CVP
RUN git clone https://github.com/rkm0959/Inequality_Solving_with_CVP
RUN ln -s /home/sage/collection_lattice_tools/Inequality_Solving_with_CVP/solver.sage /home/sage/collection_lattice_tools/inequ_cvp_solve.sage
## load("/home/sage/collection_lattice_tools/inequ_cvp_solve.sage")

### nneonneo/pwn-stuff (including math/solvelinmod.py)
RUN git clone https://github.com/nneonneo/pwn-stuff
RUN ln -s /home/sage/collection_lattice_tools/pwn-stuff/math/solvelinmod.py /home/sage/collection_lattice_tools/solvelinmod.py


### coppersmith package config
WORKDIR /home/sage/
ENV PYTHONPATH=/home/sage/coppersmith/:$PYTHONPATH
ENV COPPERSMITHFPLLLPATH=/home/sage/coppersmith/fplll/fplll/
ENV COPPERSMITHFLATTERPATH=/home/sage/coppersmith/flatter/build/bin/


### follow cryptohack-docker (https://github.com/cryptohack/cryptohack-docker)
ENV PWNLIB_NOTERM=true

ENV BLUE='\033[0;34m'
ENV YELLOW='\033[1;33m'
ENV RED='\e[91m'
ENV NOCOLOR='\033[0m'

ENV BANNER=" \n\
${BLUE}----------------------------------${NOCOLOR} \n\
${YELLOW}CryptoHack Docker Container${NOCOLOR} \n\
\n\
${RED}After Jupyter starts, visit http://127.0.0.1:8888${NOCOLOR} \n\
${BLUE}----------------------------------${NOCOLOR} \n\
"

#CMD ["echo -e $BANNER && sage -n jupyter --NotebookApp.token='' --no-browser --ip='0.0.0.0' --port=8888"]
CMD ["/home/sage/sage/sage", "--python", "coppersmith/example.py"]
