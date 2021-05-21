FROM debian:testing

RUN apt-get update && apt-get install -y git gcc g++ gfortran python3 python3-pip libopenmpi-dev cmake pkg-config && pip3 install conan
RUN conan profile new default --detect
RUN conan profile update settings.compiler.libcxx=libstdc++11 default
RUN conan remote add gitlab https://git.imp.fu-berlin.de/api/v4/packages/conan
RUN conan user zenkechr -r gitlab -p J4sqpLdyBaEhZxH89spJ
RUN git clone https://oauth2:LxnsQNy9iN6i6uNA2Dxo@git.imp.fu-berlin.de/ag-klein/FiniteVolumeSolver.git --branch Docker/Divider2d
RUN mkdir /FiniteVolumeSolver/build
WORKDIR /FiniteVolumeSolver/build
RUN conan install ../ -o AMReX:dim=2 --build missing
ARG INCUBATOR_VER=unknown
RUN git pull
RUN cmake ../ -DCMAKE_BUILD_TYPE=Release
RUN make AMReX.EB.Divider2D -j4

FROM debian:testing
RUN apt-get update
RUN apt-get install -y fish gcc g++ gfortran python3 python3-pip libopenmpi3 libpython3.9 libgomp1 vim openmpi-bin
RUN pip3 install h5py 
RUN pip3 install matplotlib==3.3.3 --force-reinstall
RUN pip3 install cython
RUN pip3 install yt==3.6.1
RUN useradd amrex -d /home/amrex
RUN mkdir /home/amrex
RUN chown amrex:amrex -R /home/amrex
RUN chsh -s `which fish` amrex
COPY --from=0 /FiniteVolumeSolver/build/examples/AMReX.EB.Divider2D /home/amrex/
COPY amrex/ /home/amrex/amrex/
COPY Divider /home/amrex/
COPY VisualizeHdf5.py /home/amrex/
ARG INCUBATOR_VER2=unknown
COPY Divider2D.py /home/amrex/
RUN chown amrex:amrex -R /home/amrex/*
USER amrex
WORKDIR /home/amrex
ENV OMPI_MCA_btl_vader_single_copy_mechanism=none
ENV OMP_NUM_THREADS=1
