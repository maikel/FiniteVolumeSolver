FROM debian:testing

RUN apt-get update && apt-get install -y git gcc g++ gfortran python3 python3-pip libopenmpi-dev cmake pkg-config && pip3 install conan
RUN conan profile new default --detect
RUN conan profile update settings.compiler.libcxx=libstdc++11 default
RUN conan remote add gitlab https://git.imp.fu-berlin.de/api/v4/packages/conan
RUN conan user guttula -r gitlab -p VwArmxLep3iHxEsshJEU
RUN git clone https://oauth2:LxnsQNy9iN6i6uNA2Dxo@git.imp.fu-berlin.de/ag-klein/FiniteVolumeSolver.git --branch feature/AxialConservative
RUN mkdir /FiniteVolumeSolver/build
WORKDIR /FiniteVolumeSolver/build
RUN conan install ../ -o AMReX:dim=2 --build missing
RUN cmake ../ -DCMAKE_BUILD_TYPE=Release
RUN make AMReX.EB.Divider2D

FROM debian:testing
RUN apt-get update
RUN apt-get install -y fish python3 libopenmpi3
RUN apt-get install -y libpython3.9 libgomp1
RUN useradd amrex -d /home/amrex
RUN mkdir /home/amrex
RUN chown amrex:amrex -R /home/amrex
RUN chsh -s `which fish` amrex
COPY --from=0 /FiniteVolumeSolver/build/examples/AMReX.EB.Divider2D /home/amrex/
COPY --from=0 /FiniteVolumeSolver/examples/AMReX/EB/2D/Divider /home/amrex/
COPY --from=0 /FiniteVolumeSolver/extra/Divider/* /home/amrex/
COPY --from=0 /FiniteVolumeSolver/examples/AMReX/EB/2D/Divider2D.py /home/amrex/
RUN chown amrex:amrex -R /home/amrex/*
RUN apt-get install -y vim 
RUN apt-get install -y openmpi-bin libopenmpi-dev
USER amrex
WORKDIR /home/amrex