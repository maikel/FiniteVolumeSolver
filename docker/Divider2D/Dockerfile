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