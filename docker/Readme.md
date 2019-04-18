# Using Docker with FiniteVolumeSolver

We have several docker images to test and build the code with various compilers.

Currently supported Compilers are 

  * gcc 7 and 8
  * LLVM Clang 5, 6, 7 and 8
  * AppleClang from in Xcode 10

If you just want to run a simulation or try something with a certain compiler 
you other would not have access to, you can simply run one of the docker images 
which are stored here.

You have to login to this gitlab instance to run a docker container and build 
the source code of this project. You can login via the command

```bash 
docker login git.imp.fu-berlin.de:5000
```

and your usual gitlab credentials. Afterwards you can use images which are 
stored in this repository, for example

```bash
git clone git@git.imp.fu-berlin.de:ag-klein/FiniteVolumeSolver.git FiniteVolumeSolver/
docker run -it -v FiniteVolumeSolver/:/FiniteVolumeSolver git.imp.fu-berlin.de:5000/ag-klein/finitevolumesolver/amrex:2d_gcc8
root@da4ed06a478c:/# mkdir build/
root@da4ed06a478c:/build/# cd build/
root@da4ed06a478c:/build/# cmake /FiniteVolumeSolver -DCMAKE_BUILD_TYPE=Release
root@da4ed06a478c:/build/# make AMReX.EB.Ramp
root@da4ed06a478c:/build/# ./examples/AMReX.EB.Ramp
```
A list of available images, with links, can be found in the [Container Registry](https://git.imp.fu-berlin.de/ag-klein/FiniteVolumeSolver/container_registry).

## Building docker images on your own

You do not have to fetch the docker image from this docker repository since you can build 
them easily with the provided Dockerfiles in this git repository.

To build an image locally you can invoke docker via

```
docker build -t MyImageTag -f amrex_3d-gcc8 .
```

and then simply run `MyImageTag`

```
docker run -it MyImageTag
```