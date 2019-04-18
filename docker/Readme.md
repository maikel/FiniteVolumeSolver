# Using Docker with FiniteVolumeSolver

We have several docker images to test building the code  with various compilers.

Currently supported Compilers are 

  * gcc 7 and 8
  * LLVM Clang 5, 6, 7 and 8
  * AppleClang from in Xcode 10

If you just want to run a simulation or try something with a certain compiler you can simply run 
one of the docker images which are stored in the gitlab. 

To run a docker container and build this FiniteVolumeSolver repository there you have first to 
login with your docker into this gitlab instance

```bash 
docker login git.imp.fu-berlin.de:5000
```

and then you can use images which are stored in this repository

```bash
git clone git@git.imp.fu-berlin.de:ag-klein/FiniteVolumeSolver.git FiniteVolumeSolver/
docker run -it -v FiniteVolumeSolver/:/FiniteVolumeSolver git.imp.fu-berlin.de:5000/ag-klein/finitevolumesolver/amrex:2_clang5
root@da4ed06a478c:/# mkdir build/
root@da4ed06a478c:/build/# cd build/
root@da4ed06a478c:/build/# cmake /FiniteVolumeSolver -DCMAKE_BUILD_TYPE=Release
root@da4ed06a478c:/build/# make AMReX.EB.Ramp
root@da4ed06a478c:/build/# ./examples/AMReX.EB.Ramp
```
