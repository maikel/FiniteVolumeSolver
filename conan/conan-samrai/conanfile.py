from conans import ConanFile, CMake, tools


class SamraiConan(ConanFile):
    name = "SAMRAI"
    version = "3.11.1"
    license = "https://raw.githubusercontent.com/LLNL/SAMRAI/feature/blt/LICENSE"
    url = "https://git.imp.fu-berlin.de/denzler/conanpackages"
    description = "Structured Adaptive Mesh Refinement Application Infrastructure - a scalable C++ framework for block-structured AMR application development"
    settings = "os", "compiler", "build_type", "arch"
    options = {"shared": [True, False]}
    default_options = "shared=False" # set openmpi fortran options if openMPI is build: , "OpenMPI:fortran=mpifh"
    generators = "cmake"

    def requirements(self):
        self.requires.add("HDF5/1.10@finite-volume/stable")
        # Use system mpi
        # self.requires.add("OpenMPI/3.1.3@finite-volume/stable")

    def source(self):
        self.run("git clone https://github.com/LLNL/SAMRAI.git --branch feature/blt --single-branch --depth=1")
        self.run("cd SAMRAI && git submodule update --init --recursive")
        # Garantee proper linkage
        tools.replace_in_file("SAMRAI/CMakeLists.txt", "project(SAMRAI C CXX Fortran)",
                              '''project(SAMRAI C CXX Fortran)
include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup()''')

    def _configure_cmake(self):
        cmake = CMake(self)
        cmake.configure(source_folder="SAMRAI", defs={
            "HDF5_USE_STATIC_LIBRARIES": "ON",
            "BUILD_SHARED_LIBS": "ON" if self.options.shared else "OFF",
            "ENABLE_EXAMPLES": "OFF",
            "ENABLE_TESTS": "OFF",
            "ENABLE_OPENMP": "OFF",
            "ENABLE_VALGRIND": "OFF",
        })

        return cmake

    def build(self):
        cmake = self._configure_cmake()
        cmake.build()

    def package(self):
        cmake = self._configure_cmake()
        cmake.install()

    def package_info(self):
        self.cpp_info.libs = tools.collect_libs(self)

