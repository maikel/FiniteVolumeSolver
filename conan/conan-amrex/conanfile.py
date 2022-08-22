from conans import ConanFile, CMake, tools
from datetime import datetime


class AmrexConan(ConanFile):
    name = "AMReX"
    version = datetime.today().strftime('%y.%m.%d')
    license = "https://raw.githubusercontent.com/AMReX-Codes/amrex/master/license.txt"
    url = "https://github.com/AMReX-Codes/amrex"
    description = "A software framework for massively parallel, block-structured adaptive mesh refinement (AMR) applications"
    settings = "os", "compiler", "build_type", "arch"
    options = {"shared": [True, False],
               "eb": [True, False],
               "dim": [1, 2, 3],
               "mpi": [True, False],
               "omp": [True, False]}
    default_options = {"mpi": True,
                       "shared": False,
                       "eb": True,
                       "dim": 3,
                       "omp": True}

    generators = "cmake"

    # def requirements(self):
        # Use system MPI
        # self.requires.add("OpenMPI/[>=3.0]@finite-volume/stable")

    def source(self):
        self.run("git clone https://github.com/AMReX-Codes/amrex.git --branch development --single-branch --depth=1")
        # Garantee proper linkage
        tools.replace_in_file("amrex/CMakeLists.txt", "project(AMReX)",
                              '''project (AMReX)
include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup()''')


    def _configure_cmake(self):
        cmake = CMake(self)
        cmake.configure(source_folder="amrex", defs={
            "BUILD_SHARED_LIBS": "ON" if self.options.shared else "OFF",
            "ENABLE_MPI": "ON" if self.options.mpi else "OFF",
            "DIM": self.options.dim,
            "ENABLE_EB": "ON" if self.options.eb else "OFF",
            "ENABLE_OMP": "ON" if self.options.omp else "OFF"
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
