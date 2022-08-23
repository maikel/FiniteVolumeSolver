import os
from conans import ConanFile, CMake, tools


class Hdf5Conan(ConanFile):
    name = "HDF5"
    version = "1.10"
    license = "https://support.hdfgroup.org/ftp/HDF5/releases/COPYING"
    url = "<Package recipe repository url here, for issues about the package>"
    description = "High-performance data management and storage suite"
    settings = "os", "compiler", "build_type", "arch"
    options = {"shared": [True, False]}
    default_options = "shared=False"
    generators = "cmake"

    def requirements(self):
        self.requires.add("zlib/[>=1.2.11]")
        # Use system mpi 
        # self.requires.add("OpenMPI/3.1.3@finite-volume/stable")

    def source(self):
        self.run("git clone https://github.com/HDFGroup/hdf5.git  --branch hdf5_1_10 --single-branch --depth=1")
        # Garantee proper linkage
        tools.replace_in_file("hdf5/CMakeLists.txt", "project (HDF5 C)",
                              '''project (HDF5 C)
include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup()''')

    def _configure_cmake(self):
        cmake = CMake(self)
        cmake_args = {
            "BUILD_SHARED_LIBS": "ON" if self.options.shared else "OFF",
            "BUILD_TESTING": "OFF",
            "BUILD_EXAMPLES": "OFF",
            "HDF5_BUILD_CPP_LIB": "OFF",
            "HDF5_BUILD_TOOLS": "ON",
            "HDF_ENABLE_PARALLEL": "ON",
        }
        cmake.configure(source_folder="hdf5", defs=cmake_args)
        return cmake

    def build(self):
        cmake = self._configure_cmake()
        cmake.build()

    def package(self):
        cmake = self._configure_cmake()
        cmake.install()

    def package_info(self):
        self.cpp_info.libs = tools.collect_libs(self)

        if not self.options.shared and not self.settings.os == "Windows":
            self.cpp_info.libs.append("dl")

