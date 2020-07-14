from conans import ConanFile, CMake, tools
from datetime import datetime


class VcConan(ConanFile):
    name = "Vc"
    version = "1.4.1"
    license = "BSD 3-Clause New or Revised License"
    url = "https://github.com/VcDevel/Vc"
    description = "SIMD Vector Classes for C++"
    settings = "os", "compiler", "build_type", "arch"
    options = {"shared": [True, False]}
    default_options = {"shared": False}
    generators = "cmake"

    def source(self):
        # self.run("git clone https://github.com/VcDevel/Vc.git --branch 1.4 --single-branch --depth=1")
        git = tools.Git(folder="Vc")
        git.clone("https://github.com/VcDevel/Vc.git", "1.4")
        # Garantee proper linkage
        tools.replace_in_file("Vc/CMakeLists.txt", "set(CMAKE_MODULE_PATH \"${CMAKE_SOURCE_DIR}/cmake\")",
                              '''set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
include("${CMAKE_BINARY_DIR}/conanbuildinfo.cmake")
conan_basic_setup()''')


    def _configure_cmake(self):
        cmake = CMake(self)
        cmake.configure(source_folder="Vc", defs={
            "BUILD_SHARED_LIBS": "ON" if self.options.shared else "OFF",
            "BUILD_TESTING": "OFF"
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


