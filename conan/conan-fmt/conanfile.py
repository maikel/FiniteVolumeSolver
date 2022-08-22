from conans import ConanFile, CMake, tools
from datetime import datetime


class AmrexConan(ConanFile):
    name = "fmt"
    version = "9.0.0" # datetime.today().strftime('%y.%m.%d')
    # license = ""
    url = "https://github.com/fmtlib/fmt.git"
    description = "{fmt} is an open-source formatting library providing a fast and safe alternative to C stdio and C++ iostreams."
    settings = "os", "compiler", "build_type", "arch"
    options = {"shared": [True, False]}
    default_options = {"shared": False}

    generators = "cmake"

    def source(self):
        git = tools.Git(folder="fmt")
        git.clone(self.url, self.version)
        # Garantee proper linkage
        tools.replace_in_file("fmt/CMakeLists.txt", "project(FMT CXX)",
                              '''project(FMT CXX)
include("${CMAKE_BINARY_DIR}/conanbuildinfo.cmake")
conan_basic_setup()''')

    def _configure_cmake(self):
        cmake = CMake(self)
        cmake.configure(source_folder="fmt", defs={
            "BUILD_SHARED_LIBS": "ON" if self.options.shared else "OFF",
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
