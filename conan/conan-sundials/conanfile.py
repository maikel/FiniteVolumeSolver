from conans import ConanFile, CMake, tools


class SundialsConan(ConanFile):
    name = "SUNDIALS"
    version = "4.0.1"
    license = "https://github.com/LLNL/sundials/blob/v4.0.1/LICENSE"
    url = "https://git.imp.fu-berlin.de/denzler/conanpackages"
    description = "SUNDIALS is a SUite of Nonlinear and DIfferential/ALgebraic equation Solvers"
    settings = "os", "compiler", "build_type", "arch"
    options = {"shared": [True, False]}
    default_options = "shared=False"
    generators = "cmake"

    def source(self):
        self.run("git clone https://github.com/LLNL/sundials.git --branch v4.0.1 --single-branch --depth=1")   
        # Garantee proper linkage
        tools.replace_in_file("sundials/CMakeLists.txt", "PROJECT(sundials C)",
                              '''PROJECT(sundials C)
include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup()''')

    def _configure_cmake(self):
        cmake = CMake(self)
        cmake.configure(source_folder="sundials", defs={
            "BUILD_SHARED_LIBS": "ON" if self.options.shared else "OFF",
            "EXAMPLES_ENABLE_C": "OFF",
            "EXAMPLES_ENABLE_CXX": "OFF",
            "EXAMPLES_INSTALL": "OFF",
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

