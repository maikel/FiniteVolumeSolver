from conans import ConanFile, CMake, tools


class ReactiveeulersolverConan(ConanFile):
    name = "ReactiveEulerSolver"
    version = "0.1"
    license = "<Put the package license here>"
    url = "https://git.imp.fu-berlin.de/sfb1029/ReactiveEulerSolver"
    description = "<Description of Reactiveeulersolver here>"
    settings = "os", "compiler", "build_type", "arch"
    options = {"shared": [True, False]}
    default_options = "shared=False"
    generators = "cmake"

    def requirements(self):
        self.requires.add("OpenMPI/3.1.3@reactive-euler/stable")
        self.requires.add("HDF5/1.8.21@reactive-euler/stable")
        self.requires.add("SAMRAI/3.11.1@reactive-euler/stable")
        self.requires.add("SUNDIALS/4.0.1@reactive-euler/testing")
        self.requires.add("boost/1.70.0@conan/stable")

    def source(self):
        self.run("git clone https://git.imp.fu-berlin.de/sfb1029/ReactiveEulerSolver.git ReactiveEulerSolver")
        self.run("cd ReactiveEulerSolver && git checkout feature/conan2")
        self.run("cd ReactiveEulerSolver && git submodule update --init --recursive")

    def _configure_cmake(self):
        cmake = CMake(self)
        cmake.configure(source_folder="ReactiveEulerSolver", defs={
            "CMAKE_BUILD_TYPE": "Debug",
            "FUB_WITH_OPENMP": "On",
            "FUB_WITH_SAMRAI": "On",
            "FUB_MEX_COMPONENTS": "Off",
            "FUB_EXAMPLES": "Off"
        })

        return cmake

    def build(self):
        cmake = self._configure_cmake()
        cmake.build()

    def package(self):
        cmake = self._configure_cmake()
        cmake.install()

    def package_info(self):
        tools.collect_libs(self)
