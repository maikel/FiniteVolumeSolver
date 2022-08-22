import os

from conans import ConanFile, CMake, tools, RunEnvironment


class Hdf5TestConan(ConanFile):
    settings = "os", "compiler", "build_type", "arch"
    generators = "cmake"
    exports = "sample.h5"

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def imports(self):
        self.copy("*.dll", dst="bin", src="bin")
        self.copy("*.dylib*", dst="bin", src="lib")
        self.copy('*.so*', dst='bin', src='lib')

    def test(self):
        run_env = RunEnvironment(self)
        hdf5_file = os.path.join(self.source_folder, "sample.h5")
        os.chdir("bin")
        with tools.environment_append(run_env.vars):
            self.run(".%sexample %s" % (os.sep, hdf5_file))
