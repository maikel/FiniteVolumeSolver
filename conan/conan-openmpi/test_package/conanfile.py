import os

from conans import ConanFile, CMake, tools, RunEnvironment


class OpenMPITestConan(ConanFile):
    settings = "os", "compiler", "build_type", "arch"
    generators = "cmake"

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def test(self):
        mpiexec = os.path.join(os.environ['MPI_BIN'], 'mpiexec')
        with open('hostfile', 'w') as f:
            f.write('localhost slots=2\n')
        command = '%s --oversubscribe -np 2 --hostfile hostfile %s' % (
            mpiexec, os.path.join("bin", "example"))
        with tools.environment_append(RunEnvironment(self).vars):
            if self.settings.os == "Windows":
                self.run(command)
            elif self.settings.os == "Macos":
                self.run("DYLD_LIBRARY_PATH=%s %s" %
                         (os.environ.get('DYLD_LIBRARY_PATH', ''), command))
            else:
                self.run("LD_LIBRARY_PATH=%s %s" %
                         (os.environ.get('LD_LIBRARY_PATH', ''), command))
