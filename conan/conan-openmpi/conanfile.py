import os
from conans import ConanFile, tools, AutoToolsBuildEnvironment


class OpenMPIConan(ConanFile):
    name = "OpenMPI"
    version = "3.1.3"
    url = "<Package recipe repository url here, for issues about the package>"
    description = "A High Performance Message Passing Library"
    license = "https://www.open-mpi.org/community/license.php"
    settings = "os", "arch", "compiler", "build_type"
    options = {"shared": [True, False],
               "fPIC": [True, False],
               "oshmem": [True, False],
               "fortran": ['yes', 'mpifh', 'usempi', 'usempi80', 'no']}
    default_options = "shared=False", "fPIC=False", "oshmem=False", "fortran=no"

    def config(self):
        del self.settings.compiler.libcxx

    def requirements(self):
        self.requires.add("zlib/[>=1.2.11]@conan/stable")

    def source(self):
        version_tokens = self.version.split('.')
        version_short = 'v%s.%s' % (version_tokens[0], version_tokens[1])
        source_url = "https://download.open-mpi.org/release/open-mpi"
        tools.get("{0}/{1}/{2}-{3}.tar.bz2".format(source_url, version_short, self.name.lower(), self.version))
        extracted_dir = self.name.lower() + "-" + self.version
        os.rename(extracted_dir, "openmpi")

    def build(self):
        with tools.chdir("openmpi"):
            env_build = AutoToolsBuildEnvironment(self)
            env_build.fpic = self.options.fPIC

            args = ['prefix=%s' % self.package_folder]
            if self.settings.build_type == 'Debug':
                args.append('--enable-debug')
            if self.options.shared:
                args.extend(['--enable-shared', '--disable-static'])
            else:
                args.extend(['--enable-static', '--disable-shared'])
            args.append('--with-pic' if self.options.fPIC else '--without-pic')
            args.append('--enable-mpi-fortran=%s' % str(self.options.fortran))
            args.append('--with-zlib=%s' % self.deps_cpp_info['zlib'].rootpath)
            args.append('--with-zlib-libdir=%s' % self.deps_cpp_info['zlib'].lib_paths[0])
            if not self.options.oshmem:
                args.append('--disable-oshmem')

            env_build.configure(args=args)
            env_build.make()
            env_build.make(args=['install'])

    def package(self):
        self.copy(pattern="LICENSE", src='sources')

    def package_info(self):
        self.cpp_info.libs = ['mpi', 'open-rte', 'open-pal']
        if self.settings.os == "Linux":
            self.cpp_info.libs.extend(['dl', 'pthread', 'rt', 'util'])
        self.env_info.MPI_HOME = self.package_folder
        self.env_info.OPAL_PREFIX = self.package_folder
        mpi_bin = os.path.join(self.package_folder, 'bin')
        self.env_info.MPI_BIN = mpi_bin
        self.env_info.PATH.append(mpi_bin)
