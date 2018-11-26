# Copyright (c) 2018 Maikel Nadolski
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# FindTBB
# ------------
#
# Find the TBB packaged libraries.
#
# This will define the following variables:
#
#   TBB_FOUND    - True if the system has the Foo library
#   TBB_VERSION  - The version of the Foo library which was found
#
# and the following imported targets:
#
#   TBB::tbb  - The main tbb library

find_package(PkgConfig)
pkg_check_modules(PC_TBB QUIET tbb)

set(TBB_SEARCH_DIR ${TBB_ROOT})

find_path(TBB_INCLUDE_DIR
  NAMES tbb/tbb.h
  HINTS ${TBB_SEARCH_DIR}
  PATHS ${PC_TBB_INCLUDE_DIRS}
  PATH_SUFFIXES include)

find_library(TBB_LIBRARY 
  NAMES tbb
  HINTS ${TBB_SEARCH_DIR}
  PATHS ${PC_TBB_LIBRARY_DIRS} 
  ENV LIBRARY_PATH
  PATH_SUFFIXES lib)

set(TBB_VERSION ${PC_TBB_VERSION})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(TBB
  FOUND_VAR TBB_FOUND
  REQUIRED_VARS 
    TBB_LIBRARY
    TBB_INCLUDE_DIR
  VERSION_VAR TBB_VERSION
)

if (TBB_FOUND AND NOT TARGET TBB::tbb)
  add_library(TBB::tbb UNKNOWN IMPORTED)
  set_target_properties(TBB::tbb PROPERTIES
      INTERFACE_COMPILE_DEFINITIONS  "${PC_TBB_CFLAGS_OTHER}"
      INTERFACE_INCLUDE_DIRECTORIES  "${TBB_INCLUDE_DIR}"
      IMPORTED_LOCATION              "${TBB_LIBRARY}"
  )
endif()