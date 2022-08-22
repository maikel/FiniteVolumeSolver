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

# FindSundials
# ------------
#
# Find the Sundials packaged libraries.
#
# This will define the following variables:
#
#   Sundials_FOUND    - True if the system has the Foo library
#   Sundials_VERSION  - The version of the Foo library which was found
#
# and the following imported targets:
#
#   Sundials::cvode  - The CVode library

find_package(PkgConfig)
pkg_check_modules(PC_Sundials QUIET Sundials)

set(Sundials_SEARCH_DIR ${Sundials_ROOT})

find_path(Sundials_INCLUDE_DIR
  NAMES cvode/cvode.h
  HINTS ${Sundials_SEARCH_DIR}
  PATHS ${PC_Sundials_INCLUDE_DIRS}
  PATH_SUFFIXES include)

find_library(Sundials_cvode_LIBRARY 
  NAMES sundials_cvode
  HINTS ${Sundials_SEARCH_DIR}
  PATHS ${PC_Sundials_LIBRARY_DIRS} 
  ENV LIBRARY_PATH
  PATH_SUFFIXES lib)

set(Sundials_VERSION ${PC_Sundials_VERSION})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Sundials
  FOUND_VAR Sundials_FOUND
  REQUIRED_VARS 
    Sundials_LIBRARY
    Sundials_INCLUDE_DIR
  VERSION_VAR Sundials_VERSION
)

if (Sundials_FOUND AND NOT TARGET Sundials::cvode)
  add_library(Sundials::cvode UNKNOWN IMPORTED)
  set_target_properties(Sundials::cvode PROPERTIES
      INTERFACE_COMPILE_DEFINITIONS  "${PC_Sundials_CFLAGS_OTHER}"
      INTERFACE_INCLUDE_DIRECTORIES  "${Sundials_INCLUDE_DIR}"
      IMPORTED_LOCATION              "${Sundials_cvode_LIBRARY}"
  )
endif()

set(Sundials_FOUND Sundials_cvode_FOUND)