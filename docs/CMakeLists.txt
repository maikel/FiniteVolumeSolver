cmake_minimum_required(VERSION 3.13)
project(FiniteVolumeSolverDoxygens)

################################################################################
# Doxygen target

find_package(Doxygen)
if (DOXYGEN_FOUND)
    # set input and output files
    set(DOXYGEN_IN ${PROJECT_SOURCE_DIR}/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    # request to configure the file
    set(FUB_SOURCE_DIR "${PROJECT_SOURCE_DIR}/../")
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    message("-- Doyxgen: Generated a Doxyfile at ${DOXYGEN_OUT}.")

    # note the option ALL which allows to build the docs together with the application
    add_custom_target( doc ALL
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM )
else (DOXYGEN_FOUND)
  message("-- Doyxgen: Doxygen needs to be installed to generate the doxygen documentation")
endif (DOXYGEN_FOUND)