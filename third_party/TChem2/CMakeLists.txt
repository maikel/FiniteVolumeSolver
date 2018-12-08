cmake_minimum_required(VERSION 3.10)

project(TChem C)

add_library(TChem
    src/TC_chg.c
    src/TC_edit.c
    src/TC_for.c
    src/TC_getthc9.c
    src/TC_info.c
    src/TC_init.c
    src/TC_jac.c
    src/TC_kmodint.c
    src/TC_mlms.c
    src/TC_rr.c
    src/TC_spec.c
    src/TC_src.c
    src/TC_thermo.c
    src/TC_utils.c)
target_include_directories(TChem 
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  PRIVATE
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>
)
install(
  TARGETS TChem
  EXPORT ReactiveEulerSolver
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
install(
  FILES ${CMAKE_CURRENT_SOURCE_DIR}/include/*.h
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include
)