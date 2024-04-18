cmake_minimum_required(VERSION 3.5)

project(core_modules)
enable_language(Fortran)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
SET(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(flags " -O3")
    set(debug_flags "-bounds-check")
endif()

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${debug_flags}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${flags}")

#find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(FFTW3  REQUIRED)
include_directories(${FFTW3_INCLUDE_DIR})
find_package(NETCDFF REQUIRED)
include_directories(${NETCDFF_INCLUDE_DIR})

list(APPEND LIBS ${LAPACK_LIBRARIES})
list(APPEND LIBS ${FFTW3_LIBRARIES})
list(APPEND LIBS ${NETCDFF_LIBRARIES})

SET(SRC_DIR ${CMAKE_SOURCE_DIR}/src/)
list(APPEND SRC "${SRC_DIR}module_constants.f90")
list(APPEND SRC "${SRC_DIR}module_physical_constants.f90")
list(APPEND SRC "${SRC_DIR}module_error.f90")
list(APPEND SRC "${SRC_DIR}module_util.f90")
list(APPEND SRC "${SRC_DIR}module_linalg.f90")
list(APPEND SRC "${SRC_DIR}module_fftw3.f90")
list(APPEND SRC "${SRC_DIR}module_interp.f90")
list(APPEND SRC "${SRC_DIR}module_special_functions.f90")
list(APPEND SRC "${SRC_DIR}module_quadrature.f90")
list(APPEND SRC "${SRC_DIR}module_spherical_harmonics.f90")
list(APPEND SRC "${SRC_DIR}module_SEM_1D.f90")
list(APPEND SRC "${SRC_DIR}module_random_fields.f90")

#list(APPEND SRC "${SRC_DIR}module_spherical_model.f90")
#list(APPEND SRC "${SRC_DIR}module_spherical_model_new.f90")
#list(APPEND SRC "${SRC_DIR}module_PREM.f90")
#list(APPEND SRC "${SRC_DIR}module_PREM_new.f90")
#list(APPEND SRC "${SRC_DIR}module_DECK.f90")
#list(APPEND SRC "${SRC_DIR}module_spherical_mesh.f90")


option(BUILD_SHARED_LIBS "Build using shared libraries" ON)
add_library(Core ${SRC})
target_link_libraries(Core PRIVATE ${LIBS})


set(TEST_DIR ${CMAKE_SOURCE_DIR}/src/tests/)

add_executable(test ${TEST_DIR}test.f90)
target_link_libraries(test PRIVATE ${LIBS} Core)

add_executable(test_spherical_harmonics ${TEST_DIR}test_spherical_harmonics.f90)
target_link_libraries(test_spherical_harmonics PRIVATE ${LIBS} Core)

add_executable(test_fftw3 ${TEST_DIR}test_fftw3.f90)
target_link_libraries(test_fftw3 PRIVATE ${LIBS} Core)

add_executable(test_quadrature ${TEST_DIR}test_quadrature.f90)
target_link_libraries(test_quadrature PRIVATE ${LIBS} Core)

add_executable(test_interp ${TEST_DIR}test_interp.f90)
target_link_libraries(test_interp PRIVATE ${LIBS} Core)


add_executable(test_SEM_1D ${TEST_DIR}test_SEM_1D.f90)
target_link_libraries(test_SEM_1D PRIVATE ${LIBS} Core)

add_executable(test_random_field ${TEST_DIR}test_random_field.f90)
target_link_libraries(test_random_field PRIVATE ${LIBS} Core)







