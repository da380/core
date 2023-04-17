cmake_minimum_required(VERSION 3.5)


project(core)
enable_language(Fortran)

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(dialect " -O3")
    set(bounds "-fbounds-check")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(dialect "-stand f08 -free -implicitnone")
    set(bounds "-check bounds")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
    set(dialect "-Mfreeform -Mdclchk -Mstandard -Mallocatable=03")
    set(bounds "-C")
  endif()



set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")


find_package(LAPACK REQUIRED)
find_package(FFTW3 COMPONENTS double single threads REQUIRED)
include_directories(${FFTW3_INCLUDE_DIR})
find_package(NETCDFF REQUIRED)
include_directories(${NETCDFF_INCLUDE_DIR})


SET(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_SOURCE_DIR}/../mod)


list(APPEND SRC "module_constants.f90")
list(APPEND SRC "module_physical_constants.f90")
list(APPEND SRC "module_error.f90")
list(APPEND SRC "module_util.f90")
list(APPEND SRC "module_fftw3.f90")
list(APPEND SRC "module_interp.f90")
list(APPEND SRC "module_special_functions.f90")
list(APPEND SRC "module_quadrature.f90")
list(APPEND SRC "module_spherical_harmonics.f90")
list(APPEND SRC "module_spherical_model.f90")
list(APPEND SRC "module_PREM.f90")
list(APPEND SRC "module_DECK.f90")
list(APPEND SRC "module_spherical_mesh.f90")
#list(APPEND SRC "module_random_fields.f90")
#list(APPEND SRC "module_LAPACK.f90")
#list(APPEND SRC "module_SEM_1D.f90")


list(APPEND LIBS ${LAPACK_LIBRARIES})
list(APPEND LIBS ${FFTW3_LIBRARIES})
list(APPEND LIBS ${NETCDFF_LIBRARIES})

set(TEST_DIR ${CMAKE_SOURCE_DIR}/tests/)

add_executable(test ${TEST_DIR}test.f90 ${SRC})
target_link_libraries(test PRIVATE ${LIBS})


add_executable(test_spherical_harmonics ${TEST_DIR}test_spherical_harmonics.f90 ${SRC})
target_link_libraries(test_spherical_harmonics PRIVATE ${LIBS})







