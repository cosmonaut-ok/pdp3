cmake_minimum_required(VERSION 2.8.4)
project(pdp3)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    Beam.cpp
    Beam.h
    Boundary_conditions.cpp
    Boundary_conditions.h
    boundary_dirichlet.cpp
    boundary_dirichlet.h
    Boundary_Maxwell_conditions.cpp
    Boundary_Maxwell_conditions.h
    boundary_neumann.cpp
    boundary_neumann.h
    Bunch.cpp
    Bunch.h
    charge_density.cpp
    charge_density.h
    Constant.h
    current.cpp
    current.h
    E_field.cpp
    E_field.h
    field.cpp
    field.h
    Fourier.cpp
    Fourier.h
    Geometry.cpp
    Geometry.h
    H_field.cpp
    H_field.h
    input_output_class.cpp
    input_output_class.h

    Load_init_param.cpp
    Load_init_param.h
    Main.cpp
    Particles.cpp
    Particles.h
    particles_list.cpp
    particles_list.h
    particles_struct.cpp
    particles_struct.h
    pdp3_time.h
    PML.cpp
    PML.h
    Poisson.cpp
    Poisson.h
    Poisson_dirichlet.cpp
    poisson_dirichlet.h
    Poisson_neumann.cpp
    poisson_neumann.h
    resource.h
    resource1.h
    Test_class.cpp
    Test_class.h
    Time.cpp
    tinystr.cpp
    tinystr.h
    tinyxml.cpp
    tinyxml.h
    tinyxmlerror.cpp
    tinyxmlparser.cpp
    Triple.cpp
    Triple.h
    #wrapper.cpp
    wrapper.h)

add_executable(pdp3 ${SOURCE_FILES})
