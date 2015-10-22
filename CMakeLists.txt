cmake_minimum_required (VERSION 2.8)
project (ssa)
include_directories ("${PROJECT_SOURCE_DIR}")
include_directories ("${PROJECT_SOURCE_DIR}/yaml-0.1.6/include")
include_directories ("${PROJECT_SOURCE_DIR}/yaml-0.1.6/src")

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DHAVE_CONFIG_H")
	
add_library(yaml STATIC yaml-0.1.6/src/api.c
                        yaml-0.1.6/src/dumper.c
                        yaml-0.1.6/src/emitter.c
                        yaml-0.1.6/src/loader.c
                        yaml-0.1.6/src/parser.c
                        yaml-0.1.6/src/reader.c
                        yaml-0.1.6/src/scanner.c
                        yaml-0.1.6/src/writer.c)

add_executable(demo_ssa demo/demo_sim.cc rdmini/rdmodel.cc rdmini/yamlview.cc )
target_link_libraries(demo_ssa yaml)
