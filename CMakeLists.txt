cmake_minimum_required(VERSION 3.14)

project (ArlequinProject LANGUAGES C CXX)


set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# this option is useful for static analyzers
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


# Finding the pz package
find_package(NeoPZ REQUIRED HINTS ../NeoPZ_install ../neopz_install)

# creating a common library for all the targets in this repo
add_library(arlequin SHARED "")
target_link_libraries(arlequin PUBLIC NeoPZ::pz)
target_include_directories(arlequin PUBLIC ${PZ_INCLUDE_DIRS})

add_subdirectory(arlequin)
get_target_property(lib_files arlequin SOURCES)

source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${lib_files})  


# creating the projects
add_subdirectory(targets)

#To enable unit test using Catch2 UnitTest framework
# option(BUILD_TESTS "Whether the Catch2 UnitTests will be linked in" OFF)
# if(BUILD_TESTS)
#   include(cmake/add_unit_test.cmake)
#   include(cmake/EnableCatch2.cmake)
#   enable_catch2()
#   add_subdirectory(test)
#   ENABLE_TESTING()
#   include(CTest)
# endif()