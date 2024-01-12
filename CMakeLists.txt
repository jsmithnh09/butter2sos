cmake_minimum_required(VERSION 3.26)

project(butter2sos
        VERSION 1.0.0
        DESCRIPTION "Butterworth IIR design"
        LANGUAGES C
)

enable_testing()

# generating the library as a DLL and for static linking.
add_subdirectory(lib)
add_subdirectory(programs)

# creating version header.
add_custom_target(mkversion
        ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_LIST_DIR}/mkversion.cmake)
add_dependencies(butter2sos mkversion)
add_dependencies(butterband mkversion)