#-------------------------------------------
# [LIB]BUTTER static library for programs.
#-------------------------------------------
add_library(design STATIC butter2sos_design.c)
target_compile_features(design PUBLIC cxx_std_11)
if (WIN32 AND MSVC)
    target_compile_options(design PUBLIC /W4)
else()
    target_compile_options(design PUBLIC -Wall -O2)
    target_link_libraries(design m)
endif()

#--------------------------------------------------
# Generating the DLL for the program.
#--------------------------------------------------
add_library(butter SHARED butter2sos_design.c design.def)
target_compile_features(butter PUBLIC cxx_std_11)
if (WIN32 AND MSVC)
    target_compile_options(butter PUBLIC /W4)
else()
    target_compile_options(butter PUBLIC -Wall -O2)
    target_link_libraries(butter m)
endif()
install(TARGETS butter DESTINATION ${CMAKE_SOURCE_DIR}/src/pybutter)
