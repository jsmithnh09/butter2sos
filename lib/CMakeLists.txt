#-------------------------------------------
# BUTTERLIB static library for programs.
#-------------------------------------------
add_library(design STATIC butter2sos_design.c)
target_compile_features(design PUBLIC cxx_std_11)
if (WIN32 AND MSVC)
    target_compile_options(design PUBLIC /W4)
else()
    target_compile_options(design PUBLIC -Wall -O2)
endif()

#--------------------------------------------------
# Generating the DLL for the program.
#--------------------------------------------------
add_library(butterlib SHARED butter2sos_design.c design.def)
target_compile_features(butterlib PUBLIC cxx_std_11)
if (WIN32 AND MSVC)
    target_compile_options(butterlib PUBLIC /W4)
else()
    target_compile_options(butterlib PUBLIC -Wall -O2)
endif()
install(TARGETS butterlib RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/src/pybutter)