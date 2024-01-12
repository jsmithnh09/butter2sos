#-----------------------------------------------
# BUTTER2SOS/BUTTERBAND single band program
#-----------------------------------------------
add_executable(butter2sos butter2sos.c)
add_executable(butterband butterband.c)

target_link_libraries(butter2sos PRIVATE design)
target_link_libraries(butterband PRIVATE design)

target_include_directories(butter2sos PUBLIC ../lib)
target_include_directories(butterband PUBLIC ../lib)

install(TARGETS butter2sos butterband
    RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/cmake-bin)