#######################################################
#  MKVERSION.CMAKE updates the version header file.
#######################################################

file(STRINGS  "${CMAKE_CURRENT_LIST_DIR}/VERSION"  BUTTER_SOS_VERSION)
configure_file(
  "${CMAKE_CURRENT_LIST_DIR}/lib/version.h.in"
  "${CMAKE_CURRENT_LIST_DIR}/lib/version.h"
  )