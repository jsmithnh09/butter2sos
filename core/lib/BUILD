cc_shared_library(
    name = "butterlib",
    shared_lib_name = "butterlib.so",
    deps = [":design"],
    visibility = ["//main:__pkg__"],
)

cc_library(
    name = "design",
    srcs = ["butter2sos_design.c"],
    hdrs = ["butter2sos_design.h"],
)