cc_shared_library(
    name = "butterlib",
    shared_lib_name = "butterlib.dll",
    deps = [":design"],
)

cc_library(
    name = "design",
    srcs = ["butter2sos_design.c"],
    hdrs = ["butter2sos_design.h"],
    visibility = ["//main:__pkg__"],
)