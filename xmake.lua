add_rules("mode.release")
add_requires("cnpy", "eigen", "openmp")

target("HW6")
    set_kind("binary")
    add_files("*.cpp")
    add_packages("cnpy")
    add_packages("eigen")
    add_packages("openmp")

    -- set working directory
    set_rundir("$(projectdir)")

