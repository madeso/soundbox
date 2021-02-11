set(src
    main.cc
    
    soundbox.cc soundbox.h
    types.h
)

add_executable(workbench ${src})
target_link_libraries(
    workbench
    PUBLIC fmt::fmt
    PRIVATE project_options # project_warnings
)
source_group("" FILES ${src})
