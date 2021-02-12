set(src
    main.cc
)

add_executable(example ${src})
target_link_libraries(
    example
    PUBLIC fmt::fmt soundbox
    PRIVATE project_options project_warnings
)
source_group("" FILES ${src})
