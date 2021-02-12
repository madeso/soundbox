set(src
    soundbox/soundbox.cc soundbox/soundbox.h
    soundbox/types.h
)

add_library(soundbox STATIC ${src})
target_link_libraries(
    soundbox
    PUBLIC fmt::fmt
    PRIVATE project_options project_warnings
)
target_include_directories(soundbox
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)
source_group("" FILES ${src})
