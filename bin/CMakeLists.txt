add_executable(main main.c)
target_link_libraries(main lbm m)
if (${CMAKE_BUILD_TYPE} MATCHES "Release")
    target_compile_options(main PRIVATE -DRELEASE_MODE)
endif ()

message(STATUS "Copying files...")
configure_file(${CONFIG_FILE} ${CMAKE_BINARY_DIR}/bin/config.txt COPYONLY)

add_executable(display display.c)
target_link_libraries(display lbm m)