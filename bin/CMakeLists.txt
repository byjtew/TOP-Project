add_executable(main main.c)
target_link_libraries(main lbm m)

message(STATUS "Copying files...")
configure_file(${CONFIG_FILE} ${CMAKE_BINARY_DIR}/bin/config.txt COPYONLY)

add_executable(display display.c)
target_link_libraries(display lbm m)