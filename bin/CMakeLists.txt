add_executable(main main.c)
target_link_libraries(main lbm m)

add_executable(display display.c)
target_link_libraries(display lbm m)