add_executable(main main.c)
target_link_libraries(main PUBLIC lbm m OpenMP::OpenMP_C Threads::Threads)
if (${CMAKE_BUILD_TYPE} MATCHES "Release" OR ${CMAKE_BUILD_TYPE} MATCHES "RelWithDebInfo")
    target_compile_options(main PRIVATE -DRELEASE_MODE)
endif ()

add_executable(display display.c)
target_link_libraries(display lbm m)