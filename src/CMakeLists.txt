FILE(GLOB sources ${SRC_DIR}/*.c)
FILE(GLOB headers ${INC_DIR}/*.h)
add_library(lbm STATIC ${sources} ${headers})
set_target_properties(lbm PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
target_include_directories(lbm PUBLIC "${INC_DIR}")