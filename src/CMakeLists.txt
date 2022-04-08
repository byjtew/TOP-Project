message(STATUS "INC_DIR: ${INC_DIR}")

add_library(lbm STATIC lbm_comm.c lbm_config.c lbm_init.c lbm_phys.c lbm_struct.c)
