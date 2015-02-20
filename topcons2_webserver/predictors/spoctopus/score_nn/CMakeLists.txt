add_executable (pred pred.c nets.c)
target_link_libraries (pred m)
install (TARGETS pred DESTINATION spoctopus/score_nn)
