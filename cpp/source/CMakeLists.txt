# add executable main
add_executable(main main.cpp)
# link hdf5 library
target_link_libraries(main ${HDF5_LIBRARIES})
# link mcnptools library
target_link_libraries(main mcnptools)
