# Gravitational_Field
This code is written for the calculation of graviational potential in aspherical bodies using the particle relabelling transformation method. There are several examples provided and it is designed to be built using cmake.

This package can be installed easily via:
1. git clone 
2. cd gplspec
3. cmake -S . -B build
4. cmake --build build/ 

For faster compilation (on machines that are designed to handle this) adding the flag -jN to the 4th step, where N is the number of cores to compile on, can significantly speed things up. 
