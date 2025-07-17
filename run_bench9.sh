#!/bin/bash

if [ $# -ne 1 ]; then
        echo "please specify 1 command line argument"
		exit 1
fi

varname="bench"
varname="${varname}$1"
# echo "$varname" 

cmake --build build/
./build/bin/"$varname"