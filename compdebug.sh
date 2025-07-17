#!/bin/bash
cmake --build build/ 
time (./build/bin/bench1 && ./build/bin/bench2 && ./build/bin/bench3 && ./build/bin/bench4 && ./build/bin/bench5 && ./build/bin/bench6)

time ./build/bin/bench7
cd work
python3 Bench1Plot.py &
python3 Bench2Plot.py &
python3 Bench3Plot.py &
python3 Bench4Plot.py &
python3 Bench5Plot.py &
python3 Plot6Error.py 1 1 &
python3 Bench7Plot.py &
cd ..