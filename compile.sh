cd src/ilp
cmake -H. -Bbuild -DCMAKE_BUILD_TYPE=Release
cmake --build build
cd ../heuristic
g++ heuristic.cpp -o heuristic
cd ../..
