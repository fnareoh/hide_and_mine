cd src/ilp
cmake -H. -Bbuild -DCMAKE_BUILD_TYPE=Release
cmake --build build
cd ../heuristic
g++ -O3 heuristic.cpp -o heuristic
#g++ -g -Wall -O heuristic.cpp -o heuristic
cd ../..
