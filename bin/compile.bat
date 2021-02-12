@echo off
cd ..\cpp\source
g++ -static -Wall main.cpp vectors.cpp matrix.cpp operators.cpp h_transform.cpp algo.cpp -o ../../bin/main
pause