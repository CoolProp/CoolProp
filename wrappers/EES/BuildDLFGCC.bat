REM ******* compile all the sources from CoolProp***************
g++ -c -Wall -I../../CoolProp ../../CoolProp/*.cpp
g++ -c -Wall -I../../CoolProp main.cpp
g++ -shared -o COOLPROP_EES.dlf *.o -Wl
erase *.o