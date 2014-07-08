copy ..\..\..\wrappers\MATLAB\example.m .
copy ..\..\..\wrappers\MATLAB\*rops.c .
matlab -wait -nojvm -nodesktop -r "MATLABBuilder; quit"
cl /c /I../../../CoolProp /EHsc CoolProp_wrap.cxx
matlab -nojvm -nodesktop -r "Example; quit" -logfile Output.txt
cl /c /I../../../CoolProp /EHsc CoolProp_wrap.cxx
erase *Props*.c