call c:\gccdist\supp\setup-gcc.bat
make
nmppc PPC603gnu/CoolProp.out > exports.txt
for %%f in (PPC603gnu/*.o) do nmppc PPC603gnu/%%f > %%~nf.txt