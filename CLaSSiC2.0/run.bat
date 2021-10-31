::No spaces!!
@echo off
cd data
del *.dat
cd ..

setlocal
set dt=1e-14
set steps=1e6
set J=2
set lambda=0
set B=0
set anisotropy=0
set T=0
set init=0
set angle=45
set structure=single


for %%t in (0 10 20) do ^
model.exe -dt %dt% -steps %steps% -J %J% -lambda %lambda% -B %B% ^
-anisotropy %anisotropy% -T %T% -init %init% -angle %%t -structure %structure% echo %%t
endlocal
