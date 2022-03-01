::No spaces!!
@echo off
if exist data\ (
    cd data
    del *.dat
    del *.csv
    cd ..
)

setlocal
set dt=1e-15
set steps=1e6
set J=2
set lambda=1e-3
set B=0
set anisotropyAxis=0
set anisotropyPlane=0
set T=0
REM 0: Z-direction, 1: angle with z-axis, 2: small angle with z-axis, 3: rotor, 4: Z-direction alternating
set init=3
set angle=30
set mode=2
REM single, line, square, triangle, kagome
set structure=line
REM always have even number of cells
set nCellsX=2 
set periodicBoundary=true

for %%i in (0,10,50,70) do ^
model.exe -dt %dt% -steps %steps% -J %J% -lambda %lambda% -B %B% ^
-anisotropyAxis %anisotropyAxis% -anisotropyPlane %anisotropyPlane% ^
-T %%i -init %init% -angle %%i -mode %mode% -nCellsX %nCellsX% -periodicBoundary %periodicBoundary% -structure %structure%
endlocal
