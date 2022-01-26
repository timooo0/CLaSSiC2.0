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
set T=1
REM 0: Z-direction, 1: angle with z-axis, 2: small angle with z-axis, 3: rotor, 4: Z-direction alternating
set init=2
set angle=45
set mode=2
REM single, line, square, triangle, kagome
set structure=square
REM always have even number of cells
set nCellsX=2 

for %%i in (10) do ^
model.exe -dt %dt% -steps %steps% -J %J% -lambda %lambda% -B %B% ^
-anisotropyAxis %anisotropyAxis% -anisotropyPlane %anisotropyPlane% ^
-T %T% -init %init% -angle %angle% -mode %mode% -nCellsX %nCellsX% -structure %structure%
endlocal
