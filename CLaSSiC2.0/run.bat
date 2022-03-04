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
set J=-2
set lambda=0
set B=0
set anisotropyAxis=0
set anisotropyPlane=10
set T=0
REM 0: Z-direction, 1: angle with z-axis, 2: small angle with z-axis, 3: rotor, 4: Z-direction alternating
set init=7
set angle=0
set mode=0
REM single, line, square, triangle, kagome
set structure=kagome
REM always have even number of cells
set nCellsX=12 
set periodicBoundary=true
set stabilize=false

for %%i in (0) do ^
model.exe -dt %dt% -steps %steps% -J %J% -lambda %lambda% -B %B% ^
-anisotropyAxis %anisotropyAxis% -anisotropyPlane %anisotropyPlane% ^
-T %T% -init %init% -angle %angle% -mode %mode% -nCellsX %nCellsX% ^
-periodicBoundary %periodicBoundary% -structure %structure% -stabilize %stabilize%
endlocal
