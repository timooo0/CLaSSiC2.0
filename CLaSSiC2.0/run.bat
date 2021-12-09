::No spaces!!
@echo off
cd data
del *.dat
cd ..

setlocal
set dt=1e-15
set steps=1e6
set J=2
set lambda=1e-3
set B=0
set anisotropyAxis=0
set anisotropyPlane=0
set T=1
@REM 0: Z-direction, 1: angle with z-axis, 2: small angle with z-axis, 3: rotor, 4: Z-direction alternating
set init=2
set angle=45
set mode=-3
set structure=square
set nCellsX=30

for %%i in (10) do ^
model.exe -dt %dt% -steps %steps% -J %J% -lambda %lambda% -B %%i ^
-anisotropyAxis %anisotropyAxis% -anisotropyPlane %anisotropyPlane% ^
-T %T% -init %init% -angle %angle% -mode %mode% -structure %structure% -nCellsX %nCellsX%
endlocal
