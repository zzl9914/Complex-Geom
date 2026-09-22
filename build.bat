@echo off
setlocal EnableExtensions
cd /d "%~dp0"

set "VCVARS="
if exist "D:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvars64.bat" (
  set "VCVARS=D:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvars64.bat"
)
if not defined VCVARS if exist "C:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvars64.bat" (
  set "VCVARS=C:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvars64.bat"
)
if not defined VCVARS if exist "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" (
  set "VCVARS=C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"
)
if not defined VCVARS (
  echo Could not find vcvars64.bat
  exit /b 1
)

call "%VCVARS%" >nul
if errorlevel 1 exit /b 1

if not exist bin mkdir bin

cl /nologo /utf-8 /std:c++17 /EHsc /O2 /MD /W3 /bigobj /openmp ^
  /DUNICODE /D_UNICODE /DGLFW_INCLUDE_NONE ^
  /I. /Isrc /Ithird_party\imgui /Ithird_party\imgui\backends /Ithird_party\glfw\include ^
  src\main.cpp scomplex.cpp ^
  third_party\imgui\imgui.cpp ^
  third_party\imgui\imgui_draw.cpp ^
  third_party\imgui\imgui_tables.cpp ^
  third_party\imgui\imgui_widgets.cpp ^
  third_party\imgui\backends\imgui_impl_glfw.cpp ^
  third_party\imgui\backends\imgui_impl_opengl3.cpp ^
  /Fobin\ /Fe:bin\ComplexGeom.exe ^
  /link /LIBPATH:third_party\glfw\lib-vc2022 glfw3.lib opengl32.lib user32.lib gdi32.lib shell32.lib winmm.lib

if errorlevel 1 exit /b 1
echo Built bin\ComplexGeom.exe
endlocal
