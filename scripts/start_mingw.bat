@echo off
chcp 65001 > nul
echo ========================================
echo    MinGW Activation - C++ Ready
echo ========================================
echo.

set "MINGW_PATH=C:\Users\ИТ Космос\winlibs-i686-posix-dwarf-gcc-16.0.0-snapshot20251026-mingw-w64ucrt-14.0.0-r1\mingw32\bin"
set "PATH=%MINGW_PATH%;%PATH%"

echo MinGW activated!
g++ --version

echo.
echo 🔧 Testing C++ compilation...
echo #include ^<iostream^> > test.cpp
echo extern "C" int add_numbers(int a, int b) { return a + b; } >> test.cpp
g++ -shared -o test_lib.dll test.cpp

if %errorlevel% equ 0 (
    echo  C++ compilation successful!
    del test.cpp test_lib.dll 2>nul
) else (
    echo C++ compilation failed
)

echo.
echo ========================================
echo     C++ DEVELOPMENT READY!
echo ========================================
echo.
cmd /k
