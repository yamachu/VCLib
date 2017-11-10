cd /d %~dp0
if %ARCH% == x86 (
    cd src
    python ..\tools\waf configure
    python ..\tools\waf -v
    ren build\vclib.dll x86_vclib.dll
) else (
    cd src
    python ..\tools\waf configure
    python ..\tools\waf -v
    ren build\vclib.dll x64_vclib.dll
)
cd ..\