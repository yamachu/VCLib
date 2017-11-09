cd /d %~dp0
if %ARCH% == x86 (
    cd src
    python ..\tools\waf configure --msvc_version="msvc 15" --msvc_target="x86"
    python ..\tools\waf -v
    ren build\vclib.dll x86_vclib.dll
) else (
    cd src
    python ..\tools\waf configure --msvc_version="msvc 15" --msvc_target="x86_amd64"
    python ..\tools\waf -v
    ren build\vclib.dll x64_vclib.dll
)
cd ..\