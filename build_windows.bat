cd /d %~dp0
if %ARCH% == x86 (
    cd src
    if %PREC% == DOUBLE (
        python ..\tools\waf configure --msvc_version="msvc 15.0" --msvc_target="x86" --double
        python ..\tools\waf -v
        ren build\vclib_double.dll x86_vclib_double.dll
    ) else (
        python ..\tools\waf configure --msvc_version="msvc 15.0" --msvc_target="x86"
        python ..\tools\waf -v
        ren build\vclib_float.dll x86_vclib_float.dll
    )
) else (
    cd src
    if %PREC% == DOUBLE (
        python ..\tools\waf configure --msvc_version="msvc 15.0" --msvc_target="x86_amd64" --double
        python ..\tools\waf -v
        ren build\vclib_double.dll x64_vclib_double.dll
    ) else (
        python ..\tools\waf configure --msvc_version="msvc 15.0" --msvc_target="x86_amd64"
        python ..\tools\waf -v
        ren build\vclib_float.dll x64_vclib_float.dll
    )
)
cd ..\