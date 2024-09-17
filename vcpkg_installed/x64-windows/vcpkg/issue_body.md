Package: nodesoup:x64-windows@2023-06-12

**Host Environment**

- Host: x64-windows
- Compiler: MSVC 19.41.34120.0
-    vcpkg-tool version: 2024-02-07-8a83681f921b10d86ae626fd833c253f4f8c355b
    vcpkg-scripts version: 345ac44ab 2024-02-09 (7 months ago)

**To Reproduce**

`vcpkg install `
**Failure logs**

```
-- Note: nodesoup only supports static library linkage. Building static library.
-- Downloading https://github.com/olvb/nodesoup/archive/3158ad082bb0cd1abee75418b12b35522dbca74f.tar.gz -> olvb-nodesoup-3158ad082bb0cd1abee75418b12b35522dbca74f.tar.gz...
[DEBUG] To include the environment variables in debug output, pass --debug-env
[DEBUG] Disabling metrics because vcpkg.disable-metrics exists
[DEBUG] Trying to load bundleconfig from C:\vcpkg\vcpkg-bundle.json
[DEBUG] Failed to open: C:\vcpkg\vcpkg-bundle.json
[DEBUG] Bundle config: readonly=false, usegitregistry=false, embeddedsha=nullopt, deployment=Git, vsversion=nullopt
[DEBUG] Feature flag 'binarycaching' unset
[DEBUG] Feature flag 'compilertracking' unset
[DEBUG] Feature flag 'registries' unset
[DEBUG] Feature flag 'versions' unset
[DEBUG] Feature flag 'dependencygraph' unset
Downloading https://github.com/olvb/nodesoup/archive/3158ad082bb0cd1abee75418b12b35522dbca74f.tar.gz
warning: Download failed -- retrying after 1000ms
warning: Download failed -- retrying after 2000ms
warning: Download failed -- retrying after 4000ms
error: Failed to download from mirror set
error: https://github.com/olvb/nodesoup/archive/3158ad082bb0cd1abee75418b12b35522dbca74f.tar.gz: WinHttpSendRequest failed with exit code 12002
error: https://github.com/olvb/nodesoup/archive/3158ad082bb0cd1abee75418b12b35522dbca74f.tar.gz: WinHttpSendRequest failed with exit code 12002
error: https://github.com/olvb/nodesoup/archive/3158ad082bb0cd1abee75418b12b35522dbca74f.tar.gz: WinHttpSendRequest failed with exit code 12002
error: https://github.com/olvb/nodesoup/archive/3158ad082bb0cd1abee75418b12b35522dbca74f.tar.gz: WinHttpSendRequest failed with exit code 12002
[DEBUG] D:\a\_work\1\s\src\vcpkg\base\downloads.cpp(1030): 
[DEBUG] Time in subprocesses: 0us
[DEBUG] Time in parsing JSON: 4us
[DEBUG] Time in JSON reader: 0us
[DEBUG] Time in filesystem: 111us
[DEBUG] Time in loading ports: 0us
[DEBUG] Exiting after 1.5 min (91244753us)

CMake Error at scripts/cmake/vcpkg_download_distfile.cmake:32 (message):
      
      Failed to download file with error: 1
      If you are using a proxy, please check your proxy setting. Possible causes are:
      
      1. You are actually using an HTTP proxy, but setting HTTPS_PROXY variable
         to `https://address:port`. This is not correct, because `https://` prefix
         claims the proxy is an HTTPS proxy, while your proxy (v2ray, shadowsocksr
         , etc..) is an HTTP proxy. Try setting `http://address:port` to both
         HTTP_PROXY and HTTPS_PROXY instead.
      
      2. If you are using Windows, vcpkg will automatically use your Windows IE Proxy Settings
         set by your proxy software. See https://github.com/microsoft/vcpkg-tool/pull/77
         The value set by your proxy might be wrong, or have same `https://` prefix issue.
      
      3. Your proxy's remote server is out of service.
      
      If you've tried directly download the link, and believe this is not a temporary
      download server failure, please submit an issue at https://github.com/Microsoft/vcpkg/issues
      to report this upstream download server failure.
      

Call Stack (most recent call first):
  scripts/cmake/vcpkg_download_distfile.cmake:270 (z_vcpkg_download_distfile_show_proxy_and_fail)
  scripts/cmake/vcpkg_from_github.cmake:106 (vcpkg_download_distfile)
  C:/Users/kcscs/AppData/Local/vcpkg/registries/git-trees/43d30f8e5e0cd4b45a590027d0c9ea5036f0ba40/portfile.cmake:3 (vcpkg_from_github)
  scripts/ports.cmake:170 (include)



```
**Additional context**

<details><summary>vcpkg.json</summary>

```
{
  "dependencies": [
    "eigen3",
    "matplotplusplus"
  ]
}

```
</details>
