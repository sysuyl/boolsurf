add_subdirectory(ext/CDT/CDT)

add_library(boolsurf
  boolsurf_utils.h
  boolsurf.h boolsurf.cpp
  boolsurf_io.h boolsurf_io.cpp
)

set_target_properties(boolsurf PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(boolsurf PUBLIC ${CMAKE_SOURCE_DIR}/libs)

# warning flags
if(APPLE)
  target_compile_options(boolsurf PUBLIC -Wall -Wconversion -Wno-sign-conversion -Wno-implicit-float-conversion)
endif(APPLE)
if(MSVC)
  # target_compile_options(boolsurf PUBLIC /D_CRT_SECURE_NO_WARNINGS /wd4018 /wd4244 /wd4305 /wd4800 /wd4267)
  target_compile_options(boolsurf PUBLIC /D_CRT_SECURE_NO_WARNINGS)
  target_compile_options(boolsurf PUBLIC /EHsc)
endif(MSVC)