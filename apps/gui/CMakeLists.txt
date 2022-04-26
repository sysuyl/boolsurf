add_executable(gui  gui.cpp)

set_target_properties(gui  PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(gui  PRIVATE ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(gui  yocto yocto_gui boolsurf)

if(WIN32)
  target_compile_definitions(gui PUBLIC "NOMINMAX")
  target_compile_definitions(gui PUBLIC "WIN32_LEAN_AND_MEAN")
endif()
