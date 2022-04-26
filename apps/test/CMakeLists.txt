add_executable(test  test.cpp)

set_target_properties(test  PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(test  PRIVATE ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(test  yocto yocto_gui boolsurf)

if(WIN32)
  target_compile_definitions(test PUBLIC "NOMINMAX")
  target_compile_definitions(test PUBLIC "WIN32_LEAN_AND_MEAN")
endif()
