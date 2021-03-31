add_executable(gui  gui.cpp)

set_target_properties(gui  PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(gui  PRIVATE ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(gui  yocto yocto_gui boolsurf)

