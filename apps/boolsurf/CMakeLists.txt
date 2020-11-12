add_executable(boolsurf  boolsurf.cpp myfile.cpp)

set_target_properties(boolsurf  PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(boolsurf  PRIVATE ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(boolsurf  yocto yocto_gui)

