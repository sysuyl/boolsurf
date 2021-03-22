add_executable(gui  gui.cpp)
add_executable(test  test.cpp)

set_target_properties(gui  PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DMY_DEBUG")
target_include_directories(gui  PRIVATE ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(gui  yocto yocto_gui boolsurf)


set_target_properties(test  PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DMY_DEBUG")
target_include_directories(test  PRIVATE ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(test  yocto yocto_gui boolsurf)

