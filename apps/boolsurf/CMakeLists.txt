add_subdirectory(ext/CDT/CDT)

add_executable(gui  gui.cpp boolsurf.h boolsurf.cpp boolsurf_io.h boolsurf_io.cpp)
add_executable(test  test.cpp boolsurf.h boolsurf.cpp boolsurf_io.h boolsurf_io.cpp)

set_target_properties(gui  PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DMY_DEBUG")
target_include_directories(gui  PRIVATE ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(gui  yocto yocto_gui)


set_target_properties(test  PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DMY_DEBUG")
target_include_directories(test  PRIVATE ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(test  yocto yocto_gui)

