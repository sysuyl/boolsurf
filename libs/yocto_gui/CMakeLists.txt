if(YOCTO_OPENGL)
  find_package(OpenGL REQUIRED)

if (WIN32)
  add_library(yocto_gui 
    yocto_opengl.h 
    yocto_opengl.cpp 
    yocto_shade.h 
    yocto_shade.cpp 
    yocto_window.h 
    yocto_window.cpp 
    yocto_imgui.h 
    yocto_imgui.cpp 
    ext/imgui/imgui.cpp 
    ext/imgui/imgui_draw.cpp 
    ext/imgui/imgui_widgets.cpp
    ext/imgui/imgui_impl_glfw.cpp 
    ext/imgui/imgui_impl_opengl3.cpp 
    ext/glad/glad.c)
else()
 add_library(yocto_gui 
    yocto_opengl.h 
    yocto_opengl.cpp 
    yocto_font.h 
    yocto_font.cpp 
    yocto_shade.h 
    yocto_shade.cpp 
    yocto_window.h 
    yocto_window.cpp 
    yocto_imgui.h 
    yocto_imgui.cpp 
    ext/imgui/imgui.cpp 
    ext/imgui/imgui_draw.cpp 
    ext/imgui/imgui_widgets.cpp
    ext/imgui/imgui_impl_glfw.cpp 
    ext/imgui/imgui_impl_opengl3.cpp 
    ext/glad/glad.c)
endif()
  
  target_include_directories(yocto_gui PRIVATE .)
  target_include_directories(yocto_gui PRIVATE ./ext)
  target_include_directories(yocto_gui PRIVATE ${FREETYPE_INCLUDE_DIRS})
  target_include_directories(yocto_gui PUBLIC ${CMAKE_SOURCE_DIR}/libs)
  target_compile_definitions(yocto_gui PRIVATE -DIMGUI_IMPL_OPENGL_LOADER_GLAD)
  target_compile_definitions(yocto_gui PUBLIC -DYOCTO_OPENGL)
  set_target_properties(yocto_gui PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

  target_link_libraries(yocto_gui PUBLIC yocto)

  target_include_directories(yocto_gui PRIVATE ${CMAKE_SOURCE_DIR}/exts/glfw/include)
  target_link_libraries(yocto_gui PUBLIC glfw ${OPENGL_gl_LIBRARY})

  if(UNIX AND NOT APPLE)
    target_link_libraries(yocto_gui PUBLIC dl)
  endif(UNIX AND NOT APPLE)
  
  if(WIN32)
    #target_include_directories(yocto_gui PRIVATE ${CMAKE_SOURCE_DIR}/exts/freetype2/include)
    #link_directories(yocto_gui PRIVATE ${CMAKE_SOURCE_DIR}/exts/freetype2/release dll/win64)
  else()
    find_package(Freetype REQUIRED)
    target_include_directories(yocto_gui PRIVATE ${FREETYPE_INCLUDE_DIRS})
    target_link_libraries(yocto_gui PRIVATE ${FREETYPE_LIBRARIES})
  endif()
  
  if(UNIX AND NOT APPLE)
    target_link_libraries(yocto_gui PRIVATE dl)
  endif(UNIX AND NOT APPLE)
endif(YOCTO_OPENGL)

# warning flags
if(APPLE)
  target_compile_options(yocto_gui PUBLIC -Wall -Wconversion -Wno-sign-conversion -Wno-implicit-float-conversion)
endif(APPLE)
if(MSVC)
target_compile_options(yocto_gui PUBLIC /D_CRT_SECURE_NO_WARNINGS /wd4018 /wd4244 /wd4305 /wd4800 /wd4267)
target_compile_options(yocto_gui PUBLIC /EHsc)
endif(MSVC)
