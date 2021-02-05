#pragma once

#include <yocto_gui/yocto_opengl.h>

#include <unordered_map>
using namespace yocto;

/// Holds all state information relevant to a character as loaded using FreeType
struct opengl_char {
  vec2i size       = {0, 0};  // size of glyph
  vec2i bearing    = {0, 0};  // offset from baseline to left/top of glyph
  uint  advance    = 0;       // horizontal offset to advance to next glyph
  uint  texture_id = 0;
};

struct opengl_font {
  float                                 size       = 0;
  ogl_shape*                            quad       = new ogl_shape{};
  ogl_program*                          program    = new ogl_program{};
  std::unordered_map<char, opengl_char> characters = {};
};

void init_glfont(opengl_font& font, const string& filename, float size);

void draw_glfont(const opengl_font& font, const string& text, float x, float y,
    float scale, const vec3f& color);

void draw_glfont(const opengl_font& font, const string& text,
    const frame2f& frame, float scale, float ratio, const vec3f& color,
    float alpha);
