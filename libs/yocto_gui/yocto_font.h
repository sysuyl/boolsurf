#pragma once

#ifdef _WIN32
#undef near
#undef far
#endif

#include <yocto_gui/yocto_opengl.h>

#include <unordered_map>
namespace yocto {

/// Holds all state information relevant to a character as loaded using FreeType
struct opengl_char {
  vec2i        size    = {0, 0};  // size of glyph
  vec2i        bearing = {0, 0};  // offset from baseline to left/top of glyph
  uint         advance = 0;       // horizontal offset to advance to next glyph
  ogl_texture* texture = new ogl_texture{};
};

struct opengl_font {
  uint                size       = 0;
  ogl_shape*          quad       = new ogl_shape{};
  ogl_program*        program    = new ogl_program{};
  vector<opengl_char> characters = {};
};

void init_font(opengl_font* font, const string& filename, uint size);

void draw_text(const opengl_font* font, const string& text, float x, float y,
    float scale, const vec3f& color, float alpha = 1);

inline void draw_text(const opengl_font* font, const vector<string>& lines,
    float x, float y, float scale, const vec3f& color, float alpha = 1) {
  float down_step = 2 * font->characters.at('I').size.y * scale / font->size;
  for (auto& line : lines) {
    draw_text(font, line, x, y, scale, color, alpha);
    y -= down_step;
  }
}

}  // namespace yocto
