#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#endif

#include "ext/glad/glad.h"
// We only need glPixelStorei();

#include <ft2build.h>
#include FT_FREETYPE_H

#include "yocto_font.h"

// https://learnopengl.com/In-Practice/Text-Rendering

auto vertex = R"(
  #version 330 core
  layout(location = 0) in vec2 vertex;
  out vec2     TexCoords;
  uniform vec2 scale  = vec2(1, 1);
  uniform vec2 center = vec2(0.0, 0.0);

  void main() {
    TexCoords = vec2(vertex.x, 1 - vertex.y);
    vec2 p = vertex;
    p *= scale;
    p += center;
    gl_Position = vec4(p.x, p.y, 0, 1);
  }
)";

auto fragment = R"(
  #version 330 core
  in vec2 TexCoords;
  out vec4 out_color;

  uniform sampler2D char_texture;
  uniform vec3 color = vec3(1, 1, 1);
  uniform float alpha = 1;

  void main() {    
      float sampled = texture(char_texture, TexCoords).r;
      out_color = vec4(color, alpha * sampled);
  }
)";

namespace yocto {
void init_font(opengl_font* font, const string& filename, uint size) {
  auto error    = string{};
  auto errorlog = string{};
  if (!set_program(font->program, vertex, fragment, error, errorlog)) {
    printf("%s: %s\n", __FUNCTION__, error.c_str());
    printf("%s: %s\n", __FUNCTION__, errorlog.c_str());
  }

  font->size = size;

  // FreeType
  FT_Library ft;
  if (FT_Init_FreeType(&ft))
    printf("ERROR::FREETYPE: Could not init FreeType Library\n");

  // Load font as face
  FT_Face face;
  if (FT_New_Face(ft, filename.c_str(), 0, &face))
    printf("ERROR::FREETYPE: Failed to load font %s\n", filename.c_str());

  // Set size to load glyphs as
  FT_Set_Pixel_Sizes(face, 0, font->size);

  // Disable byte-alignment restriction
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

  // Load first 128 characters of ASCII set
  font->characters.resize(128);
  for (int c = 0; c < 128; c++) {
    // Load character glyph
    if (FT_Load_Char(face, c, FT_LOAD_RENDER)) {
      printf(
          "%s: FreeType failed to load glyph \"%c\"\n", __FUNCTION__, (char)c);
      continue;
    }

    // Add entry to hash map.
    auto& character = font->characters[(char)c];
    character.size  = {
        int(face->glyph->bitmap.width), int(face->glyph->bitmap.rows)};
    character.bearing = {face->glyph->bitmap_left, face->glyph->bitmap_top};
    character.advance = (uint)face->glyph->advance.x;

    // Generate texture.
    set_texture(character.texture, character.size, 1,
        face->glyph->bitmap.buffer, true, true, false);
  }

  // Destroy FreeType
  FT_Done_Face(face);
  FT_Done_FreeType(ft);

  // Set quad shape.
  static const auto positions = vector<vec3f>{{0, 0}, {1, 0}, {0, 1}, {1, 1}};
  set_vertex_buffer(font->quad, positions, 0);
  font->quad->elements = ogl_element_type::triangle_strip;
}

void draw_text(const opengl_font* font, const string& text, float x, float y,
    float scale, const vec3f& color, float alpha) {
  auto& shader = font->program;

  set_ogl_blending(true);
  bind_program(shader);
  set_uniform(shader, "color", color);
  set_uniform(shader, "alpha", alpha);

  float aspect = 1.0f;
  {
    auto size = get_ogl_viewport_size();
    aspect    = float(size.x) / size.y;
  }
  scale /= font->size;

  for (int i = 0; i < text.size(); ++i) {
    auto& ch = font->characters[(int)text[i]];

    float xpos = x + ch.bearing.x * scale;
    float ypos = y - (ch.size.y - ch.bearing.y) * scale;

    float width  = ch.size.x * scale;
    float height = ch.size.y * scale;

    auto center  = vec2f{xpos, ypos};
    auto scaling = vec2f{width / aspect, height};

    set_uniform(shader, "center", center);
    set_uniform(shader, "scale", scaling);

    set_uniform(shader, "char_texture", ch.texture, 0);
    draw_shape(font->quad);

    // Advance cursors for next glyph (advance is number of 1/64 pixels)
    x += (ch.advance >> 6) * scale / aspect;
  }
}

}  // namespace yocto
