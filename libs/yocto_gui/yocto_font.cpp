#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#endif
#include "ext/glad/glad.h"
//
#include <ft2build.h>
#include FT_FREETYPE_H
#include FT_OUTLINE_H

#include "yocto_font.h"
using namespace yocto;

// From opengl tutorial <3
// https://learnopengl.com/In-Practice/Text-Rendering

#if 1
auto vertex = R"(
  #version 330 core
  layout(location = 0) in vec2 vertex;
  out vec2     TexCoords;
  uniform vec2 scale    = vec2(1, 1);
  uniform vec2 center = vec2(0.0, 0.0);
  float        aspect = 1;

  void main() {
    TexCoords = vec2(vertex.x, 1 - vertex.y);
    vec2 p = vertex;
    p *= scale;
    p.y *= -1;
    p += center;
    gl_Position = vec4(p.x, p.y, 0, 1);
  }
)";

auto fragment = R"(
  #version 330 core
  in vec2 TexCoords;
  out vec4 out_color;

  uniform sampler2D text;
  uniform vec3 color = vec3(1, 1, 1);
  uniform float alpha = 1;

  void main()
  {    
      float sampled = texture(text, TexCoords).r;
      out_color = vec4(color, alpha * sampled);
  }
)";
#else
auto vertex = R"(
  #version 330 core
  layout (location = 0) in vec4 vertex; // <vec2 pos, vec2 tex>
  out vec2 TexCoords;

  uniform mat4 projection;

  void main()
  {
      gl_Position = projection * vec4(vertex.xy, 0.0, 1.0);
      TexCoords = vertex.zw;
  } 
)";

auto fragment = R"(
  #version 330 core
  in vec2 TexCoords;
  out vec4 color;

  uniform sampler2D text;
  uniform vec3 textColor;

  void main()
  {    
      vec4 sampled = vec4(1.0, 1.0, 1.0, texture(text, TexCoords).r);
      color = vec4(textColor, 1.0) * sampled;
  })";
#endif

void init_glfont(opengl_font& font, const string& filename, float size) {
  auto error    = string{};
  auto errorlog = string{};
  if (!set_program(font.program, vertex, fragment, error, errorlog)) {
    printf("%s: %s\n", __FUNCTION__, error.c_str());
    printf("%s: %s\n", __FUNCTION__, errorlog.c_str());
  }

  font.size = size;

  // FreeType
  FT_Library ft;
  if (FT_Init_FreeType(&ft))
    printf("ERROR::FREETYPE: Could not init FreeType Library\n");

  // Load font as face
  FT_Face face;
  if (FT_New_Face(ft, filename.c_str(), 0, &face))
    printf("ERROR::FREETYPE: Failed to load font %s\n", filename.c_str());

  // Set size to load glyphs as
  FT_Set_Pixel_Sizes(face, 0, (int)font.size);

  // Disable byte-alignment restriction
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

  // Load first 128 characters of ASCII set
  for (int c = 0; c < 128; c++) {
    // Load character glyph
    if (FT_Load_Char(face, c, FT_LOAD_RENDER)) {
      printf("%s: ERROR::FREETYTPE: Failed to load Glyph\n", __FUNCTION__);
      continue;
    }

    // Add entry to hash map.
    auto& character = font.characters[(char)c];
    character.size  = {
        int(face->glyph->bitmap.width), int(face->glyph->bitmap.rows)};
    character.bearing = {face->glyph->bitmap_left, face->glyph->bitmap_top};
    character.advance = (GLuint)face->glyph->advance.x;

    // Generate texture
    assert(glGetError() == GL_NO_ERROR);
    auto& texture = character.texture_id;
    glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, face->glyph->bitmap.width,
        face->glyph->bitmap.rows, 0, GL_RED, GL_UNSIGNED_BYTE,
        face->glyph->bitmap.buffer);

    // Set texture options
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    assert(glGetError() == GL_NO_ERROR);
  }
  glBindTexture(GL_TEXTURE_2D, 0);

  // Destroy FreeType
  FT_Done_Face(face);
  FT_Done_FreeType(ft);

  // clang-format off
  static const auto positions = vector<vec3f>{
    {0, 0}, {1, 0}, {1, 1}, {0, 1},
  };
  static const auto triangles = vector<vec3i>{
    {0, 1, 3}, {3, 2, 1}
  };
  // clang-format on
  set_vertex_buffer(font.quad, positions, 0);
  set_vertex_buffer(font.quad, positions, 1);
  set_index_buffer(font.quad, triangles);
}

void draw_glfont_dicoso(const opengl_font& font, const string& text, float x,
    float y, float scale, const vec3f& color) {
  unsigned int VAO, VBO = 0;
  glGenVertexArrays(1, &VAO);
  glGenBuffers(1, &VBO);
  glBindVertexArray(VAO);
  glBindBuffer(GL_ARRAY_BUFFER, VBO);
  glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 6 * 4, NULL, GL_DYNAMIC_DRAW);
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(float), 0);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);

  bind_program(font.program);
  set_uniform(font.program, "textColor", color);
  glActiveTexture(GL_TEXTURE0);
  glBindVertexArray(VAO);

  // iterate through all characters
  for (int i = 0; i < text.size(); ++i) {
    const auto& ch   = font.characters.at(text[i]);
    float       xpos = x + ch.bearing.x * scale;
    float       ypos = y - (ch.size.y - ch.bearing.y) * scale;

    float w = ch.size.x * scale;
    float h = ch.size.y * scale;
    // update VBO for each character
    float vertices[6][4] = {{xpos, ypos + h, 0.0f, 0.0f},
        {xpos, ypos, 0.0f, 1.0f}, {xpos + w, ypos, 1.0f, 1.0f},

        {xpos, ypos + h, 0.0f, 0.0f}, {xpos + w, ypos, 1.0f, 1.0f},
        {xpos + w, ypos + h, 1.0f, 0.0f}};
    // render glyph texture over quad
    glBindTexture(GL_TEXTURE_2D, ch.texture_id);
    // update content of VBO memory
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(vertices), vertices);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    // render quad
    glDrawArrays(GL_TRIANGLES, 0, 6);
    // now advance cursors for next glyph (note that advance is number of 1/64
    // pixels)
    x += (ch.advance >> 6) *
         scale;  // bitshift by 6 to get value in pixels (2^6 = 64)
  }
  glBindVertexArray(0);
  glBindTexture(GL_TEXTURE_2D, 0);
}

void draw_glfont(const opengl_font& font, const string& text, float x, float y,
    float scale, const vec3f& color) {
  auto shader = font.program;

  //  auto bottom_left = frame.o - frame.x - frame.y;
  // glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  set_ogl_blending(true);
  bind_program(shader);
  set_uniform(shader, "color", color);
  //  set_uniform(shader, "alpha", alpha);
  assert(glGetError() == GL_NO_ERROR);
  glActiveTexture(GL_TEXTURE0);

  for (int i = 0; i < text.size(); ++i) {
    auto& ch = font.characters.at(text[i]);

    printf("\ncharacter[%d]: %c\n", i, text[i]);
    printf("size: %d, %d\n", ch.size.x, ch.size.y);
    printf("bearing: %d, %d\n", ch.bearing.x, ch.bearing.y);
    printf("advance: %d\n", (int)ch.advance);

    float xpos = x + ch.bearing.x * scale;
    float ypos = y - (ch.size.y - ch.bearing.y) * scale;

    float width  = ch.size.x * scale;
    float height = ch.size.y * scale;

    // float vertices[6][4] = {{xpos, ypos + h}, {xpos, ypos}, {xpos + w, ypos},
    //     {xpos, ypos + h}, {xpos + w, ypos}, {xpos + w, ypos + h}};
    auto center  = vec2f{xpos, ypos};
    auto scaling = vec2f{width, height};

    // auto scale = vec2f{width * 0.5f, height * 0.5f} * dxy;
    set_uniform(shader, "center", center);
    set_uniform(shader, "scale", scaling);
    // auto smat = mat2f{{scale.x, 1}, {1, scale.y}};
    // smat.x.x /= ratio;
    // auto mat = rot * smat;
    // set_uniform(shader, "mat", mat);
    // set_uniform(shader, "aspect", ratio);

    glBindTexture(GL_TEXTURE_2D, ch.texture_id);
    draw_shape(font.quad);
    assert(glGetError() == GL_NO_ERROR);

    // Advance cursors for next glyph (advance is number of 1/64 pixels)
    // line_pos.x += (ch.advance / 64.0);
    x += (ch.advance >> 6) * scale;
  }
}
