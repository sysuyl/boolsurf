#include <yocto_gui/yocto_font.h>
#include <yocto_gui/yocto_shade.h>

inline void set_patch_shape(
    shade_shape* shape, const bool_mesh& mesh, const vector<int>& faces) {
  auto positions = vector<vec3f>(faces.size() * 3);
  for (int i = 0; i < faces.size(); i++) {
    auto [a, b, c]       = mesh.triangles[faces[i]];
    positions[3 * i + 0] = mesh.positions[a];
    positions[3 * i + 1] = mesh.positions[b];
    positions[3 * i + 2] = mesh.positions[c];
  }
  set_positions(shape, positions);
  set_instances(shape, {});
  shape->shape->elements = ogl_element_type::triangles;
}

inline void set_polygon_shape(shade_scene* scene, const bool_mesh& mesh,
    mesh_polygon& polygon, int index) {
  if (!polygon.polyline_shape) {
    polygon.polyline_shape             = add_instance(scene);
    polygon.polyline_shape->material   = add_material(scene);
    polygon.polyline_shape->depth_test = ogl_depth_test::always;
  }
  if (!polygon.polyline_shape->shape) {
    polygon.polyline_shape->shape = add_shape(scene);
  }
  polygon.polyline_shape->material->color = get_color(index);

  if (polygon.segments.empty()) return;
  auto positions = vector<vec3f>(polygon.segments.size() + 1);
  for (int i = 0; i < polygon.segments.size(); i++) {
    auto& segment = polygon.segments[i];
    positions[i]  = eval_position(mesh, {segment.face, segment.start});
  }
  {
    auto& segment    = polygon.segments.back();
    positions.back() = eval_position(mesh, {segment.face, segment.end});
  }

  set_positions(polygon.polyline_shape->shape, positions);
  polygon.polyline_shape->shape->shape->elements = ogl_element_type::line_strip;
  set_instances(polygon.polyline_shape->shape, {}, {});
}
// inline void draw_triangulation(
//     ogl_texture* texture, int face, vec2i size = {2048, 2048}) {
//   auto& triangles = debug_triangles[face];
//   auto& positions = debug_nodes[face];
//   auto& indices   = debug_indices[face];
//   if (positions.empty()) return;

//   static opengl_font* font = nullptr;
//   if (!font) {
//     font = new opengl_font{};
//     init_font(font, "data/Menlo-Regular.ttf", 100);
//   }

//   static ogl_shape* faces = nullptr;
//   if (!faces) {
//     faces = new ogl_shape{};
//   }
//   set_vertex_buffer(faces, positions, 0);
//   set_index_buffer(faces, triangles);
//   faces->elements = ogl_element_type::triangles;

//   static ogl_shape* edges = nullptr;
//   if (!edges) {
//     edges = new ogl_shape{};
//   }
//   set_vertex_buffer(edges, positions, 0);
//   auto lines = vector<vec2i>();
//   for (auto& t : triangles) {
//     lines.push_back({t.x, t.y});
//     lines.push_back({t.y, t.z});
//     lines.push_back({t.z, t.x});
//   }
//   set_index_buffer(edges, lines);
//   edges->elements = ogl_element_type::lines;

//   static ogl_shape* points = nullptr;
//   if (!points) {
//     points = new ogl_shape{};
//   }
//   set_vertex_buffer(points, positions, 0);
//   points->elements   = ogl_element_type::points;
//   points->point_size = 5;

//   static ogl_program* program = nullptr;
//   if (!program) {
//     program          = new ogl_program{};
//     auto vertex_code = R"(
//       #version 330
//       layout(location = 0) in vec2 positions;

//       void main() {
//         vec2 position = (positions - vec2(0.5)) * 1.5;
//         gl_Position = vec4(position, 0, 1);
//       }
//     )";

//     auto fragment_code = R"(
//       #version 330
//       out vec4 frag_color;
//       uniform vec3 color = vec3(1,1,1);
//       uniform float alpha = 1;
//       void main() {
//         frag_color = vec4(color, alpha);
//     })";

//     set_program(program, vertex_code, fragment_code, true);
//   }

//   if (texture->texture_id == 0 || size != texture->size) {
//     set_texture(texture, size, 4, (float*)nullptr, true, true, false, false);
//   }

//   static ogl_framebuffer* framebuffer = nullptr;
//   if (!framebuffer) {
//     framebuffer = new ogl_framebuffer{};
//     set_framebuffer(framebuffer, size);
//     set_framebuffer_texture(framebuffer, texture, 0);
//   }

//   bind_framebuffer(framebuffer);
//   set_ogl_viewport(size);
//   clear_ogl_framebuffer({0, 0, 0.1, 1}, true);
//   set_ogl_depth_test(ogl_depth_test::always);

//   bind_program(program);
//   set_uniform(program, "alpha", 1.0f);
//   set_uniform(program, "color", vec3f{1, 1, 1});
//   draw_shape(edges);

//   set_uniform(program, "color", vec3f{1, 0, 0});
//   draw_shape(points);

//   set_uniform(program, "alpha", 0.5f);
//   set_uniform(program, "color", vec3f{0.5, 0.5, 0.5});
//   draw_shape(faces);

//   auto text_size = 0.0004;
//   for (int i = 0; i < positions.size(); i++) {
//     auto text   = to_string(i);
//     auto coords = (positions[i] - vec2f{0.5f, 0.5f}) * 1.5f;
//     draw_text(font, text, coords.x, coords.y, text_size, {0.8, 0.4, 0.1});
//   }

//   {
//     auto lines = vector<string>{};
//     auto color = vec3f{0.8, 0.4, 0.1};

//     lines += "face: "s + to_string(face);
//     lines += ""s;
//     lines += "triangles"s;
//     for (int i = 0; i < triangles.size(); i++) {
//       auto [a, b, c] = triangles[i];
//       auto text      = "(" + to_string(a) + ", " + to_string(b) + ", " +
//                   to_string(c) + ")";
//       lines += text;
//     }
//     lines += ""s;

//     lines += "positions"s;
//     for (int i = 0; i < positions.size(); i++) {
//       auto [a, b] = positions[i];
//       auto text   = "[" + to_string(indices[i]) + "] ";
//       text += to_string(i) + ": (" + to_string(a) + ", " + to_string(b) +
//       ")"; lines += text;
//     }

//     draw_text(font, lines, 0.1, 0.9, text_size, color);
//   }

//   unbind_framebuffer();
// }

// inline void save_triangulation(const string& filename, int face) {
//   auto texture = new ogl_texture{};
//   draw_triangulation(texture, face);
//   auto img = get_texture(texture);

//   // flip verticllay
//   for (int y = 0; y < img.imsize().y / 2; y++)
//     for (int x = 0; x < img.imsize().x; x++)
//       std::swap(img[{x, y}], img[{x, img.imsize().y - y - 1}]);

//   auto error = ""s;
//   if (!save_image(filename, img, error)) {
//     printf("%s: %s\n", __FUNCTION__, error.c_str());
//   }
// }

[[nodiscard]] shade_instance* draw_sphere(shade_scene* scene,
    const bool_mesh& mesh, shade_material* material, const vector<vec3f>& pos,
    float dim) {
  auto sphere = make_sphere(4, dim);

  auto shape = add_shape(scene, {}, {}, {}, sphere.quads, sphere.positions,
      sphere.normals, sphere.texcoords, {});
  set_instances(shape, pos);
  return add_instance(scene, identity3x4f, shape, material, false);
}

[[nodiscard]] shade_instance* draw_segment(shade_scene* scene,
    const bool_mesh& mesh, shade_material* material, const vec3f& start,
    const vec3f& end, float radius = 0.0006f) {
  auto cylinder = make_uvcylinder({4, 1, 1}, {radius, 1});
  for (auto& p : cylinder.positions) {
    p.z = p.z * 0.5 + 0.5;
  }

  auto shape = add_shape(scene);
  set_quads(shape, cylinder.quads);
  set_positions(shape, cylinder.positions);
  set_normals(shape, cylinder.normals);
  set_texcoords(shape, cylinder.texcoords);
  set_instances(shape, {start}, {end});
  return add_instance(scene, identity3x4f, shape, material, false);
}

#if 0

[[nodiscard]] shade_instance* draw_mesh_point(shade_scene* scene,
    const bool_mesh& mesh, shade_material* material, const mesh_point& point,
    float dim) {
  auto pos = eval_position(mesh.triangles, mesh.positions, point);
  return draw_sphere(scene, mesh, material, {pos}, dim);
}

void update_path_shape(shade_shape* shape, const bool_mesh& mesh,
    const geodesic_path& path, float radius, float offset = 0,
    bool thin = true) {
  auto positions = path_positions(
      path, mesh.triangles, mesh.positions, mesh.adjacencies);

  if (thin) {
    set_positions(shape, positions);
    shape->shape->elements = ogl_element_type::line_strip;
    set_instances(shape, {}, {});
    return;
  }

  auto froms = vector<vec3f>();
  auto tos   = vector<vec3f>();
  froms.reserve(positions.size() - 1);
  tos.reserve(positions.size() - 1);
  for (int i = 0; i < positions.size() - 1; i++) {
    auto from = positions[i];
    auto to   = positions[i + 1];
    if (from == to) continue;
    froms.push_back(from);
    tos.push_back(to);
  }

  auto cylinder = make_uvcylinder({4, 1, 1}, {radius, 1});
  for (auto& p : cylinder.positions) {
    p.z = p.z * 0.5 + 0.5;
  }

  set_quads(shape, cylinder.quads);
  set_positions(shape, cylinder.positions);
  set_normals(shape, cylinder.normals);
  set_texcoords(shape, cylinder.texcoords);
  set_instances(shape, froms, tos);
}

void update_path_shape(shade_shape* shape, const bool_mesh& mesh,
    const mesh_path& path, float radius, float offset = 0, bool thin = false) {
  auto positions = vector<vec3f>(path.points.size());
  for (int i = 0; i < positions.size(); i++) {
    positions[i] = eval_position(
        mesh.triangles, mesh.positions, path.points[i]);
  }

  if (offset > 0) {
    // auto mesh_points = convert_mesh_path(mesh.triangles, mesh.adjacencies,
    // path.strip, path.lerps, path.start, path.end);
    auto pos              = positions;
    int  num_subdivisions = 8;
    for (int i = 0; i < num_subdivisions; i++) {
      auto pos = positions;
      for (int i = 0; i < pos.size(); i++) {
        auto a       = (i - 1 + (int)pos.size()) % pos.size();
        auto c       = (i + 1) % pos.size();
        positions[i] = (positions[a] + positions[c]) / 2;
      }
    }
    for (int i = 0; i < pos.size(); i++) {
      auto a  = (i - 1 + (int)pos.size()) % pos.size();
      auto b  = i;
      auto c  = (i + 1) % pos.size();
      auto n0 = eval_normal(mesh.triangles, mesh.normals, path.points[b]);
      auto v  = pos[b] - pos[a];
      auto w  = pos[c] - pos[b];
      positions[i] += offset * normalize(cross(n0, v)) / 2;
      positions[i] += offset * normalize(cross(n0, w)) / 2;
    }
  }

  if (thin) {
    set_positions(shape, positions);
    shape->shape->elements = ogl_element_type::line_strip;
    set_instances(shape, {}, {});
    return;
  }

  auto froms = vector<vec3f>();
  auto tos   = vector<vec3f>();
  froms.reserve(positions.size() - 1);
  tos.reserve(positions.size() - 1);
  for (int i = 0; i < positions.size() - 1; i++) {
    auto from = positions[i];
    auto to   = positions[i + 1];
    if (from == to) continue;
    froms.push_back(from);
    tos.push_back(to);
  }

  auto cylinder = make_uvcylinder({4, 1, 1}, {radius, 1});
  for (auto& p : cylinder.positions) {
    p.z = p.z * 0.5 + 0.5;
  }

  set_quads(shape, cylinder.quads);
  set_positions(shape, cylinder.positions);
  set_normals(shape, cylinder.normals);
  set_texcoords(shape, cylinder.texcoords);
  set_instances(shape, froms, tos);
}


[[nodiscard]] shade_instance* draw_mesh_segment(shade_scene* scene,
    const bool_mesh& mesh, shade_material* material,
    const mesh_segment& segment, float radius = 0.0012f) {
  auto start = mesh_point{segment.face, segment.start};
  auto end   = mesh_point{segment.face, segment.end};

  draw_mesh_point(scene, mesh, material, start, radius);
  draw_mesh_point(scene, mesh, material, end, radius);

  auto pos_start = eval_position(mesh.triangles, mesh.positions, start);
  auto pos_end   = eval_position(mesh.triangles, mesh.positions, end);

  return draw_segment(scene, mesh, material, pos_start, pos_end, radius / 2);
}

[[nodiscard]] vector<shade_instance*> draw_arrangement(shade_scene* scene,
    const bool_mesh& mesh, const vector<shade_material*>& material,
    const vector<mesh_point>& points, vector<cell_polygon>& cells) {
  auto instances = vector<shade_instance*>{};
  for (auto p = 0; p < cells.size(); p++) {
    auto& polygon = cells[p];
    auto  mat     = material[p % material.size()];
    auto  path    = mesh_path{};
    for (auto n = 0; n < polygon.points.size() - 1; n++) {
      auto& start = points[polygon.points[n]];
      auto& end   = points[polygon.points[n + 1]];

      auto geo_path = compute_geodesic_path(mesh, start, end);
      append(path.points,
          convert_mesh_path(mesh.triangles, mesh.adjacencies, geo_path.strip,
              geo_path.lerps, geo_path.start, geo_path.end)
              .points);

      // auto segments = mesh_segments(mesh.triangles, geo_path.strip,
      //    geo_path.lerps, geo_path.start, geo_path.end);
      // update_mesh_polygon(polygon, segments);
    }
    auto shape = add_shape(scene);
    // TODO: Make this proportional to avg_edge_length
    float offset = 0.002f;
    update_path_shape(shape, mesh, path, 0.0010f, offset);
    instances.push_back(add_instance(scene, identity3x4f, shape, mat, false));
  }
  return instances;
}
#endif
