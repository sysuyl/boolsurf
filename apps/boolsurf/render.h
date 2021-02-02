#include <yocto_gui/yocto_shade.h>

[[nodiscard]] shade_instance* draw_sphere(shade_scene* scene,
    const bool_mesh& mesh, shade_material* material, const vector<vec3f>& pos,
    float dim) {
  auto sphere = make_sphere(4, dim);

  auto shape = add_shape(scene, {}, {}, {}, sphere.quads, sphere.positions,
      sphere.normals, sphere.texcoords, {});
  set_instances(shape, pos);
  return add_instance(scene, identity3x4f, shape, material, false);
}

[[nodiscard]] shade_instance* draw_mesh_point(shade_scene* scene,
    const bool_mesh& mesh, shade_material* material, const mesh_point& point,
    float dim) {
  auto pos = eval_position(mesh.triangles, mesh.positions, point);
  return draw_sphere(scene, mesh, material, {pos}, dim);
}

[[nodiscard]] shade_instance* draw_path(shade_scene* scene,
    const bool_mesh& mesh, shade_material* material, const geodesic_path& path,
    float radius) {
  auto shape = add_shape(scene);
  update_path_shape(shape, mesh, path, radius);
  add_instance(scene, identity3x4f, shape, material, false);
  return add_instance(scene, identity3x4f, shape, material, false);
}

[[nodiscard]] shade_instance* draw_intersections(shade_scene* scene,
    const bool_mesh& mesh, shade_material* material, const vector<int>& isecs) {
  auto pos = vector<vec3f>(isecs.size());
  for (auto i = 0; i < isecs.size(); i++) {
    auto v = mesh.triangles[isecs[i]];
    pos[i] = (mesh.positions[v.x] + mesh.positions[v.y] + mesh.positions[v.z]) /
             3.0f;
  }

  return draw_sphere(scene, mesh, material, pos, 0.0015f);
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

[[nodiscard]] void set_patch_shape(shade_shape* shape, const bool_mesh& mesh,
    const vector<int>& faces, const float distance) {
  auto positions = vector<vec3f>(faces.size() * 3);
  for (int i = 0; i < faces.size(); i++) {
    auto [a, b, c]       = mesh.triangles[faces[i]];
    positions[3 * i + 0] = mesh.positions[a] + distance * mesh.normals[a];
    positions[3 * i + 1] = mesh.positions[b] + distance * mesh.normals[b];
    positions[3 * i + 2] = mesh.positions[c] + distance * mesh.normals[c];
  }
  set_positions(shape, positions);
  set_instances(shape, {});
  shape->shape->elements = ogl_element_type::triangles;
}