#pragma once

#include <yocto/yocto_mesh.h>
#include <yocto/yocto_shape.h>

#include <unordered_set>

using namespace yocto;
using std::unordered_set;

struct bool_mesh {
  vector<vec3i>        triangles   = {};
  vector<vec3i>        adjacencies = {};
  vector<vec3f>        positions   = {};
  vector<vec3f>        normals     = {};
  dual_geodesic_solver dual_solver = {};
};

struct mesh_polygon {
  vector<mesh_point> points = {};
  mesh_path          path   = {};
};

inline bool is_closed(const mesh_polygon& polygon) {
  if (polygon.points.size() < 3) return false;
  return (polygon.points.front().face == polygon.points.back().face) &&
         (polygon.points.front().uv == polygon.points.back().uv);
}

inline void update_mesh_polygon(mesh_polygon& polygon, const mesh_path& path) {
  if (!polygon.path.points.size())
    polygon.path.points.push_back(polygon.points.front());

  polygon.path.points.insert(
      polygon.path.points.end(), path.points.begin() + 1, path.points.end());
}

inline vector<int> polygon_strip(const mesh_polygon& polygon) {
  auto strip = vector<int>(polygon.path.points.size());
  for (auto i = 0; i < polygon.path.points.size(); i++)
    strip[i] = polygon.path.points[i].face;
  return strip;
}

inline vector<vec2i> polygon_segments_from_face(
    const mesh_polygon& polygon, int face) {
  auto segments = vector<vec2i>();
  for (auto i = 1; i < polygon.path.points.size(); i++)
    if (polygon.path.points[i].face == face)
      segments.push_back(vec2i{i - 1, i});
  return segments;
}

inline vector<int> intersect_polygons(
    const mesh_polygon& left, const mesh_polygon right) {
  auto left_faces  = polygon_strip(left);
  auto right_faces = polygon_strip(right);

  std::sort(left_faces.begin(), left_faces.end());
  std::sort(right_faces.begin(), right_faces.end());

  auto intersections = vector<int>();

  auto i = 0;
  auto j = 0;
  while (i < left_faces.size() && j < right_faces.size()) {
    if (left_faces[i] == right_faces[j]) {
      intersections.push_back(left_faces[i]);
      i++;
      j++;
    } else if (left_faces[i] < right_faces[j])
      i++;
    else
      j++;
  }

  return intersections;
}

inline bool_mesh init_mesh(const generic_shape* shape) {
  auto mesh        = bool_mesh{};
  mesh.triangles   = shape->triangles;
  mesh.positions   = shape->positions;
  mesh.normals     = shape->normals;
  mesh.adjacencies = face_adjacencies(mesh.triangles);

  // Fit shape in [-1, 1]^3
  auto bbox = invalidb3f;
  for (auto& pos : mesh.positions) bbox = merge(bbox, pos);
  for (auto& pos : mesh.positions) pos = (pos - center(bbox)) / max(size(bbox));

  mesh.dual_solver = make_dual_geodesic_solver(
      mesh.triangles, mesh.positions, mesh.adjacencies);
  return mesh;
}

inline geodesic_path compute_geodesic_path(
    const bool_mesh& mesh, const mesh_point& start, const mesh_point& end) {
  auto path = geodesic_path{};
  if (start.face == end.face) {
    path.start = start;
    path.end   = end;
    path.strip = {start.face};
    return path;
  }

  auto strip = strip_on_dual_graph(
      mesh.dual_solver, mesh.triangles, mesh.positions, end.face, start.face);
  path = shortest_path(
      mesh.triangles, mesh.positions, mesh.adjacencies, start, end, strip);
  return path;
}

// TODO(giacomo): Expose this function in yocto_mesh.h
static int find_in_vec(const vec3i& vec, int x) {
  for (auto i = 0; i < 3; i++)
    if (vec[i] == x) return i;
  return -1;
}
// TODO(giacomo): Expose this function in yocto_mesh.h
inline int find_adjacent_triangle(
    const vec3i& triangle, const vec3i& adjacent) {
  for (int i = 0; i < 3; i++) {
    auto k = find_in_vec(adjacent, triangle[i]);
    if (k != -1) {
      if (find_in_vec(adjacent, triangle[mod3(i + 1)]) != -1) {
        return i;
      } else {
        return mod3(i + 2);
      }
    }
  }
  // assert(0 && "input triangles are not adjacent");
  return -1;
}

struct mesh_segment {
  vec2f start = {};
  vec2f end   = {};
  int   face  = -1;
};

inline vector<mesh_segment> mesh_segments(const vector<vec3i>& triangles,
    const vector<int>& strip, const vector<float>& lerps,
    const mesh_point& start, const mesh_point& end) {
  auto result = vector<mesh_segment>(strip.size());
  
  for (int i = 0; i < strip.size(); ++i) {
    vec2f start_uv;
    if (i == 0) {
      start_uv = start.uv;
    } else {
      vec2f uvw[3] = {{0, 0}, {1, 0}, {0, 1}};
      auto  k      = find_adjacent_triangle(
          triangles[strip[i]], triangles[strip[i - 1]]);
      auto a   = uvw[k];
      auto b   = uvw[mod3(k + 1)];
      start_uv = lerp(a, b, lerps[i - 1]);
    }

    vec2f end_uv;
    if (i == strip.size() - 1) {
      end_uv = end.uv;
    } else {
      vec2f uvw[3] = {{0, 0}, {1, 0}, {0, 1}};
      auto  k      = find_adjacent_triangle(
          triangles[strip[i]], triangles[strip[i + 1]]);
      auto a   = uvw[k];
      auto b   = uvw[mod3(k + 1)];
      start_uv = lerp(a, b, lerps[i]);
    }

    result[i] = {start_uv, end_uv, strip[i]};
  }
  return result;
}
