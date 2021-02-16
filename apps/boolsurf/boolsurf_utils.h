#pragma once

#include <yocto/yocto_mesh.h>
#include <yocto/yocto_shape.h>

#include <cassert>
#include <deque>
#include <unordered_set>

#include "ext/CDT/CDT/include/CDT.h"
#include "ext/delaunator.cpp"

using namespace yocto;
using namespace std;

struct bool_mesh {
  vector<vec3i>        triangles   = {};
  vector<vec3i>        adjacencies = {};
  vector<vec3f>        positions   = {};
  vector<vec3f>        normals     = {};
  dual_geodesic_solver dual_solver = {};
  vector<vec3i>        tags        = {};
};

struct mesh_segment {
  vec2f start = {};
  vec2f end   = {};
  int   face  = -1;
};

struct mesh_polygon {
  vector<int>          points      = {};
  vector<mesh_segment> segments    = {};
  vector<int>          inner_faces = {};
  vector<int>          outer_faces = {};

  shade_instance* polyline_shape = nullptr;
  shade_instance* inner_shape    = nullptr;
  shade_instance* outer_shape    = nullptr;
};

struct hashgrid_segment {
  int polygon = -1;
  int segment = -1;

  vec2f start = {};
  vec2f end   = {};
};

struct intersection {
  int   vertex_id = -1;
  float lerp      = -1.0f;
};

struct cell_polygon {
  vector<int>          points    = {};
  vector<mesh_segment> segments  = {};
  vector<int>          embedding = {};
};

// Vector append and concatenation
template <typename T>
inline void operator+=(vector<T>& a, const vector<T>& b) {
  a.insert(a.end(), b.begin(), b.end());
}
template <typename T>
inline void operator+=(vector<T>& a, const T& b) {
  a.push_back(b);
}
template <typename T>
inline vector<T> operator+(const vector<T>& a, const vector<T>& b) {
  auto c = a;
  c += b;
  return b;
}

inline vec3f eval_position(const bool_mesh& mesh, const mesh_point& point) {
  return eval_position(mesh.triangles, mesh.positions, point);
}

// (marzia) Not Used
inline bool is_closed(const mesh_polygon& polygon) {
  if (polygon.points.size() < 3) return false;
  return (polygon.points.front() == polygon.points.back());
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

template <class T>
inline int find_idx(const vector<T>& vec, const T& x) {
  for (auto i = 0; i < vec.size(); i++)
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

inline vector<mesh_segment> mesh_segments(const vector<vec3i>& triangles,
    const vector<int>& strip, const vector<float>& lerps,
    const mesh_point& start, const mesh_point& end) {
  auto result = vector<mesh_segment>{};
  result.reserve(strip.size());

  for (int i = 0; i < strip.size(); ++i) {
    vec2f start_uv;
    if (i == 0) {
      start_uv = start.uv;
    } else {
      vec2f uvw[3] = {{0, 0}, {1, 0}, {0, 1}};
      auto  k      = find_adjacent_triangle(
          triangles[strip[i]], triangles[strip[i - 1]]);
      auto a   = uvw[mod3(k)];
      auto b   = uvw[mod3(k + 1)];
      start_uv = lerp(a, b, 1 - lerps[i - 1]);
    }

    vec2f end_uv;
    if (i == strip.size() - 1) {
      end_uv = end.uv;
    } else {
      vec2f uvw[3] = {{0, 0}, {1, 0}, {0, 1}};
      auto  k      = find_adjacent_triangle(
          triangles[strip[i]], triangles[strip[i + 1]]);
      auto a = uvw[k];
      auto b = uvw[mod3(k + 1)];
      end_uv = lerp(a, b, lerps[i]);
    }
    if (start_uv == end_uv) continue;
    result.push_back({start_uv, end_uv, strip[i]});
  }
  return result;
}

// From yocto_mesh.h + small update
inline vec2f intersect_segments(const vec2f& start1, const vec2f& end1,
    const vec2f& start2, const vec2f& end2) {
  if (end1 == start2) return zero2f;
  if (end2 == start1) return one2f;
  if (start2 == start1) return zero2f;
  if (end2 == end1) return one2f;

  auto a = end1 - start1;    // direction of line a
  auto b = start2 - end2;    // direction of line b, reversed
  auto d = start2 - start1;  // right-hand side

  auto det = a.x * b.y - a.y * b.x;
  if (det == 0) return {-1, -1};

  auto r = (d.x * b.y - d.y * b.x) / det;
  auto s = (a.x * d.y - a.y * d.x) / det;
  return {r, s};
}

struct hashgrid_polyline {
  int           polygon  = -1;
  vector<vec2f> points   = {};
  vector<int>   vertices = {};
};

inline pair<int, vec2f> add_vertex_new(
    bool_mesh& mesh, const mesh_point& point) {
  float eps = 0.00001;
  auto  uv  = point.uv;
  auto  tr  = mesh.triangles[point.face];
  if (uv.x < eps && uv.y < eps) return {tr.x, {0, 0}};
  if (uv.x > 1 - eps && uv.y < eps) return {tr.y, {1, 0}};
  if (uv.y > 1 - eps && uv.x < eps) return {tr.z, {1, 1}};
  auto vertex = (int)mesh.positions.size();
  auto pos    = eval_position(mesh.triangles, mesh.positions, point);
  mesh.positions.push_back(pos);
  return {vertex, uv};
}

inline unordered_map<int, vector<hashgrid_polyline>> compute_hashgrid(
    const vector<mesh_polygon>& polygons, const vector<vector<int>>& vertices) {
  auto hashgrid = unordered_map<int, vector<hashgrid_polyline>>{};

  for (auto polygon_id = 0; polygon_id < polygons.size(); polygon_id++) {
    auto& polygon   = polygons[polygon_id];
    int   last_face = -1;

    if (polygon.segments.empty()) continue;
    int first_face = polygon.segments[0].face;
    int s_first    = -1;

    for (auto s = 0; s < polygon.segments.size(); s++) {
      auto& segment = polygon.segments[s];

      if (segment.face == first_face && s_first == -1) continue;
      if (s_first == -1) s_first = s;

      auto& entry = hashgrid[segment.face];
      if (segment.face != last_face) {
        auto& polyline   = entry.emplace_back();
        polyline.polygon = polygon_id;

        auto ss = s - 1;
        assert(ss >= 0);
        polyline.vertices.push_back(vertices[polygon_id][ss]);
        polyline.points.push_back(segment.start);

        polyline.vertices.push_back(vertices[polygon_id][s]);
        polyline.points.push_back(segment.end);
      } else {
        auto& polyline = entry.back();
        assert(segment.end != polyline.points.back());
        polyline.points.push_back(segment.end);
        polyline.vertices.push_back(vertices[polygon_id][s]);
      }
      auto& polyline = entry.back();
      if (polyline.points.size() >= 2) {
        assert(polyline.points.back() != polyline.points.end()[-2]);
      }
      last_face = segment.face;
    }

    // Ripetiamo perche' la prima polyline non la calcoliamo al primo giro.
    assert(last_face != -1);
    for (auto s = 0; s < s_first; s++) {
      auto& segment  = polygon.segments[s];
      auto& entry    = hashgrid[segment.face];
      auto& polyline = entry.back();
      assert(segment.face == last_face);
      polyline.points.push_back(segment.end);
      polyline.vertices.push_back(vertices[polygon_id][s]);
    }
  }
  return hashgrid;
}

inline int add_vertex(bool_mesh& mesh, const mesh_point& point) {
  float eps = 0.00001;
  auto  uv  = point.uv;
  auto  tr  = mesh.triangles[point.face];
  if (uv.x < eps && uv.y < eps) return tr.x;
  if (uv.x > 1 - eps && uv.y < eps) return tr.y;
  if (uv.y > 1 - eps && uv.x < eps) return tr.z;
  auto vertex = (int)mesh.positions.size();
  auto pos    = eval_position(mesh.triangles, mesh.positions, point);
  mesh.positions.push_back(pos);
  return vertex;
}

// TODO: put in utils
template <typename T>
inline void insert(vector<T>& vec, size_t i, const T& x) {
  vec.insert(vec.begin() + i, x);
}

inline bool check_tags(const bool_mesh& mesh) {
  for (int i = 0; i < mesh.triangles.size(); i++) {
    auto face = i;
    auto tr   = mesh.triangles[face];
    if (tr == vec3i{0, 0, 0}) continue;
    for (int k = 0; k < 3; k++) {
      auto neighbor = mesh.adjacencies[face][k];
      if (neighbor == -1) continue;
      auto n0 = mesh.adjacencies[face];
      auto n1 = mesh.adjacencies[neighbor];
      auto kk = find_in_vec(mesh.adjacencies[neighbor], face);
      assert(kk != -1);

      auto tags0 = mesh.tags[face];
      auto tags1 = mesh.tags[neighbor];
      auto tag0  = tags0[k];
      auto tag1  = tags1[kk];
      assert(tag0 == -tag1);
    }
  }
  return true;
}
inline void compute_intersections(
    unordered_map<int, vector<hashgrid_polyline>>& hashgrid, bool_mesh& mesh) {
  for (auto& [face, polylines] : hashgrid) {
    // Check for polyline self interesctions
    for (auto p0 = 0; p0 < polylines.size(); p0++) {
      auto& poly = polylines[p0];

      int num_added = 0;
      for (int s0 = 0; s0 < poly.points.size() - 2; s0++) {
        auto& start0 = poly.points[s0];
        auto& end0   = poly.points[(s0 + 1) % poly.points.size()];
        for (int s1 = s0 + 2; s1 < poly.points.size(); s1++) {
          auto& start1 = poly.points[s1];
          auto& end1   = poly.points[(s1 + 1) % poly.points.size()];

          auto l = intersect_segments(start0, end0, start1, end1);
          if (l.x <= 0.0f || l.x >= 1.0f || l.y <= 0.0f || l.y >= 1.0f) {
            continue;
          }

          auto uv     = lerp(start1, end1, l.y);
          auto vertex = add_vertex(mesh, {face, uv});

          insert(poly.points, s0 + 1, uv);
          insert(poly.vertices, s0 + 1, vertex);
          insert(poly.points, s1 + 2, uv);
          insert(poly.vertices, s1 + 2, vertex);
          num_added += 1;
          s1 += 2;
        }
        s0 += num_added;
      }
    }

    // Check for intersections between different polylines
    for (auto p0 = 0; p0 < polylines.size() - 1; p0++) {
      for (auto p1 = p0 + 1; p1 < polylines.size(); p1++) {
        auto& poly0     = polylines[p0];
        auto& poly1     = polylines[p1];
        int   num_added = 0;
        for (int s0 = 0; s0 < poly0.points.size() - 1; s0++) {
          auto& start0 = poly0.points[s0];
          auto& end0   = poly0.points[(s0 + 1)];
          for (int s1 = 0; s1 < poly1.points.size() - 1; s1++) {
            auto& start1 = poly1.points[s1];
            auto& end1   = poly1.points[(s1 + 1)];
            auto  l      = intersect_segments(start0, end0, start1, end1);
            if (l.x <= 0.0f || l.x >= 1.0f || l.y <= 0.0f || l.y >= 1.0f) {
              continue;
            }

            auto uv     = lerp(start1, end1, l.y);
            auto vertex = add_vertex(mesh, {face, uv});

            insert(poly0.points, s0 + 1, uv);
            insert(poly0.vertices, s0 + 1, vertex);
            insert(poly1.points, s1 + 1, uv);
            insert(poly1.vertices, s1 + 1, vertex);
            num_added += 1;
            s1 += 1;
          }
          s0 += num_added;
        }
      }
    }
  }
}

inline vec2i make_edge_key(const vec2i& edge) {
  if (edge.x > edge.y) return {edge.y, edge.x};
  return edge;
};

inline tuple<vec2i, float> get_mesh_edge(
    const vec3i& triangle, const vec2f& uv) {
  if (uv.y == 0)
    return {vec2i{triangle.x, triangle.y}, uv.x};  // point on edge(xy)
  else if (uv.x == 0)
    return {vec2i{triangle.z, triangle.x}, 1.0f - uv.y};  // point on edge (xz)
  else if (fabs(uv.x + uv.y - 1.0f) < 0.0001)
    return {vec2i{triangle.y, triangle.z}, uv.y};  // point on edge (yz)
  else
    return {zero2i, -1};
}

inline vec2i get_edge(const vec3i& triangle, int k) {
  if (k == 0)
    return {triangle.x, triangle.y};
  else if (k == 1)
    return {triangle.y, triangle.z};
  else if (k == 2)
    return {triangle.z, triangle.x};
  else {
    assert(0);
    return {-1, -1};
  }
}

// (Previous implementation) Simple Delaunay Triangulation
inline vector<vec3i> triangulate(const vector<vec2f>& nodes) {
  auto coords = vector<double>();
  coords.reserve(nodes.size() * 2);
  for (auto& node : nodes) {
    coords.push_back(node.x);
    coords.push_back(node.y);
  }

  auto dt        = delaunator::Delaunator(coords);
  auto triangles = vector<vec3i>();
  triangles.reserve(dt.triangles.size() / 3);
  for (int i = 0; i < dt.triangles.size(); i += 3) {
    auto verts = vec3i{(int)dt.triangles[i], (int)dt.triangles[i + 2],
        (int)dt.triangles[i + 1]};

    // Check collinearity
    auto& a    = nodes[verts.x];
    auto& b    = nodes[verts.y];
    auto& c    = nodes[verts.z];
    auto  area = cross(b - a, c - a);
    if (fabs(area) < 0.00001) {
      printf("heyyyy\n");
      continue;
    }
    // if (fabs(orientation) < 0.00001) {
    //   continue;
    // }

    triangles.push_back(verts);
  }

  // // Area of whole triangle must be 1.
  // auto real_area = cross(nodes[1] - nodes[0], nodes[2] - nodes[0]);
  // assert(fabs(real_area - 1) < 0.001);

  // // Check total area.
  // auto area = 0.0f;
  // for (auto& tr : triangles) {
  //   area += cross(nodes[tr.y] - nodes[tr.x], nodes[tr.z] - nodes[tr.x]);
  // }
  // assert(fabs(area - real_area) < 0.001);

  return triangles;
}

// Constrained Delaunay Triangulation
inline vector<vec3i> constrained_triangulation(
    vector<vec2f> nodes, const vector<vec2i>& edges) {
  for (auto& n : nodes) n *= 1e9;

  auto cdt = CDT::Triangulation<double>(
      CDT::FindingClosestPoint::ClosestRandom);
  cdt.insertVertices(
      nodes.begin(), nodes.end(),
      [](const vec2f& point) -> double { return point.x; },
      [](const vec2f& point) -> double { return point.y; });
  cdt.insertEdges(
      edges.begin(), edges.end(), [](const vec2i& edge) { return edge.x; },
      [](const vec2i& edge) { return edge.y; });

  cdt.eraseOuterTriangles();
  auto triangles = vector<vec3i>();
  triangles.reserve(cdt.triangles.size());

  for (auto& tri : cdt.triangles) {
    auto verts = vec3i{
        (int)tri.vertices[0], (int)tri.vertices[1], (int)tri.vertices[2]};

    // Check collinearity
    auto& a           = nodes[verts.x];
    auto& b           = nodes[verts.y];
    auto& c           = nodes[verts.z];
    auto  orientation = cross(b - a, c - b);
    if (fabs(orientation) < 0.00001) {
      printf("Collinear\n");
      continue;
    }
    triangles.push_back(verts);
  }
  return triangles;
}

inline void update_face_edgemap(unordered_map<vec2i, vec2i>& face_edgemap,
    const vec2i& edge, const int face) {
  auto key = make_edge_key(edge);

  auto it = face_edgemap.find(key);
  if (it == face_edgemap.end()) {
    //   auto& faces = face_edgemap[key];
    // }
    //   if (faces.x == -1) {
    //     assert(faces.y == -1);
    //     faces.x = face;
    face_edgemap.insert(it, {key, {face, -1}});
  } else {
    // assert(faces.y == -1);
    it->second.y = face;
  }
}

inline vector<vec3i> compute_face_tags(
    const bool_mesh& mesh, const vector<mesh_polygon>& polygons) {
  auto tags = vector<vec3i>(mesh.triangles.size(), zero3i);
  for (auto p = 1; p < polygons.size(); p++) {
    for (auto f : polygons[p].inner_faces) {
      for (auto k = 0; k < 3; k++) {
        if (tags[f][k] == 0) {
          tags[f][k] = -p;
          break;
        }
      }
    }

    for (auto f : polygons[p].outer_faces) {
      for (auto k = 0; k < 3; k++) {
        if (tags[f][k] == 0) {
          tags[f][k] = p;
          break;
        }
      }
    }
  }
  return tags;
}

template <typename F>
vector<int> flood_fill(const bool_mesh& mesh, const vector<int>& start,
    const int polygon, F&& check) {
  auto visited = vector<bool>(mesh.adjacencies.size(), false);

  auto result = vector<int>();
  auto stack  = start;

  while (!stack.empty()) {
    auto face = stack.back();
    stack.pop_back();

    if (visited[face]) continue;
    visited[face] = true;

    result.push_back(face);

    for (auto neighbor : mesh.adjacencies[face]) {
      if (neighbor < 0 || visited[neighbor])
        continue;
      else if (check(face, -polygon) && check(neighbor, -polygon))
        // Check if "face" is not inner and "neighbor" is outer
        stack.push_back(neighbor);
      else if (check(neighbor, polygon))
        stack.push_back(neighbor);
    }
  }

  return result;
}

template <typename F>
vector<int> flood_fill(
    const bool_mesh& mesh, const vector<int>& start, F&& check) {
  auto visited = vector<bool>(mesh.adjacencies.size(), false);

  auto result = vector<int>();
  auto stack  = start;

  while (!stack.empty()) {
    auto face = stack.back();
    stack.pop_back();

    if (visited[face]) continue;
    visited[face] = true;

    result.push_back(face);

    for (auto neighbor : mesh.adjacencies[face]) {
      if (neighbor < 0 || visited[neighbor]) continue;
      if (check(face, neighbor)) stack.push_back(neighbor);
    }
  }

  return result;
}

static auto debug_result  = vector<int>();
static auto debug_visited = vector<bool>{};
static auto debug_stack   = vector<int>{};
static auto debug_restart = true;

template <typename F>
void flood_fill_debug(
    const bool_mesh& mesh, const vector<int>& start, F&& check) {
  int face = -1;
  if (debug_stack.empty()) {
    debug_restart = true;
    return;
  }
  while (!debug_stack.empty()) {
    auto f = debug_stack.back();
    debug_stack.pop_back();
    if (debug_visited[f]) continue;
    face = f;
    break;
  }
  if (face == -1) return;

  debug_visited[face] = true;

  debug_result.push_back(face);

  auto tag = mesh.tags[face];
  auto adj = mesh.adjacencies[face];
  printf("\nfrom %d: tag(%d %d %d) adj(%d %d %d)\n", face, tag[0], tag[1],
      tag[2], adj[0], adj[1], adj[2]);

  for (auto neighbor : mesh.adjacencies[face]) {
    if (neighbor < 0 || debug_visited[neighbor]) continue;
    auto tag = mesh.tags[neighbor];
    auto adj = mesh.adjacencies[neighbor];
    if (check(face, neighbor)) {
      debug_stack.push_back(neighbor);
      printf("ok   %d: tag(%d %d %d) adj(%d %d %d)\n", neighbor, tag[0], tag[1],
          tag[2], adj[0], adj[1], adj[2]);
    }
    printf("no   %d: tag(%d %d %d) adj(%d %d %d)\n", neighbor, tag[0], tag[1],
        tag[2], adj[0], adj[1], adj[2]);
  }
}

// // Polygon operations (from previous implementation)
// inline void polygon_and(const vector<cell_polygon>& cells,
//     vector<int>& cell_ids, const int polygon) {
//   for (auto i = 0; i < cells.size(); i++) {
//     auto& label = cells[i].embedding[polygon];
//     cell_ids[i] = cell_ids[i] && label;
//   }
// }

// inline void polygon_or(const vector<cell_polygon>& cells, vector<int>&
// cell_ids,
//     const int polygon) {
//   for (auto i = 0; i < cells.size(); i++) {
//     auto& label = cells[i].embedding[polygon];
//     cell_ids[i] = cell_ids[i] || label;
//   }
// }

// inline void polygon_not(const vector<cell_polygon>& cells,
//     vector<int>& cell_ids, const int polygon) {
//   for (auto i = 0; i < cells.size(); i++) {
//     auto& label = cells[i].embedding[polygon];
//     cell_ids[i] = !label;
//   }
// }

// inline void polygon_common(
//     const vector<cell_polygon>& cells, vector<int>& cell_ids, const int
//     num)
//     {
//   if (num < 1) return;

//   for (auto i = 0; i < cells.size(); i++) {
//     auto  sum   = 0;
//     auto& label = cells[i].embedding;
//     for (auto& l : label) sum += l;
//     cell_ids[i] = sum >= num;
//   }
//   return;
// }
