#pragma once

#include <yocto/yocto_mesh.h>
#include <yocto/yocto_shape.h>

using namespace yocto;

struct bezier_mesh {
  vector<vec3i>        triangles   = {};
  vector<vec3i>        adjacencies = {};
  vector<vec3f>        positions   = {};
  vector<vec3f>        normals     = {};
  dual_geodesic_solver dual_solver = {};
};

struct polygon {
  vector<mesh_point>    points = {};
  vector<geodesic_path> paths  = {};
};

inline bool is_updated(const polygon& polyg) {
  return (polyg.points.size() - 1) == polyg.paths.size();
}

inline bool is_closed(const polygon& polyg) {
  if (polyg.points.size() < 3) return false;
  return (polyg.points.front().face == polyg.points.back().face) &&
         (polyg.points.front().uv == polyg.points.back().uv);
}

inline bezier_mesh init_bezier_mesh(const generic_shape* shape) {
  auto mesh        = bezier_mesh{};
  mesh.triangles   = shape->triangles;
  mesh.normals     = shape->normals;
  mesh.positions   = shape->positions;
  mesh.adjacencies = face_adjacencies(mesh.triangles);
  mesh.dual_solver = make_dual_geodesic_solver(
      mesh.triangles, mesh.positions, mesh.adjacencies);
  return mesh;
}

inline geodesic_path compute_geodesic_path(
    const bezier_mesh& mesh, const mesh_point& start, const mesh_point& end) {
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