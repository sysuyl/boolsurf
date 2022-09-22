#include "boolsurf.h"

#include <cinolib/homotopy_basis.h>
#include <yocto/yocto_parallel.h>

#include <cassert>
#include <deque>

#include "ext/CDT/CDT/include/CDT.h"
namespace yocto {

constexpr auto adjacent_to_nothing = -2;

static bool_state* global_state = nullptr;

#define DEBUG_DATA 0
#if DEBUG_DATA
#define add_debug_triangle(face, triangle) debug_triangles()[face] = triangle
#else
#define add_debug_triangle(face, triangle) ;
#endif
#if DEBUG_DATA
#define add_debug_edge(face, edge) debug_edges()[face] = edge
#else
#define add_debug_edge(face, edge) ;
#endif
#if DEBUG_DATA
#define add_debug_node(face, node) debug_nodes()[face] = node
#else
#define add_debug_node(face, node) ;
#endif
#if DEBUG_DATA
#define add_debug_index(face, index) debug_indices()[face] = index
#else
#define add_debug_index(face, index) ;
#endif

static int scope_timer_indent = 0;
scope_timer::scope_timer(const string& msg) {
  if (scope_timer_indent == 0) printf("       \n");
  printf("[timer]");
  printf(" %.*s", scope_timer_indent, "|||||||||||||||||||||||||");
  // printf("%d", scope_timer_indent);
  printf("%s started\n", msg.c_str());
  scope_timer_indent += 1;
  message    = msg;
  start_time = get_time_();
}

scope_timer::~scope_timer() {
  scope_timer_indent -= 1;
  if (start_time < 0) return;
  auto elapsed = get_time_() - start_time;
  printf("[timer]");
  printf(" %.*s", scope_timer_indent, "|||||||||||||||||||||||||");
  // printf("%d", scope_timer_indent);
  printf("%s %s\n", message.c_str(), format_duration(elapsed).c_str());
}

// Build adjacencies between faces (sorted counter-clockwise)
static vector<vec3i> face_adjacencies_fast(const vector<vec3i>& triangles) {
  auto get_edge = [](const vec3i& triangle, int i) -> vec2i {
    auto x = triangle[i], y = triangle[i < 2 ? i + 1 : 0];
    return x < y ? vec2i{x, y} : vec2i{y, x};
  };

  auto adjacencies = vector<vec3i>{triangles.size(), vec3i{-1, -1, -1}};
  auto edge_map    = hash_map<vec2i, int>();
  edge_map.reserve((size_t)(triangles.size() * 1.5));
  for (int i = 0; i < triangles.size(); ++i) {
    for (int k = 0; k < 3; ++k) {
      auto edge = get_edge(triangles[i], k);
      auto it   = edge_map.find(edge);
      if (it == edge_map.end()) {
        edge_map[edge] = i;
      } else {
        auto neighbor     = it->second;
        adjacencies[i][k] = neighbor;
        for (int kk = 0; kk < 3; ++kk) {
          auto edge2 = get_edge(triangles[neighbor], kk);
          if (edge2 == edge) {
            adjacencies[neighbor][kk] = i;
            break;
          }
        }
      }
    }
  }
  return adjacencies;
}

void init_mesh(bool_mesh& mesh) {
  if (mesh.quads.size()) {
    mesh.triangles = quads_to_triangles(mesh.quads);
    mesh.quads.clear();
  }

  mesh.normals        = compute_normals(mesh);
  mesh.adjacencies    = face_adjacencies_fast(mesh.triangles);
  mesh.triangle_rings = vertex_to_triangles(
      mesh.triangles, mesh.positions, mesh.adjacencies);
  mesh.num_triangles = (int)mesh.triangles.size();
  mesh.num_positions = (int)mesh.positions.size();

  // Fit shape in [-1, +1]^3
  auto bbox = invalidb3f;
  for (auto& pos : mesh.positions) bbox = merge(bbox, pos);
  for (auto& pos : mesh.positions) pos = (pos - center(bbox)) / max(size(bbox));

  mesh.bbox     = bbox;
  mesh.bbox.min = (mesh.bbox.min - center(bbox)) / max(size(bbox));
  mesh.bbox.max = (mesh.bbox.max - center(bbox)) / max(size(bbox));
  mesh.bvh      = make_triangles_bvh(mesh.triangles, mesh.positions, {});

  mesh.dual_solver = make_dual_geodesic_solver(
      mesh.triangles, mesh.positions, mesh.adjacencies);

  mesh.graph = make_geodesic_solver(
      mesh.triangles, mesh.adjacencies, mesh.positions);

  mesh.triangles.reserve(mesh.triangles.size() * 2);
  mesh.adjacencies.reserve(mesh.adjacencies.size() * 2);
}

std::tuple<vector<int>, mesh_point, mesh_point> handle_short_strips(
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<int>& strip, const mesh_point& start, const mesh_point& end) {
  if (strip.size() == 1) {
    return {strip, start, end};
  } else if (strip.size() == 2) {
    auto [inside, b2f] = point_in_triangle(triangles, positions, start.face,
        eval_position(triangles, positions, end));
    if (inside) {
      auto new_end = mesh_point{start.face, b2f};
      return {{start.face}, start, new_end};
    }
    std::tie(inside, b2f) = point_in_triangle(triangles, positions, end.face,
        eval_position(triangles, positions, start));

    if (inside) {
      auto new_start = mesh_point{end.face, b2f};
      return {{end.face}, new_start, end};
    }

    return {strip, start, end};
  }
  return {{-1}, {}, {}};
}

vec3f flip_bary_to_adjacent_tri(const vector<vec3i>& adjacencies,
    const int tid0, const int tid1, const vec3f& bary) {
  if (tid0 == tid1) return bary;
  auto new_bary = zero3f;
  auto k1       = find_in_vec(adjacencies[tid1], tid0);
  auto k0       = find_in_vec(adjacencies[tid0], tid1);
  if (k1 == -1) {
    std::cout << "Error, faces are not adjacent" << std::endl;
    return zero3f;
  }
  new_bary[k1]           = bary[(k0 + 1) % 3];
  new_bary[(k1 + 1) % 3] = bary[k0];
  new_bary[(k1 + 2) % 3] = bary[(k0 + 2) % 3];

  return new_bary;
}

static vec3f get_bary(const vec2f& uv) {
  return vec3f{1 - uv.x - uv.y, uv.x, uv.y};
}

std::tuple<vector<int>, mesh_point, mesh_point> cleaned_strip(
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<int>& strip,
    const mesh_point& start, const mesh_point& end) {
  vector<int> cleaned = strip;

  auto start_entry = 0, end_entry = (int)strip.size() - 1;
  auto b3f           = zero3f;
  auto new_start     = start;
  auto new_end       = end;
  auto [is_vert, kv] = point_is_vert(end);
  auto [is_edge, ke] = point_is_edge(end);
  if (strip.size() <= 2)
    return handle_short_strips(triangles, positions, strip, start, end);
  // Erasing from the bottom
  if (is_vert) {
    auto vid      = triangles[end.face][kv];
    auto curr_tid = strip[end_entry - 1];
    kv            = find_in_vec(triangles[curr_tid], vid);
    while (kv != -1) {
      cleaned.pop_back();
      --end_entry;
      if (end_entry == 1) break;
      // see comment below
      auto curr_tid = strip[end_entry - 1];
      kv            = find_in_vec(triangles[curr_tid], vid);
    }

    kv = find_in_vec(triangles[cleaned.back()], vid);
    assert(kv != -1);
    b3f[kv] = 1;
    new_end = mesh_point{cleaned.back(), vec2f{b3f.y, b3f.z}};  // updating end
  } else if (is_edge) {
    if (end.face != strip.back()) {
      assert(adjacencies[end.face][ke] == strip.back());

      if (end.face == strip[end_entry - 1]) cleaned.pop_back();
    } else if (adjacencies[end.face][ke] == strip[end_entry - 1])
      cleaned.pop_back();

    b3f = flip_bary_to_adjacent_tri(
        adjacencies, end.face, cleaned.back(), get_bary(end.uv));

    new_end = mesh_point{cleaned.back(), vec2f{b3f.y, b3f.z}};  // updating end
  }
  std::tie(is_vert, kv) = point_is_vert(start);
  std::tie(is_edge, ke) = point_is_vert(start);

  if (is_vert) {
    auto vid      = triangles[start.face][kv];
    auto curr_tid = strip[start_entry + 1];
    kv            = find_in_vec(triangles[curr_tid], vid);
    while (kv != -1) {
      cleaned.erase(cleaned.begin());
      ++start_entry;
      if (start_entry > end_entry - 1) break;
      auto curr_tid = strip[start_entry + 1];
      kv            = find_in_vec(triangles[curr_tid], vid);
    }
    kv = find_in_vec(triangles[cleaned[0]], vid);
    assert(kv != -1);
    b3f       = zero3f;
    b3f[kv]   = 1;
    new_start = mesh_point{cleaned[0], vec2f{b3f.y, b3f.z}};  // udpdating start

  } else if (is_edge) {
    if (start.face != strip[0]) {
      assert(adjacencies[start.face][ke] == strip[0]);
      if (start.face == strip[1]) cleaned.erase(cleaned.begin());
    } else if (adjacencies[start.face][ke] == strip[1]) {
      cleaned.erase(cleaned.begin());
    }
    b3f = flip_bary_to_adjacent_tri(
        adjacencies, start.face, cleaned[0], get_bary(start.uv));
    new_start = {cleaned[0], vec2f{b3f.y, b3f.z}};  // updating start
  }
  return {cleaned, new_start, new_end};
}

void remove_loops_from_strip(vector<int>& strip) {
  auto faces      = unordered_map<int, int>{};
  faces[strip[0]] = 0;
  auto result     = vector<int>(strip.size());
  result[0]       = strip[0];
  auto index      = 1;
  for (auto i = 1; i < strip.size(); ++i) {
    if (faces.count(strip[i]) != 0) {
      // printf("fixing %d (%d)\n", i, strip[i]);
      auto t = faces[strip[i]];
      index  = t + 1;
      continue;
    }
    faces[strip[i]] = i;
    result[index++] = strip[i];
  }
  result.resize(index);
  strip = result;
}

void reset_mesh(bool_mesh& mesh) {
  mesh.triangles.resize(mesh.num_triangles);
  mesh.positions.resize(mesh.num_positions);
  mesh.adjacencies.resize(mesh.num_triangles);
  mesh.dual_solver.graph.resize(mesh.num_triangles);
  mesh.triangulated_faces.clear();

  auto get_triangle_center = [](const vector<vec3i>&  triangles,
                                 const vector<vec3f>& positions,
                                 int                  face) -> vec3f {
    vec3f pos[3] = {positions[triangles[face].x], positions[triangles[face].y],
        positions[triangles[face].z]};
    auto  l0     = length(pos[0] - pos[1]);
    auto  p0     = (pos[0] + pos[1]) / 2;
    auto  l1     = length(pos[1] - pos[2]);
    auto  p1     = (pos[1] + pos[2]) / 2;
    auto  l2     = length(pos[2] - pos[0]);
    auto  p2     = (pos[2] + pos[0]) / 2;
    return (l0 * p0 + l1 * p1 + l2 * p2) / (l0 + l1 + l2);
  };

  for (auto& [face, _] : mesh.triangulated_faces) {
    for (int k = 0; k < 3; k++) {
      auto neighbor = mesh.adjacencies[face][k];
      if (neighbor == -1) continue;
      auto kk = find_adjacent_triangle(
          mesh.triangles[neighbor], mesh.triangles[face]);

      // Fix adjacencies and dual_solver.
      mesh.adjacencies[neighbor][kk]              = face;
      mesh.dual_solver.graph[neighbor][kk].node   = face;
      mesh.dual_solver.graph[neighbor][kk].length = length(
          get_triangle_center(mesh.triangles, mesh.positions, neighbor) -
          get_triangle_center(mesh.triangles, mesh.positions, face));
    }
  }
}

cinolib::Trimesh<> get_cinomesh(const bool_mesh& mesh) {
  auto verts = vector<cinolib::vec3d>();
  verts.reserve(mesh.positions.size());
  for (auto& pos : mesh.positions)
    verts.push_back(cinolib::vec3d{pos.x, pos.y, pos.z});

  auto triangles = vector<vector<cinolib::uint>>();
  triangles.reserve(mesh.triangles.size());
  for (auto& tr : mesh.triangles) {
    auto triangle = vector<cinolib::uint>{(uint)tr.x, (uint)tr.y, (uint)tr.z};
    triangles.push_back(triangle);
  }

  auto cinomesh = cinolib::Trimesh<>(verts, triangles);
  return cinomesh;
}

void update_boolmesh(const cinolib::Trimesh<>& cinomesh, bool_mesh& mesh) {
  auto positions = vector<vec3f>();
  positions.reserve(cinomesh.num_verts());
  for (auto v = 0; v < cinomesh.num_verts(); v++) {
    auto pos = cinomesh.vert(v);
    positions.push_back({float(pos[0]), float(pos[1]), float(pos[2])});
  }

  auto triangles = vector<vec3i>();
  triangles.reserve(cinomesh.num_polys());
  for (auto p = 0; p < cinomesh.num_polys(); p++) {
    auto tri = cinomesh.poly_verts_id(p);
    triangles.push_back({(int)tri[0], (int)tri[1], (int)tri[2]});
  }

  mesh.positions = positions;
  mesh.triangles = triangles;
}

bool_homotopy_basis compute_homotopy_basis(bool_mesh& mesh, int root) {
  auto cinomesh = get_cinomesh(mesh);

  auto homotopy_basis_data           = cinolib::HomotopyBasisData{};
  homotopy_basis_data.detach_loops   = true;
  homotopy_basis_data.root           = root;
  homotopy_basis_data.split_strategy = cinolib::EDGE_SPLIT_STRATEGY;
  homotopy_basis(cinomesh, homotopy_basis_data);

  update_boolmesh(cinomesh, mesh);

  auto homo_basis = bool_homotopy_basis{};
  homo_basis.root = root;
  homo_basis.basis.reserve(homotopy_basis_data.loops.size());
  homo_basis.lengths.reserve(homotopy_basis_data.loops.size());

  for (auto& base : homotopy_basis_data.loops) {
    auto new_base = vector<int>(base.begin(), base.end());

    auto distance = 0.0f;
    for (auto e = 0; e < new_base.size(); e++) {
      auto start = new_base[e];
      auto end   = new_base[(e + 1) % new_base.size()];
      distance += length(mesh.positions[start] - mesh.positions[end]);
    }

    printf("Distance: %f\n", distance);
    printf("Size: %d\n", new_base.size());
    homo_basis.basis.push_back(new_base);
    homo_basis.lengths.push_back(distance);
  }

  return homo_basis;
}

vector<mesh_polygon> smooth_homotopy_basis(
    const bool_homotopy_basis& homotopy_basis, const bool_mesh& mesh,
    bool smooth_generators) {
  auto  num_basis    = mesh.homotopy_basis.basis.size();
  auto& basis        = mesh.homotopy_basis.basis;
  auto& root         = mesh.homotopy_basis.root;
  auto  smooth_basis = vector<mesh_polygon>();

  for (auto b = 0; b < basis.size(); b++) {
    auto& base  = basis[b];
    auto  strip = compute_strip_from_basis(
         base, mesh.triangle_rings, mesh.triangles, root);

    auto& first = strip.front();
    auto& last  = strip.back();

    auto& ftri = mesh.triangles[first];
    auto& ltri = mesh.triangles[last];

    auto start_point = mesh_point{first, get_uv_from_vertex(ftri, root)};
    auto end_point   = mesh_point{last, get_uv_from_vertex(ltri, root)};

    auto path = shortest_path(mesh.triangles, mesh.positions, mesh.adjacencies,
        start_point, end_point, strip);

    auto find_closest_to_vertex = [&]() {
      auto mid = (int)path.strip.size() / 2;
      return mid;
    };

    // Adjusting strip
    if (smooth_generators) {
      for (auto t = 0; t < 10; t++) {
        auto clean_result = cleaned_strip(mesh.triangles, mesh.positions,
            mesh.adjacencies, path.strip, start_point, end_point);

        path.strip = std::get<0>(clean_result);

        auto mid = find_closest_to_vertex();

        // Old first and last face of the strip
        auto first = path.strip.front();
        auto last  = path.strip.back();

        auto adj_first    = mesh.adjacencies[first];
        auto are_adjacent = find_in_vec(adj_first, last);

        // New first and last of the strip
        auto second_face = path.strip[mid + 1];
        auto first_face  = path.strip[mid];

        auto kk = find_adjacent_triangle(
            mesh.triangles[second_face], mesh.triangles[first_face]);
        auto [cc, dd] = get_triangle_uv_from_index(kk);
        auto mid_uv1  = lerp(cc, dd, 1 - 0.5f);
        auto mp1      = mesh_point{second_face, mid_uv1};

        // New last face of the strip
        auto k = find_adjacent_triangle(
            mesh.triangles[first_face], mesh.triangles[second_face]);
        auto [aa, bb] = get_triangle_uv_from_index(k);
        auto mid_uv   = lerp(aa, bb, 0.5f);
        auto mp2      = mesh_point{first_face, mid_uv};

        auto res = vector<int>{};
        res.reserve(path.strip.size() + 12);
        res.insert(res.end(), path.strip.begin() + (mid + 1), path.strip.end());

        // Add fan (only for root vertex)
        if (t == 0) {
          auto root_fan = mesh.triangle_rings[root];
          auto last_idx = find_idx(root_fan, last);

          while (true) {
            last_idx = (last_idx + 1) % root_fan.size();
            if (root_fan[last_idx] == first) break;
            res.push_back(root_fan[last_idx]);
          }
        }

        res.insert(res.end(), path.strip.begin(), path.strip.begin() + mid + 1);
        // check_triangle_strip(mesh.adjacencies, res);

        path = shortest_path(
            mesh.triangles, mesh.positions, mesh.adjacencies, mp1, mp2, res);
        start_point = mp1;
        end_point   = mp2;
      }
    }

    auto basis_shortest_segments = mesh_segments(
        mesh.triangles, path.strip, path.lerps, path.start, path.end);

    auto& smooth_base_polygon = smooth_basis.emplace_back();
    smooth_base_polygon.edges.push_back(basis_shortest_segments);
    smooth_base_polygon.length = basis_shortest_segments.size();
  }
  return smooth_basis;
}

void compute_homotopy_basis_borders(bool_mesh& mesh) {
  auto& homology_basis_borders = mesh.homotopy_basis_borders;
  auto& homology_basis         = mesh.homotopy_basis;

  for (auto b = 0; b < homology_basis.basis.size(); b++) {
    auto& base = homology_basis.basis[b];

    auto distance = 0.0;
    for (auto e = 0; e < base.size(); e++) {
      auto start = base[e];
      auto end   = base[(e + 1) % base.size()];

      distance += length(mesh.positions[end] - mesh.positions[start]);
      auto edge = vec2i{start, end};
      // auto rev_edge = vec2i{(int)end, (int)start};

      homology_basis_borders[edge] = pair<int, float>(b + 1, distance);
      // homology_basis_borders[rev_edge] = -homology_basis_borders[edge];
    }
  }
}

vector<int> sort_homotopy_basis_around_vertex(
    const bool_mesh& mesh, const bool_homotopy_basis& basis) {
  auto rotate = [](const vec3i& v, int k) {
    if (mod3(k) == 0)
      return v;
    else if (mod3(k) == 1)
      return vec3i{v.y, v.z, v.x};
    else
      return vec3i{v.z, v.x, v.y};
  };

  auto ordered_basis = vector<int>();
  for (auto& tri : mesh.triangle_rings[basis.root]) {
    auto triangle  = mesh.triangles[tri];
    auto root_idx  = find_in_vec(triangle, basis.root);
    auto rtriangle = rotate(triangle, root_idx);

    for (auto k = 0; k < 3; k++) {
      auto edge = get_mesh_edge_from_index(rtriangle, k);
      if (contains(mesh.homotopy_basis_borders, edge)) {
        auto& info = mesh.homotopy_basis_borders.at(edge);

        // Outgoing from root (negative) - Ingoing in root (positive)
        auto base_sign = (edge.x == basis.root) ? 1 : -1;
        ordered_basis.push_back(base_sign * info.first);
        printf("Triangle: %d - %d - %d\n", tri, k, base_sign * info.first);
      }
    }
  }

  return ordered_basis;
}

vector<int> compute_polygonal_schema(const vector<int>& basis) {
  auto inverse_mapping = hash_map<int, int>();
  for (auto i = 0; i < basis.size(); i++) {
    inverse_mapping[basis[i]] = i;
  }

  auto& start       = basis.front();
  auto  current_idx = inverse_mapping[-start] + 1;
  auto  current     = basis[current_idx % basis.size()];

  auto polygonal_schema = vector<int>();
  polygonal_schema.reserve(basis.size());
  polygonal_schema.push_back(start);

  while (current != start) {
    polygonal_schema.push_back(current);
    auto next_idx = inverse_mapping[-current] + 1;
    current       = basis[next_idx % basis.size()];
  }

  return polygonal_schema;
}

vector<pair<int, float>> compute_polygon_basis_intersections(
    const mesh_polygon& polygon, bool_mesh& mesh) {
  auto basis_intersections = vector<pair<int, float>>();
  if (polygon.length == 0) return basis_intersections;

  // Computing intersections between polygon and homotopy_basis
  for (auto& edge : polygon.edges) {
    for (auto& seg : edge) {
      auto [k, _] = get_edge_lerp_from_uv(seg.end);
      if (k == -1) continue;

      auto edge     = get_mesh_edge_from_index(mesh.triangles[seg.face], k);
      auto rev_edge = vec2i{edge.y, edge.x};

      if (contains(mesh.homotopy_basis_borders, edge)) {
        auto& info = mesh.homotopy_basis_borders.at(edge);
        basis_intersections.push_back(info);
      }
      if (contains(mesh.homotopy_basis_borders, rev_edge)) {
        auto& info = mesh.homotopy_basis_borders.at(rev_edge);
        info.first *= -1;

        basis_intersections.push_back(info);
      }
    }
  }
  return basis_intersections;
}

vector<int> compute_polygon_word(const vector<pair<int, float>>& isecs,
    const vector<int>& polygonal_schema) {
  // Polygon code by from polygon representation of polygonal schema
  auto polygon_word         = vector<int>();
  auto inv_polygonal_schema = hash_map<int, int>{};
  for (auto id = 0; id < polygonal_schema.size(); id++) {
    inv_polygonal_schema[polygonal_schema[id]] = id;
  }

  auto transform_point = [&](const pair<int, float>& point1) {
    auto& [base1, distance1] = point1;
    auto base_side           = inv_polygonal_schema[base1];

    auto next_side = (base_side + 1) % polygonal_schema.size();
    auto next_base = polygonal_schema[next_side];
    auto next_distance =
        (sign(base1) != sign(next_base)) ? distance1 : abs(distance1 - 1.0f);

    return pair<int, float>{next_base, next_distance};
  };

  for (auto c = 0; c < isecs.size(); c++) {
    auto point1       = isecs[c];
    auto trans_point1 = transform_point(point1);

    auto next_intersection = (c + 1) % isecs.size();
    auto point2            = isecs[next_intersection];
    point2.first           = -point2.first;

    printf("%d (%f) [%d (%f)] -> %d (%f)\n", point1.first, point1.second,
        trans_point1.first, trans_point1.second, point2.first, point2.second);

    if (point1 == point2) continue;
    if (trans_point1 == point2) continue;

    auto& [base_id1, dist1] = trans_point1;
    auto& [base_id2, dist2] = point2;

    auto base_side1 = inv_polygonal_schema[base_id1];
    auto base_side2 = inv_polygonal_schema[base_id2];

    // Reading the correct portion of polygonal schema
    auto id_dist = abs((int)base_side2 - (int)base_side1);
    // printf("From: %d to: %d\n", base_side1, base_side2);
    auto current_id = base_side1;
    while (current_id != base_side2) {
      auto current_side = polygonal_schema[current_id];
      // printf("\t%d ", current_side);
      polygon_word.push_back(current_side);
      current_id = (current_id + 1) % polygonal_schema.size();
    }
    // printf("\n");
  }
  return polygon_word;
}

mesh_polygon vectorize_generator_loop(
    bool_state& state, const mesh_polygon& generator_loop, int orientation) {
  auto generator_curve      = mesh_polygon{};
  generator_curve.is_closed = true;

  auto num_curves = 5;
  auto interval   = (int)(generator_loop.length / (num_curves * 3));

  for (auto e = 0; e < generator_loop.edges.size(); e++) {
    for (auto s = 0; s < generator_loop.edges[e].size(); s += interval) {
      auto& segment       = generator_loop.edges[e][s];
      auto  uvx           = lerp(segment.start.x, segment.end.x, 0.5f);
      auto  uvy           = lerp(segment.start.y, segment.end.y, 0.5f);
      auto  control_point = mesh_point{segment.face, vec2f{uvx, uvy}};

      generator_curve.points.push_back((int)state.points.size());
      state.points.push_back(control_point);
    }
  }

  if (orientation > 0)
    reverse(generator_curve.points.begin(), generator_curve.points.end());

  return generator_curve;
}

vector<int> compute_strip_from_basis_old(const vector<int>& base,
    const vector<vector<int>>& triangle_rings, const vector<vec3i>& triangles,
    int root) {
  auto strip      = vector<int>();
  auto root_ring  = triangle_rings[root];
  auto last_edge  = vec2i{base.back(), base.front()};
  auto first_edge = vec2i{base[0], base[1]};

  auto first_idx = -1;
  for (auto f = 0; f < root_ring.size(); f++) {
    auto triangle_verts = triangles[root_ring[f]];
    if (!edge_in_triangle(triangle_verts, last_edge)) continue;
    first_idx = f;
    break;
  }

  first_idx = (first_idx + 1) % root_ring.size();
  for (auto f = 0; f < root_ring.size(); f++) {
    auto it         = (first_idx + f) % root_ring.size();
    auto face       = root_ring[it];
    auto face_verts = triangles[face];

    if (edge_in_triangle(face_verts, first_edge)) break;
    strip.push_back(face);
  }

  for (auto e = 0; e < base.size(); e++) {
    auto next_edge = vec2i{base[e % base.size()], base[(e + 1) % base.size()]};
    auto vertex_ring = triangle_rings[base[e]];

    auto face_idx = find_idx(vertex_ring, strip.back()) + 1;
    for (auto f = 0; f < vertex_ring.size(); f++) {
      auto it         = (face_idx + f) % vertex_ring.size();
      auto face       = vertex_ring[it];
      auto face_verts = triangles[face];

      if (edge_in_triangle(face_verts, next_edge)) break;
      strip.push_back(face);
    }
  }

  if (strip.back() == strip.front()) strip.pop_back();

  while (!contains(root_ring, strip.back()))
    rotate(strip.begin(), strip.begin() + 1, strip.end());

  return strip;
}

vector<int> compute_strip_from_basis(const vector<int>& base,
    const vector<vector<int>>& triangle_rings, const vector<vec3i>& triangles,
    int root) {
  auto strip = vector<int>();

  auto edge_to_face = hash_map<vec2i, int>();
  for (auto e = 0; e < base.size(); e++) {
    auto edge        = vec2i{base[(e + 1) % base.size()], base[e]};
    auto vertex_ring = triangle_rings[base[e]];

    for (auto& face : vertex_ring) {
      auto face_verts = triangles[face];

      if (edge_in_triangle(face_verts, edge)) {
        edge_to_face[edge] = face;
      }
    }
  }

  auto last_triangle = 0;
  for (auto e = 1; e < base.size(); e++) {
    auto prev_idx = (e != 0) ? e - 1 : base.size() - 1;
    auto curr_idx = e % base.size();
    auto next_idx = (e + 1) % base.size();

    auto prev_v = base[prev_idx];
    auto curr_v = base[curr_idx];
    auto next_v = base[next_idx];

    auto prev_edge = vec2i{curr_v, prev_v};
    auto next_edge = vec2i{next_v, curr_v};

    auto prev_face = edge_to_face.at(prev_edge);
    auto next_face = edge_to_face.at(next_edge);

    auto triangle_ring = triangle_rings[curr_v];

    // printf("Triangle ring of %d:\n", curr_v);
    // for (auto& tri : triangle_ring) {
    //   printf("%d ", tri);
    // }
    // printf("\n");

    auto start = find_idx(triangle_ring, prev_face);
    auto end   = find_idx(triangle_ring, next_face);

    // printf("Start idx: %d (%d) (%d) - end idx: %d (%d) (%d)  \n", start,
    //     triangle_ring[start], prev_face, end, triangle_ring[end], next_face);

    for (auto f = 0; f < triangle_ring.size(); f++) {
      auto idx = (f + start) % triangle_ring.size();
      if (idx == end) {
        last_triangle = triangle_ring[idx];
        break;
      } else {
        strip.push_back(triangle_ring[idx]);
      }
      // printf("%d -> %d\n", idx, triangle_ring[idx]);
    }
    // printf("\n");
    // printf("\n");
  }

  // printf("Last triangle: %d\n", last_triangle);
  strip.push_back(last_triangle);

  return strip;
}

geodesic_path compute_geodesic_path(
    const bool_mesh& mesh, const mesh_point& start, const mesh_point& end) {
  auto path = geodesic_path{};
  if (start.face == end.face) {
    path.start = start;
    path.end   = end;
    path.strip = {start.face};
    return path;
  }

  auto strip = compute_strip(
      mesh.dual_solver, mesh.triangles, mesh.positions, end, start);
  path = shortest_path(
      mesh.triangles, mesh.positions, mesh.adjacencies, start, end, strip);
  return path;
}

mesh_point eval_geodesic_path(
    const bool_mesh& mesh, const geodesic_path& path, float t) {
  return eval_path_point(
      path, mesh.triangles, mesh.positions, mesh.adjacencies, t);
}

vector<mesh_segment> mesh_segments(const vector<vec3i>& triangles,
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

void recompute_polygon_segments(const bool_mesh& mesh, const bool_state& state,
    mesh_polygon& polygon, int index) {
  // Remove this to automatically close the curves
  if (polygon.points.size() == 0) return;

  if (index > 0) {
    auto& last_segment = polygon.edges.back();
    polygon.length -= last_segment.size();
    polygon.edges.pop_back();
  } else {
    polygon.length = 0;
    polygon.edges.clear();
  }

  auto faces = hash_set<int>();

  auto num_points = polygon.is_closed ? polygon.points.size()
                                      : polygon.points.size() - 1;

  for (int i = index; i < num_points; i++) {
    auto start = polygon.points[i];

    auto end = polygon.is_closed
                   ? polygon.points[(i + 1) % polygon.points.size()]
                   : polygon.points[(i + 1)];

    faces.insert(state.points[start].face);

    auto path = compute_geodesic_path(
        mesh, state.points[start], state.points[end]);
    auto threshold = 0.001f;
    for (auto& l : path.lerps) {
      l = yocto::clamp(l, 0 + threshold, 1 - threshold);
    }
    auto segments = mesh_segments(
        mesh.triangles, path.strip, path.lerps, path.start, path.end);

    polygon.edges.push_back(segments);
    polygon.length += segments.size();
  }

  if (index == polygon.points.size() - 1) {
    auto start = polygon.points[index];
    auto end   = polygon.points[0];

    auto path = compute_geodesic_path(
        mesh, state.points[start], state.points[end]);
    auto threshold = 0.001f;
    for (auto& l : path.lerps) {
      l = yocto::clamp(l, 0 + threshold, 1 - threshold);
    }
    auto segments = mesh_segments(
        mesh.triangles, path.strip, path.lerps, path.start, path.end);

    polygon.edges.push_back(segments);
    polygon.length += segments.size();
  }

  polygon.is_contained_in_single_face = (faces.size() == 1);
}

inline int num_segments(const hashgrid_polyline& polyline) {
  if (polyline.is_closed) return (int)polyline.points.size();
  return (int)polyline.points.size() - 1;
}

inline pair<vec2f, vec2f> get_segment(
    const hashgrid_polyline& polyline, int i) {
  if (polyline.is_closed) {
    return {
        polyline.points[i], polyline.points[(i + 1) % polyline.points.size()]};
  } else {
    return {polyline.points[i], polyline.points[i + 1]};
  }
}

inline vec2i get_segment_vertices(const hashgrid_polyline& polyline, int i) {
  if (polyline.is_closed) {
    return {polyline.vertices[i],
        polyline.vertices[(i + 1) % polyline.vertices.size()]};
  } else {
    return {polyline.vertices[i], polyline.vertices[i + 1]};
  }
}

// struct mesh_hashgrid : hash_map<int, vector<hashgrid_polyline>> {};

inline int add_vertex(bool_mesh& mesh, mesh_hashgrid& hashgrid,
    const mesh_point& point, int polyline_id, int vertex = -1) {
  float eps             = 0.00001;
  auto  update_polyline = [&](int v) {
    if (polyline_id < 0) return;
    auto& polyline = hashgrid[point.face][polyline_id];
    polyline.vertices.push_back(v);
    polyline.points.push_back(point.uv);
  };

  {  // Maybe collapse with original mesh vertices.
    auto uv = point.uv;
    auto tr = mesh.triangles[point.face];
    if (uv.x < eps && uv.y < eps) {
      update_polyline(tr.x);
      return tr.x;
    }
    if (uv.x > 1 - eps && uv.y < eps) {
      update_polyline(tr.y);
      return tr.y;
    }
    if (uv.y > 1 - eps && uv.x < eps) {
      update_polyline(tr.z);
      return tr.z;
    }
  }

  {  // Maybe collapse with already added vertices.
    auto& polylines = hashgrid[point.face];
    for (auto& polyline : polylines) {
      for (int i = 0; i < polyline.vertices.size(); i++) {
        if (length(point.uv - polyline.points[i]) < eps) {
          update_polyline(polyline.vertices[i]);
          return polyline.vertices[i];
        }
      }
    }
  }

  // No collapse. Add new vertex to mesh.
  if (vertex == -1) {
    vertex   = (int)mesh.positions.size();
    auto pos = eval_position(mesh.triangles, mesh.positions, point);
    mesh.positions.push_back(pos);
  }

  update_polyline(vertex);
  return vertex;
}

static mesh_hashgrid compute_hashgrid(bool_mesh& mesh,
    const vector<shape>& shapes, hash_map<int, int>& control_points) {
  _PROFILE();
  // La hashgrid associa ad ogni faccia una lista di polilinee.
  // Ogni polilinea è definita da una sequenza punti in coordinate
  // baricentriche, ognuno di essi assiociato al corrispondente vertice della
  // mesh.
  auto hashgrid = mesh_hashgrid{};

  for (auto shape_id = 0; shape_id < shapes.size(); shape_id++) {
    auto& polygons = shapes[shape_id].polygons;
    for (auto polygon_id = 0; polygon_id < polygons.size(); polygon_id++) {
      auto& polygon = polygons[polygon_id];
      if (polygon.length == 0) continue;
      if (polygon.edges.empty()) continue;

      // Open polygons are handled in update_hashgrid()
      if (!polygon.is_closed) {
        continue;
      }

      // La polilinea della prima faccia del poligono viene processata alla
      // fine (perché si trova tra il primo e l'ultimo edge)
      int  first_face   = polygon.edges[0][0].face;
      int  first_vertex = -1;
      auto indices      = vec2i{-1, -1};  // edge_id, segment_id

      int last_face   = -1;
      int last_vertex = -1;

      for (auto e = 0; e < polygon.edges.size(); e++) {
        auto& edge = polygon.edges[e];

        for (auto s = 0; s < edge.size(); s++) {
          auto& segment = edge[s];

          // Iniziamo a riempire l'hashgrid a partire da quando troviamo una
          // faccia diversa da quella iniziale del poligono (il primo tratto
          // verrà aggiunto a posteriori per evitare inconsistenza)
          if (segment.face == first_face && indices == vec2i{-1, -1}) continue;
          if (indices == vec2i{-1, -1}) indices = {e, s};

          auto& entry = hashgrid[segment.face];
          auto  ids   = vec2i{e, s};
          ids.y       = (s + 1) % edge.size();
          ids.x       = ids.y > s ? e : (e + 1) % polygon.edges.size();

          // Se la faccia del segmento che stiamo processando è diversa
          // dall'ultima salvata allora creiamo una nuova polilinea,
          // altrimenti accodiamo le nuove informazioni.
          if (segment.face != last_face) {
            auto  polyline_id = (int)entry.size();
            auto& polyline    = entry.emplace_back();
            // polyline.polygon  = polygon_id;
            polyline.shape_id   = shape_id;
            polyline.polygon_id = polygon_id;

            last_vertex = add_vertex(mesh, hashgrid,
                {segment.face, segment.start}, polyline_id, last_vertex);
            if (first_vertex == -1) first_vertex = last_vertex;

            last_vertex = add_vertex(
                mesh, hashgrid, {segment.face, segment.end}, polyline_id);

          } else {
            auto  polyline_id = (int)entry.size() - 1;
            auto& polyline    = entry.back();
            assert(segment.end != polyline.points.back());

            last_vertex = add_vertex(
                mesh, hashgrid, {segment.face, segment.end}, polyline_id);
          }

          last_face = segment.face;
        }

        if (last_vertex != -1)
          control_points[last_vertex] =
              polygon.points[(e + 1) % polygon.edges.size()];
      }

      if (indices == vec2i{-1, -1}) {
        auto& entry       = hashgrid[first_face];
        auto  polyline_id = (int)entry.size();
        auto& polyline    = entry.emplace_back();
        // polyline.polygon   = polygon_id;
        polyline.shape_id   = shape_id;
        polyline.polygon_id = polygon_id;

        polyline.is_closed = true;

        for (auto e = 0; e < polygon.edges.size(); e++) {
          auto& edge = polygon.edges[e];
          for (int s = 0; s < edge.size(); s++) {
            auto& segment = edge[s];

            last_vertex = add_vertex(
                mesh, hashgrid, {segment.face, segment.start}, polyline_id);
          }

          if (last_vertex != -1)
            control_points[last_vertex] =
                polygon.points[(e + 1) % polygon.edges.size()];
        }
      };

      // Ripetiamo parte del ciclo (fino a indices) perché il primo tratto di
      // polilinea non è stato inserito nell'hashgrid
      auto vertex = -1;
      for (auto e = 0; e <= indices.x; e++) {
        auto end_idx = (e < indices.x) ? polygon.edges[e].size() : indices.y;
        for (auto s = 0; s < end_idx; s++) {
          auto ids = vec2i{e, s};
          ids.y    = (s + 1) % polygon.edges[e].size();
          ids.x    = ids.y > s ? e : e + 1;

          auto& segment     = polygon.edges[e][s];
          auto& entry       = hashgrid[segment.face];
          auto  polyline_id = (int)entry.size() - 1;

          if (e == indices.x && s == indices.y - 1) vertex = first_vertex;

          // auto& polyline    = entry.back();
          // assert(segment.face == last_face);
          last_vertex = add_vertex(
              mesh, hashgrid, {segment.face, segment.end}, polyline_id, vertex);
        }

        if (e > 0 && last_vertex != -1)
          control_points[last_vertex] =
              polygon.points[(e + 1) % polygon.edges.size()];
      }
    }
  }
  return hashgrid;
}

// Qui vengono gestiti solo i poligoni aperti (per comodità in una funzione
// separata)
static mesh_hashgrid compute_open_shapes_hashgrid(bool_mesh& mesh,
    const vector<shape>& shapes, hash_map<int, int>& control_points) {
  _PROFILE();

  auto hashgrid = mesh_hashgrid();
  for (auto shape_id = 0; shape_id < shapes.size(); shape_id++) {
    auto& polygons = shapes[shape_id].polygons;

    auto shape_vertices = vector<vector<vector<int>>>();
    for (auto polygon_id = 0; polygon_id < polygons.size(); polygon_id++) {
      auto& polygon = polygons[polygon_id];
      if (polygon.length == 0) continue;
      if (polygon.edges.empty()) continue;
      if (polygon.is_closed) continue;

      // La polilinea della prima faccia del poligono viene processata alla
      // fine (perché si trova tra il primo e l'ultimo edge)

      int last_face   = -1;
      int last_vertex = -1;

      auto polygon_vertices = vector<vector<int>>();
      for (auto e = 0; e < polygon.edges.size(); e++) {
        auto& edge = polygon.edges[e];

        for (auto s = 0; s < edge.size(); s++) {
          auto& segment = edge[s];

          auto& entry = hashgrid[segment.face];

          // Se la faccia del segmento che stiamo processando è diversa
          // dall'ultima salvata allora creiamo una nuova polilinea,
          // altrimenti accodiamo le nuove informazioni.
          if (segment.face != last_face) {
            auto  polyline_id = (int)entry.size();
            auto& polyline    = entry.emplace_back();
            // polyline.polygon  = polygon_id;
            polyline.shape_id   = shape_id;
            polyline.polygon_id = polygon_id;

            last_vertex = add_vertex(mesh, hashgrid,
                {segment.face, segment.start}, polyline_id, last_vertex);

            last_vertex = add_vertex(
                mesh, hashgrid, {segment.face, segment.end}, polyline_id);

          } else {
            auto  polyline_id = (int)entry.size() - 1;
            auto& polyline    = entry.back();
            assert(segment.end != polyline.points.back());

            last_vertex = add_vertex(
                mesh, hashgrid, {segment.face, segment.end}, polyline_id);
          }
          last_face = segment.face;
        }

        if (last_vertex != -1)
          control_points[last_vertex] =
              polygon.points[(e + 1) % polygon.edges.size()];
      }
    }
  }

  return hashgrid;
}

[[maybe_unused]] static hash_map<int, int> compute_control_points(
    vector<mesh_polygon>&             polygons,
    const vector<vector<vector<int>>> vertices) {
  auto control_points = hash_map<int, int>();
  for (auto p = 0; p < vertices.size(); p++) {
    for (auto e = 0; e < vertices[p].size(); e++) {
      auto control_point_idx            = vertices[p][e][0];
      auto mesh_point_idx               = polygons[p].points[e];
      control_points[control_point_idx] = mesh_point_idx;
    }
  }
  return control_points;
}

void save_tree_png(const bool_state& state, string filename,
    const string& extra, bool color_shapes);

vector<mesh_cell> make_mesh_cells(vector<int>& cell_tags,
    const vector<vec3i>& adjacencies, const vector<bool>& border_tags) {
  auto result = vector<mesh_cell>{};
  cell_tags   = vector<int>(adjacencies.size(), -1);

  // consume task stack
  auto starts = vector<int>{(int)adjacencies.size() - 1};

  while (starts.size()) {
    auto start = starts.back();
    starts.pop_back();
    if (cell_tags[start] >= 0) continue;

    // pop element from task stack
    auto first_face = start;

    // static int c = 0;
    // // save_tree_png(*global_state,
    // // "data/tests/flood_fill_" + to_string(c) + ".png", "", false);
    // c += 1;

    auto  cell_id = (int)result.size();
    auto& cell    = result.emplace_back();
    cell.faces.reserve(adjacencies.size());
    auto face_stack = vector<int>{first_face};

    while (!face_stack.empty()) {
      auto face = face_stack.back();
      face_stack.pop_back();

      if (cell_tags[face] >= 0) continue;
      cell_tags[face] = cell_id;

      cell.faces.push_back(face);

      for (int k = 0; k < 3; k++) {
        auto neighbor = adjacencies[face][k];
        if (neighbor < 0) continue;

        auto neighbor_cell = cell_tags[neighbor];
        if (neighbor_cell >= 0) continue;
        if (border_tags[3 * face + k]) {
          starts.push_back(neighbor);
        } else {
          face_stack.push_back(neighbor);
        }
      }
    }  // end of while
    cell.faces.shrink_to_fit();
  }  // end of while

  return result;
}

vector<mesh_cell> make_cell_graph(bool_mesh& mesh) {
  _PROFILE();
  // Iniziamo dall'ultima faccia che sicuramente non e' stata distrutta.
  auto cells = make_mesh_cells(
      mesh.face_tags, mesh.adjacencies, mesh.borders.tags);

  {
    _PROFILE_SCOPE("tag_cell_edges");
    for (auto& [ids, faces] : mesh.polygon_borders) {
      auto& shape_id = ids.x;
      for (auto& [inner_face, outer_face] : faces) {
        if (inner_face < 0 || outer_face < 0) continue;
        auto a = mesh.face_tags[inner_face];
        auto b = mesh.face_tags[outer_face];
        cells[a].adjacency.insert({b, -shape_id});
        cells[b].adjacency.insert({a, +shape_id});
      }
    }
  }

  return cells;
}

static vector<int> find_roots(const vector<mesh_cell>& cells) {
  // Trova le celle non hanno archi entranti con segno di poligono positivo.
  auto adjacency = vector<int>(cells.size(), 0);
  for (auto& cell : cells) {
    for (auto& [adj, p] : cell.adjacency) {
      if (p > 0) adjacency[adj] += 1;
    }
  }

  auto result = vector<int>{};
  for (int i = 0; i < adjacency.size(); i++) {
    if (adjacency[i] == 0) result.push_back(i);
  }
  return result;
}

static vector<mesh_cell> compute_shape_macrograph(
    const vector<mesh_cell>& cells, int shape_id) {
  auto components      = vector<vector<int>>();
  auto shape_component = hash_map<int, int>();
  auto visited         = vector<bool>(cells.size(), false);

  for (auto c = 0; c < cells.size(); c++) {
    if (visited[c]) continue;

    auto  component_id = (int)components.size();
    auto& component    = components.emplace_back();

    auto stack = vector<int>();
    stack.push_back(c);

    while (!stack.empty()) {
      auto cell_id = stack.back();
      stack.pop_back();

      if (visited[cell_id]) continue;
      visited[cell_id] = true;
      component.push_back(cell_id);
      shape_component[cell_id] = component_id;

      auto& cell = cells[cell_id];
      for (auto& [neighbor, shape] : cell.adjacency) {
        if (yocto::abs(shape) == shape_id) continue;
        if (visited[neighbor]) continue;
        stack.push_back(neighbor);
      }
    }
  }

  auto macrograph = vector<mesh_cell>((int)components.size());
  for (auto c = 0; c < (int)components.size(); c++) {
    for (auto id : components[c]) {
      for (auto [neighbor, shape] : cells[id].adjacency) {
        if (yocto::abs(shape) != shape_id) continue;

        auto neighbor_component = shape_component.at(neighbor);
        macrograph[c].adjacency.insert({neighbor_component, shape});
      }
    }
  }

  return macrograph;
}

static void compute_cycles(const vector<mesh_cell>& cells, int node,
    vec2i parent, vector<int> visited, vector<vec2i> parents,
    vector<vector<vec2i>>& cycles) {
  // Se il nodo il considerazione è già stato completamente visitato allora
  // terminiamo la visita
  if (visited[node] == 2) return;

  // Se il nodo in considerazione non è stato completamente visitato e lo
  // stiamo rivisitando ora allora abbiamo trovato un ciclo
  if (visited[node] == 1) {
    auto  cycle   = vector<vec2i>();
    auto& current = parent;
    cycle.push_back(current);

    // Risalgo l'albero della visita fino a che non trovo lo stesso nodo e
    // salvo il ciclo individuato
    while (current.x != node) {
      auto prev = parents[current.x];

      // (marzia) check: è vero che ho un ciclo corretto se il verso
      // (entrante/uscente) è lo stesso per tutti gli archi?
      // if (sign(prev.y) != sign(current.y)) return;
      current = prev;
      cycle.push_back(current);
    }

    cycles.push_back(cycle);
    return;
  }

  // Settiamo il padre del nodo attuale e iniziamo ad esplorare i suoi vicini
  parents[node] = parent;
  visited[node] = 1;

  for (auto& [neighbor, polygon] : cells[node].adjacency) {
    // Se stiamo percorrendo lo stesso arco ma al contrario allora continuo,
    // altrimenti esploriamo il vicino
    // if (polygon > 0) continue;
    // if (neighbor == parent.x && polygon == -parent.y) continue;
    compute_cycles(cells, neighbor, {node, polygon}, visited, parents, cycles);
  }

  // Settiamo il nodo attuale come completamente visitato
  visited[node] = 2;
}

inline vector<vector<vec2i>> compute_graph_cycles(
    const vector<mesh_cell>& cells) {
  // _PROFILE();
  auto visited        = vector<int>(cells.size(), 0);
  auto parents        = vector<vec2i>(cells.size(), {0, 0});
  auto cycles         = vector<vector<vec2i>>();
  auto start_node     = 0;
  auto invalid_parent = vec2i{-1, -1};
  compute_cycles(cells, start_node, invalid_parent, visited, parents, cycles);
  return cycles;
}

hash_set<int> compute_invalid_shapes(
    const vector<mesh_cell>& cells, int num_shapes) {
  _PROFILE();
  auto invalid_shapes = hash_set<int>();
  for (auto s = 1; s < num_shapes; s++) {
    auto shape_graph = compute_shape_macrograph(cells, s);
    auto cycles      = compute_graph_cycles(shape_graph);

    for (auto cycle : cycles)
      if (cycle.size() % 2 == 1) invalid_shapes.insert(s);
  }
  return invalid_shapes;
}

inline vector<vector<int>> compute_components(
    const bool_state& state, const shape& bool_shape) {
  // Calcoliamo le componenti tra le celle presenti in una bool_shape
  // (per calcolarne i bordi in maniera più semplice)
  auto cells   = vector<int>(bool_shape.cells.begin(), bool_shape.cells.end());
  auto visited = hash_map<int, bool>();
  for (auto cell : cells) visited[cell] = false;

  auto components = vector<vector<int>>();

  for (auto cell : cells) {
    if (visited[cell]) continue;

    auto& component = components.emplace_back();

    auto stack = vector<int>();
    stack.push_back(cell);

    while (!stack.empty()) {
      auto cell_idx = stack.back();
      stack.pop_back();

      if (visited[cell_idx]) continue;
      visited[cell_idx] = true;
      component.push_back(cell_idx);

      auto& cell = state.cells[cell_idx];
      for (auto& [neighbor, shape] : cell.adjacency) {
        if (find_idx(cells, neighbor) == -1) continue;
        if (visited[neighbor]) continue;
        stack.push_back(neighbor);
      }
    }
  }
  return components;
}

static vector<vector<int>> propagate_cell_labels(bool_state& state) {
  _PROFILE();
  // Inizializziamo le label delle celle a 0.
  auto  num_shapes     = (int)state.bool_shapes.size();
  auto& cells          = state.cells;
  auto& labels         = state.labels;
  auto& invalid_shapes = state.invalid_shapes;

  labels = vector<vector<int>>(
      cells.size(), vector<int>(num_shapes, null_label));

  // TODO: Initialization of visit (is it ok?)
  auto new_start = vector<int>();
  new_start.push_back((int)cells.size() - 1);

  for (auto cell_id : new_start) {
    auto& cell      = cells[cell_id];
    labels[cell_id] = vector<int>(num_shapes, 0);
    for (auto& [neighbor_id, shape_id] : cell.adjacency) {
      if (shape_id > 0) continue;
      auto ushape_id             = yocto::abs(shape_id);
      labels[cell_id][ushape_id] = 1;
    }
  }

  auto queue   = std::deque<int>(new_start.begin(), new_start.end());
  auto visited = vector<bool>(cells.size(), false);
  for (auto& s : new_start) visited[s] = true;

  while (!queue.empty()) {
    // print("queue", queue);
    auto cell_id = queue.front();
    queue.pop_front();
    // static int c = 0;
    // save_tree_png(
    //     *global_state, "data/tests/" + to_string(c) + ".png", "", false);
    // c += 1;

    auto& cell = cells[cell_id];
    for (auto& [neighbor_id, shape_id] : cell.adjacency) {
      auto ushape_id = yocto::abs(shape_id);

      auto& neighbor_labels = labels[neighbor_id];
      auto  updated_labels  = labels[cell_id];
      updated_labels[ushape_id] += yocto::sign(shape_id);

      auto updated = false;
      for (int s = 0; s < neighbor_labels.size(); s++) {
        if (neighbor_labels[s] == null_label) {
          neighbor_labels[s] = updated_labels[s];
          updated            = true;
        }
      }

      for (int s = 0; s < neighbor_labels.size(); s++) {
        if (updated_labels[s] % 2 != neighbor_labels[s] % 2) {
          invalid_shapes.insert(s);
        }
      }

      if (updated) {
        if (!contains(queue, neighbor_id)) {
          queue.push_back(neighbor_id);
        }
      }

      visited[neighbor_id] = true;
    }
  }
  return labels;
}

static void add_polygon_intersection_points(bool_state& state,
    hash_map<int, vector<hashgrid_polyline>>& hashgrid, bool_mesh& mesh) {
  _PROFILE();
  // Calcoliamo sia le intersezioni che le self-intersections, aggiungendo i
  // vertici nuovi alla mesh.

  for (auto& [face, polylines] : hashgrid) {
    // Check for polyline self interesctions
    for (auto p0 = 0; p0 < polylines.size(); p0++) {
      auto& poly = polylines[p0];

      for (int s0 = 0; s0 < num_segments(poly) - 1; s0++) {
        auto [start0, end0] = get_segment(poly, s0);
        int num_added       = 0;

        for (int s1 = s0 + 1; s1 < num_segments(poly); s1++) {
          // Skip adjacent segments.
          if (poly.is_closed) {
            if (yocto::abs(s0 - s1) % num_segments(poly) <= 1) continue;
          } else {
            if (yocto::abs(s0 - s1) <= 1) continue;
          }

          auto [start1, end1] = get_segment(poly, s1);

          auto l = intersect_segments(start0, end0, start1, end1);
          if (l.x <= 0.0f || l.x >= 1.0f || l.y <= 0.0f || l.y >= 1.0f) {
            continue;
          }

          auto uv                      = lerp(start1, end1, l.y);
          auto point                   = mesh_point{face, uv};
          auto vertex                  = add_vertex(mesh, hashgrid, point, -1);
          state.control_points[vertex] = (int)state.points.size();
          state.isecs_generators[vertex] = {poly.shape_id, poly.shape_id};

          state.points.push_back(point);
          // printf("self-intersection: polygon %d, vertex %d\n",
          // poly.polygon,
          //     vertex);

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
        auto& poly0 = polylines[p0];
        auto& poly1 = polylines[p1];
        for (int s0 = 0; s0 < num_segments(poly0); s0++) {
          auto [start0, end0] = get_segment(poly0, s0);
          int num_added       = 0;

          for (int s1 = 0; s1 < num_segments(poly1); s1++) {
            auto [start1, end1] = get_segment(poly1, s1);

            auto l = intersect_segments(start0, end0, start1, end1);
            if (l.x <= 0.0f || l.x >= 1.0f || l.y <= 0.0f || l.y >= 1.0f) {
              continue;
            }

            auto uv     = lerp(start1, end1, l.y);
            auto point  = mesh_point{face, uv};
            auto vertex = add_vertex(mesh, hashgrid, point, -1);
            state.control_points[vertex]   = (int)state.points.size();
            state.isecs_generators[vertex] = {poly0.shape_id, poly1.shape_id};

            state.points.push_back(point);

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

triangulation_info triangulation_constraints(const bool_mesh& mesh, int face,
    const vector<hashgrid_polyline>& polylines) {
  auto info    = triangulation_info{};
  info.face    = face;
  info.nodes   = vector<vec2f>{{0, 0}, {1, 0}, {0, 1}};
  info.indices = vector<int>(
      &mesh.triangles[face][0], &mesh.triangles[face][3]);

  for (auto& polyline : polylines) {
    // Per ogni segmento della polyline, aggiorniamo triangulation_info,
    // aggiungendo nodi, indici, e edge constraints.
    for (auto i = 0; i < num_segments(polyline); i++) {
      vec2f uvs[2];
      std::tie(uvs[0], uvs[1]) = get_segment(polyline, i);
      auto vertices            = get_segment_vertices(polyline, i);
      assert(vertices[0] < mesh.positions.size());
      assert(vertices[1] < mesh.positions.size());

      // TODO(giacomo): Set to -1 or 'invalid'.
      auto local_vertices = vec2i{-7, -8};
      for (int k = 0; k < 2; k++) {
        local_vertices[k] = find_idx(info.indices, vertices[k]);
        if (local_vertices[k] == -1) {
          info.indices.push_back(vertices[k]);
          info.nodes.push_back(uvs[k]);
          local_vertices[k] = (int)info.indices.size() - 1;
        }
      }

      // Aggiungiamo l'edge ai constraints della triangolazione.
      info.edges.push_back({local_vertices[0], local_vertices[1]});

      // Extra: Se i nodi sono su un lato k != -1 di un triangolo allora li
      // salviamo nella edgemap.
      for (int j = 0; j < 2; j++) {
        auto [k, lerp] = get_edge_lerp_from_uv(uvs[j]);
        if (k != -1) {
          info.edgemap[k].push_back({local_vertices[j], lerp});
        }
      }
    }
  }
  return info;
}

void add_boundary_edge_constraints(
    array<vector<pair<int, float>>, 3>& edgemap, vector<vec2i>& edges) {
  // Aggiungiamo gli edge di vincolo per i lati del triangolo.
  for (int k = 0; k < 3; k++) {
    auto  edge        = get_triangle_edge_from_index(k);
    auto& edge_points = edgemap[k];

    // Se sul lato non ci sono punti aggiuntivi, allora lo aggiungiamo ai
    // vincoli cosi' come e'.
    if (edge_points.empty()) {
      edges.push_back(edge);
      continue;
    }

    // Ordiniamo i punti che giacciono sul lato per lerp crescente.
    if (edge_points.size() > 1) {
      sort(edge_points.begin(), edge_points.end(),
          [](auto& a, auto& b) { return a.second < b.second; });
    }

    // Creiamo primo e ultimo vincolo di questo lato.
    edges.push_back({edge.x, edge_points[0].first});
    edges.push_back({edge.y, edge_points.back().first});

    // Creiamo i rimanenti vincoli contenuti nel lato.
    for (auto i = 0; i < edge_points.size() - 1; i++) {
      edges.push_back({edge_points[i].first, edge_points[i + 1].first});
    }
  }
}

static pair<vector<vec3i>, vector<vec3i>> single_split_triangulation(
    const vector<vec2f>& nodes, const vec2i& edge) {
  // Calcoliamo la triangolazione con un singolo segmento all'interno del
  // triangolo.
  auto start_edge = get_edge_from_uv(nodes[edge.x]);
  auto end_edge   = get_edge_from_uv(nodes[edge.y]);

  auto triangles   = vector<vec3i>(3);
  auto adjacencies = vector<vec3i>(3);

  // image: libs/boolsurf/notes/sinlge-split-adjacency.jpg
  if (edge.x < 3) {
    // Se il segmento ha come inizio un punto in un lato e come fine il
    // vertice del triangolo opposto
    triangles[0] = {edge.x, end_edge.x, edge.y};
    triangles[1] = {edge.x, edge.y, end_edge.y};
    triangles.resize(2);

    adjacencies[0] = {adjacent_to_nothing, adjacent_to_nothing, 1};
    adjacencies[1] = {0, adjacent_to_nothing, adjacent_to_nothing};
    adjacencies.resize(2);

  } else if (edge.y < 3) {
    // Se il segmento ha come inizio un vertice di un triangolo e come fine un
    // punto punto nel lato opposto
    triangles[0] = {edge.y, start_edge.x, edge.x};
    triangles[1] = {edge.y, edge.x, start_edge.y};
    triangles.resize(2);

    adjacencies[0] = {adjacent_to_nothing, adjacent_to_nothing, 1};
    adjacencies[1] = {0, adjacent_to_nothing, adjacent_to_nothing};
    adjacencies.resize(2);

  } else {
    // Se il segmento ha inizio e fine su due lati del triangolo
    auto [x, y] = start_edge;
    if (start_edge.y == end_edge.x) {
      auto z       = end_edge.y;
      triangles[0] = {x, edge.x, z};
      triangles[1] = {edge.x, edge.y, z};
      triangles[2] = {edge.x, y, edge.y};

      adjacencies[0] = {adjacent_to_nothing, 1, adjacent_to_nothing};
      adjacencies[1] = {2, adjacent_to_nothing, 0};
      adjacencies[2] = {adjacent_to_nothing, adjacent_to_nothing, 1};

    } else if (start_edge.x == end_edge.y) {
      auto z       = end_edge.x;
      triangles[0] = {x, edge.x, edge.y};
      triangles[1] = {edge.x, z, edge.y};
      triangles[2] = {edge.x, y, z};

      adjacencies[0] = {adjacent_to_nothing, 1, adjacent_to_nothing};
      adjacencies[1] = {2, adjacent_to_nothing, 0};
      adjacencies[2] = {adjacent_to_nothing, adjacent_to_nothing, 1};

    } else {
      assert(0);
    }
  }

  return {triangles, adjacencies};
}

// Constrained Delaunay Triangulation
static pair<vector<vec3i>, vector<vec3i>> constrained_triangulation(
    vector<vec2f>& nodes, const vector<vec2i>& edges, int face) {
  // Questo purtroppo serve.
  for (auto& n : nodes) n *= 1e9;

  auto  result      = pair<vector<vec3i>, vector<vec3i>>{};
  auto& triangles   = result.first;
  auto& adjacencies = result.second;

  // (marzia): qui usiamo float, ma si possono usare anche i double
  using Float = float;
  auto cdt = CDT::Triangulation<Float>(CDT::FindingClosestPoint::ClosestRandom);
  cdt.insertVertices(
      nodes.begin(), nodes.end(),
      [](const vec2f& point) -> Float { return point.x; },
      [](const vec2f& point) -> Float { return point.y; });
  cdt.insertEdges(
      edges.begin(), edges.end(),
      [](const vec2i& edge) -> int { return edge.x; },
      [](const vec2i& edge) -> int { return edge.y; });

  cdt.eraseOuterTriangles();
  adjacencies.reserve(cdt.triangles.size());

  triangles.reserve(cdt.triangles.size());

  for (auto& tri : cdt.triangles) {
    auto verts = vec3i{
        (int)tri.vertices[0], (int)tri.vertices[1], (int)tri.vertices[2]};

    auto adjacency = vec3i{};
    for (auto k = 0; k < 3; k++) {
      auto neigh = tri.neighbors[k];
      if (neigh == CDT::noNeighbor)
        adjacency[k] = adjacent_to_nothing;
      else
        adjacency[k] = (int)neigh;
    }

#if DEBUG_DATA
    // TODO: serve? (marzia): Forse no!
    auto& a           = nodes[verts.x];
    auto& b           = nodes[verts.y];
    auto& c           = nodes[verts.z];
    auto  orientation = cross(b - a, c - b);
    if (fabs(orientation) < 0.00001) {
      global_state->failed = true;
      printf("[%s]: Collinear in face : %d\n", __FUNCTION__, face);
      return {};
    }
#endif

    triangles.push_back(verts);
    adjacencies.push_back(adjacency);
  }
  return result;
}

static void update_face_adjacencies(bool_mesh& mesh) {
  _PROFILE();
  // Aggiorniamo le adiacenze per i triangoli che sono stati processati
  auto border_edgemap = hash_map<vec2i, int>{};
  border_edgemap.reserve(mesh.triangulated_faces.size() * 6);

  // Per ogni triangolo processato elaboro tutti i suoi sottotriangoli
  for (auto& [face, faces] : mesh.triangulated_faces) {
    // Converto il triangolo in triplette di vertici globali.
    auto triangles = vector<vec3i>(faces.size());
    for (int i = 0; i < faces.size(); i++) {
      triangles[i] = mesh.triangles[faces[i].id];
    }

    for (int i = 0; i < faces.size(); i++) {
      // Guardo se nell'adiacenza ci sono dei triangoli mancanti
      // (segnati con adjacent_to_nothing per non confonderli i -1 già
      // presenti nell'adiacenza della mesh originale).
      auto& adj = mesh.adjacencies[faces[i].id];
      for (int k = 0; k < 3; k++) {
        if (adj[k] != adjacent_to_nothing) continue;

        // Prendo l'edge di bordo corrispondente ad un adjacent_to_nothing
        auto edge = get_mesh_edge_from_index(triangles[i], k);

        // Se è un arco della mesh originale lo processo subito
        if (edge.x < mesh.num_positions && edge.y < mesh.num_positions) {
          // Cerco il triangolo adiacente al triangolo originale su quel lato
          for (int kk = 0; kk < 3; kk++) {
            auto edge0 = get_mesh_edge_from_index(mesh.triangles[face], kk);
            if (make_edge_key(edge) == make_edge_key(edge0)) {
              // Aggiorno direttamente l'adiacenza nel nuovo triangolo e del
              // vicino
              auto neighbor = mesh.adjacencies[face][kk];
              if (neighbor == -1) continue;

              mesh.adjacencies[faces[i].id][k] = neighbor;

              auto it = find_in_vec(mesh.adjacencies[neighbor], face);
              mesh.adjacencies[neighbor][it] = faces[i].id;
            }
          }
          continue;
        }

        // Se non è un arco della mesh originale
        auto edge_key = make_edge_key(edge);
        auto it       = border_edgemap.find(edge_key);

        // Se non l'ho mai incontrato salvo in una mappa l'edge e la
        // faccia corrispondente. Se l'ho già incontrato ricostruisco
        // l'adiacenza tra il triangolo corrente e il neighbor già trovato.
        if (it == border_edgemap.end()) {
          // border_edgemap.insert(it, {edge_key, faces[i]});
          border_edgemap[edge_key] = faces[i].id;
        } else {
          auto neighbor                    = it->second;
          mesh.adjacencies[faces[i].id][k] = neighbor;
          for (int kk = 0; kk < 3; ++kk) {
            auto edge2 = get_mesh_edge_from_index(mesh.triangles[neighbor], kk);
            edge2      = make_edge_key(edge2);
            if (edge2 == edge_key) {
              mesh.adjacencies[neighbor][kk] = faces[i].id;
              break;
            }
          }
        }
      }
    }
  }
}

inline bool check_tags(
    const bool_mesh& mesh, const vector<vec3i>& border_tags) {
  for (int i = 0; i < mesh.triangles.size(); i++) {
    if (mesh.triangulated_faces.find(i) != mesh.triangulated_faces.end()) {
      continue;
    }
    auto face = i;
    auto tr   = mesh.triangles[face];
    if (tr == vec3i{0, 0, 0}) continue;
    for (int k = 0; k < 3; k++) {
      auto neighbor = mesh.adjacencies[face][k];
      if (neighbor < 0) continue;
      // auto n0 = mesh.adjacencies[face];
      // auto n1 = mesh.adjacencies[neighbor];
      auto kk = find_in_vec(mesh.adjacencies[neighbor], face);
      assert(kk != -1);

      auto tags0 = border_tags[face];
      auto tags1 = border_tags[neighbor];
      auto tag0  = tags0[k];
      auto tag1  = tags1[kk];
      assert(tag0 == -tag1);
    }
  }
  return true;
}

bool check_polygon_validity(bool_mesh& mesh, int shape_id, int polygon_id) {
  auto& polygon_face_borders = mesh.polygon_borders[{shape_id, polygon_id}];
  auto  polygon_border_tags  = vector<bool>(3 * mesh.triangles.size(), false);
  compute_polygon_border_tags(mesh, polygon_face_borders, polygon_border_tags);

  auto polygon_face_tags = vector<int>(
      vector<int>(mesh.adjacencies.size(), -1));
  auto polygon_cells = make_mesh_cells(
      polygon_face_tags, mesh.adjacencies, polygon_border_tags);

  for (auto& [inner_face, outer_face] : polygon_face_borders) {
    auto& inner_face_tag = polygon_face_tags[inner_face];
    auto& outer_face_tag = polygon_face_tags[outer_face];

    if (inner_face_tag == outer_face_tag) {
      return false;
    }
  }
  return true;
}

vector<mesh_point> compute_parallel_loop(
    bool_mesh& mesh, const mesh_polygon& polygon) {
  auto parallel_points = vector<mesh_point>();
  parallel_points.reserve(polygon.points.size());

  for (auto e = 0; e < polygon.edges.size() - 2; e++) {
    if (!polygon.edges[e].size()) continue;
    auto segment     = polygon.edges[e][0];
    auto start_point = mesh_point{segment.face, segment.start};
    auto end_point   = mesh_point{segment.face, segment.end};

    // Computing the tangent (?)
    auto start_tr = triangle_coordinates(
        mesh.triangles, mesh.positions, end_point);
    auto tangent = interpolate_triangle(
        start_tr[0], start_tr[1], start_tr[2], segment.start);
    tangent     = normalize(tangent);
    auto normal = vec2f{tangent.y, -tangent.x};

    auto path     = straightest_path(mesh, start_point, normal, 0.03f);
    path.end.uv.x = clamp(path.end.uv.x, 0.0f, 1.0f);
    path.end.uv.y = clamp(path.end.uv.y, 0.0f, 1.0f);
    parallel_points.push_back(path.end);
  }

  // Parallel loop is traced in the opposite direction
  std::reverse(parallel_points.begin(), parallel_points.end());
  return parallel_points;
}

template <typename F>
inline void parallel_for_batch(int num_threads, size_t size, F&& f) {
  auto threads    = vector<std::thread>(num_threads);
  auto batch_size = size / num_threads;
  auto batch_job  = [&](size_t k) {
    auto from = k * batch_size;
    auto to   = std::min(from + batch_size, size);
    for (auto i = from; i < to; i++) f(i);
  };

  for (auto k = 0; k < num_threads; k++) {
    threads[k] = std::thread(batch_job, k);
  }
  for (auto k = 0; k < num_threads; k++) {
    threads[k].join();
  }
}

static void triangulate(bool_mesh& mesh, const mesh_hashgrid& hashgrid) {
  _PROFILE();
  // auto mesh_triangles_size = atomic<size_t>{mesh.triangles.size()};
  auto mesh_mutex = std::mutex{};
  auto i          = 0;
  auto faces      = vector<int>(hashgrid.size());
  for (auto& [key, _] : hashgrid) {
    faces[i++] = key;
  }
  mesh.triangulated_faces.reserve(hashgrid.size());

  // for (auto& [face, polylines] : hashgrid) {
  auto f = [&](size_t index) {
    auto  face      = faces[index];
    auto& polylines = hashgrid.at(face);

    // Calcola le info per la triangolazione, i.e. (nodi ed edge
    // constraints).
    auto info = triangulation_constraints(mesh, face, polylines);

    //    add_debug_node(face, info.nodes);
    //    add_debug_index(face, info.indices);

    // Se la faccia contiene solo segmenti corrispondenti ad edge del
    // triangolo stesso, non serve nessuna triangolazione.
    if (info.nodes.size() == 3) {
      auto nodes = std::array<vec2f, 3>{vec2f{0, 0}, vec2f{1, 0}, vec2f{0, 1}};
      mesh.triangulated_faces[face] = {facet{nodes, face}};
      return;
    }

    auto triangulated_faces = vector<facet>{};
    auto triangles          = vector<vec3i>();
    auto adjacency          = vector<vec3i>();

    // Se il triangolo ha al suo interno un solo segmento allora chiamiamo
    // la funzione di triangolazione più semplice, altrimenti chiamiamo CDT
    if (info.edges.size() == 1) {
      tie(triangles, adjacency) = single_split_triangulation(
          info.nodes, info.edges[0]);
    } else {
      add_boundary_edge_constraints(info.edgemap, info.edges);
      tie(triangles, adjacency) = constrained_triangulation(
          info.nodes, info.edges, face);
    }

    // Converti triangli locali in globali.
    for (int i = 0; i < triangles.size(); i++) {
      auto& tr = triangles[i];
      auto  n  = std::array<vec2f, 3>{
          info.nodes[tr[0]], info.nodes[tr[1]], info.nodes[tr[2]]};
      triangulated_faces.push_back({n, -1});
      tr = {info.indices[tr.x], info.indices[tr.y], info.indices[tr.z]};
    }

    vector<vec4i> polygon_faces;

    // Border map: from edge (expressed in mesh vertices to shape_id)
    auto border_map = hash_map<vec2i, vec2i>{};
    for (auto& polyline : polylines) {
      auto shape_id   = polyline.shape_id;
      auto polygon_id = polyline.polygon_id;

      for (auto i = 0; i < num_segments(polyline); i++) {
        auto edge        = get_segment_vertices(polyline, i);
        border_map[edge] = vec2i{shape_id, polygon_id};
      }
    }

    for (int i = 0; i < triangles.size(); i++) {
      for (int k = 0; k < 3; k++) {
        auto edge = get_mesh_edge_from_index(triangles[i], k);
        if (auto it = border_map.find(edge); it != border_map.end()) {
          auto ids = it->second;
          polygon_faces.push_back({ids.x, ids.y, i, adjacency[i][k]});
        }
      }
    }

    add_debug_edge(face, info.edges);
    add_debug_triangle(face, triangles);

    // TODO(giacomo): Pericoloso: se resize() innesca una riallocazione, il
    // codice dopo l'unlock che sta eseguendo su un altro thread puo' fare
    // casino.
    auto mesh_triangles_old_size = 0;
    {
      auto lock               = std::lock_guard{mesh_mutex};
      mesh_triangles_old_size = (int)mesh.triangles.size();
      mesh.triangles.resize(mesh_triangles_old_size + triangles.size());
      mesh.adjacencies.resize(mesh_triangles_old_size + triangles.size());
      for (int i = 0; i < triangles.size(); i++) {
        triangulated_faces[i].id = mesh_triangles_old_size + i;
      }
      mesh.triangulated_faces[face] = triangulated_faces;
      for (auto& pf : polygon_faces) {
        if (pf.z >= 0) pf.z += mesh_triangles_old_size;
        if (pf.w >= 0) pf.w += mesh_triangles_old_size;
        mesh.polygon_borders[{pf.x, pf.y}].push_back({pf.z, pf.w});
      }
    }

    for (int i = 0; i < triangles.size(); i++) {
      mesh.triangles[mesh_triangles_old_size + i] = triangles[i];
    }
    for (int i = 0; i < triangles.size(); i++) {
      for (auto& x : adjacency[i])
        if (x != adjacent_to_nothing) x += mesh_triangles_old_size;
      mesh.adjacencies[mesh_triangles_old_size + i] = adjacency[i];
    }
  };

// parallel_for_batch(8, hashgrid.size(), f);
#if DEBUG_DATA
  for (int i = 0; i < hashgrid.size(); i++) f(i);
#else
  parallel_for(hashgrid.size(), f);
#endif
}

void compute_polygon_border_tags(bool_mesh& mesh,
    const vector<vec2i>& polygon_borders, vector<bool>& border_tags) {
  for (auto& [inner_face, outer_face] : polygon_borders) {
    if (inner_face < 0 || outer_face < 0) continue;
    auto k = find_in_vec(mesh.adjacencies[inner_face], outer_face);
    assert(k != -1);
    border_tags[3 * inner_face + k] = true;
    auto kk = find_in_vec(mesh.adjacencies[outer_face], inner_face);
    assert(kk != -1);
    border_tags[3 * outer_face + kk] = true;
  }
}

void compute_border_tags(bool_mesh& mesh, bool_state& state) {
  _PROFILE();
  mesh.borders.tags = vector<bool>(3 * mesh.triangles.size(), false);
  for (auto& [ids, faces] : mesh.polygon_borders) {
    compute_polygon_border_tags(mesh, faces, mesh.borders.tags);
  }
}

void slice_mesh(bool_mesh& mesh, bool_state& state) {
  _PROFILE();
  auto& shapes = state.bool_shapes;

  // for (auto& shape : shapes) {
  //   for (auto& polygon : shape.polygons) {
  //     printf("Points: %d\n", polygon.points.size());
  //     printf("Length: %d\n", polygon.length);
  //   }
  //   printf("\n");
  // }

  // Calcoliamo i vertici nuovi della mesh
  // auto vertices             = add_vertices(mesh, polygons);
  state.num_original_points = (int)state.points.size();

  // Calcoliamo hashgrid e intersezioni tra poligoni,
  // aggiungendo ulteriori vertici nuovi alla mesh
  auto hashgrid = compute_hashgrid(mesh, shapes, state.control_points);

  auto open_hashgrid = compute_open_shapes_hashgrid(
      mesh, shapes, state.control_points);

  // Todo (marzia): Qui c'è ripetizione di dati
  auto total_hashgrid = mesh_hashgrid();
  for (auto& [key, polylines] : hashgrid) {
    for (auto& polyline : polylines) {
      total_hashgrid[key].push_back(polyline);
    }
  }

  for (auto& [key, polylines] : open_hashgrid) {
    for (auto& polyline : polylines) {
      total_hashgrid[key].push_back(polyline);
    }
  }

  add_polygon_intersection_points(state, total_hashgrid, mesh);
  add_polygon_intersection_points(state, hashgrid, mesh);

  state.hashgrid = total_hashgrid;

  // Triangolazione e aggiornamento dell'adiacenza
  triangulate(mesh, hashgrid);
  update_face_adjacencies(mesh);

  // Calcola i border_tags per le facce triangolata.
  compute_border_tags(mesh, state);
}

void compute_cell_labels(bool_state& state, bool non_zero) {
  _PROFILE();
  propagate_cell_labels(state);
}

bool compute_cells(bool_mesh& mesh, bool_state& state, bool non_zero) {
  // Triangola mesh in modo da embeddare tutti i poligoni come mesh-edges.
  _PROFILE();
  global_state = &state;
  slice_mesh(mesh, state);

  if (global_state->failed) return false;

  // Trova celle e loro adiacenza via flood-fill.
  state.cells = make_cell_graph(mesh);

  // Calcola i label delle celle con una visita sulla loro adiacenza.
  compute_cell_labels(state, non_zero);
  return true;
}

void compute_shapes(bool_state& state) {
  // Calcoliamo le informazioni sulla shape, come le celle che ne fanno parte
  auto& shapes  = state.bool_shapes;
  auto& sorting = state.shapes_sorting;
  sorting.resize(state.bool_shapes.size());

  // Assign a polygon and a color to each shape.
  for (auto s = 0; s < state.bool_shapes.size(); s++) {
    shapes[s].color = get_color(s);
    sorting[s]      = s;
  }

  // Distribute cells to shapes.
  // La prima shape è relativa alla cella ambiente, che è root per
  // definizione
  //  shapes[0].cells = hash_set<int>(
  //      state.ambient_cells.begin(), state.ambient_cells.end());
  //  shapes[0].is_root = false;

  for (auto c = 0; c < state.cells.size(); c++) {
    auto count = 0;
    for (auto p = 1; p < state.labels[c].size(); p++) {
      if (state.labels[c][p] > 0) {
        shapes[p].cells.insert(c);
        count += 1;
      }
    }

    if (count == 0) shapes[0].cells.insert(c);
  }
}

void compute_generator_polygons(
    const bool_state& state, int shape_idx, hash_set<int>& result) {
  // Calcoliamo ricorsivamente i poligoni iniziali coinvolti nelle operazioni
  // che hanno generato la shape corrente
  auto& bool_shape = state.bool_shapes[shape_idx];

  // Se la shape non ha generatori allora corrisponde ad una shape di un
  // poligono
  if (bool_shape.generators.empty()) {
    result.insert(shape_idx);
    return;
  }

  // Calcolo i generatori per le shape che hanno generato quella corrente
  for (auto& generator : bool_shape.generators) {
    if (generator == -1) continue;
    compute_generator_polygons(state, generator, result);
  }
}

void compute_shape_borders(const bool_mesh& mesh, bool_state& state) {
  // Calcoliamo tutti i bordi di una shape
  for (auto s = 1; s < state.bool_shapes.size(); s++) {
    auto& bool_shape = state.bool_shapes[s];

    // Calcoliamo il bordo solo per le shape root dell'albero csg
    if (!bool_shape.is_root) continue;

    // Calcoliamo i poligoni iniziali coinvolti nelle operazioni che hanno
    // generato la root (ci serve successivamente per salvare nel bordo
    // solamente i punti corretti)
    auto generator_polygons = hash_set<int>();
    compute_generator_polygons(state, s, generator_polygons);
    if (!generator_polygons.size()) continue;

    auto components = compute_components(state, bool_shape);

    for (auto& component : components) {
      // Step 1: Calcoliamo gli edges che stanno sul bordo
      auto edges = hash_set<vec2i>();

      for (auto c : component) {
        auto& cell = state.cells[c];
        // Per ogni cella che compone la shape calcolo il bordo a partire
        // dalle facce che ne fanno parte

        for (auto face : cell.faces) {
          // Se è una faccia interna allora non costituirà il bordo
          for (auto k = 0; k < 3; k++) {
            if (mesh.borders.tags[3 * face + k] == false) continue;
          }

          // Per ogni lato del triangolo considero solamente quelli che sono
          // di bordo (tag != 0)
          auto& tri = mesh.triangles[face];
          for (auto k = 0; k < 3; k++) {
            auto tag = mesh.borders.tags[3 * face + k];
            if (tag == false) continue;
            auto edge     = get_mesh_edge_from_index(tri, k);
            auto rev_edge = vec2i{edge.y, edge.x};

            // Se 'edge' è già stato incontrato allora esso è un bordo tra due
            // celle che fanno parte dela stessa shape, quindi lo elimino dal
            // set.
            auto it = edges.find(rev_edge);
            if (it == edges.end())
              edges.insert(edge);
            else
              edges.erase(it);
          }
        }
      }

      // Step 2: Riordiniamo i bordi
      // Per ogni vertice salviamo il proprio successivo
      auto next_vert = hash_map<int, int>();
      for (auto& edge : edges) next_vert[edge.x] = edge.y;

      for (auto& [key, value] : next_vert) {
        // Se il valore è -1 abbiamo già processato il punto
        if (value == -1) continue;

        // Aggiungiamo un nuovo bordo
        auto border_points = vector<int>();

        auto current = key;

        while (true) {
          auto next = next_vert.at(current);
          if (next == -1) break;

          next_vert.at(current) = -1;

          // Se il vertice corrente è un punto di controllo lo aggiungo al
          // bordo
          if (contains(state.control_points, current)) {
            // Se è un punto di intersezione controlliamo che i poligoni che
            // lo hanno generato siano entrambi compresi nei poligoni che
            // hanno generato anche la shape.
            if (contains(state.isecs_generators, current)) {
              auto& isec_generators = state.isecs_generators.at(current);

              if (contains(generator_polygons, isec_generators.x) &&

                  contains(generator_polygons, isec_generators.y))
                border_points.push_back(current);
            } else
              border_points.push_back(current);
          }

          // Se un bordo è stato chiuso correttamente lo inseriamo tra i bordi
          // della shape
          if (next == key) {
            bool_shape.border_points.push_back(border_points);
            break;
          } else
            current = next;
        }
      }
    }
  }
}

bool_state compute_border_polygons(const bool_state& state) {
  auto new_state   = bool_state{};
  new_state.points = state.points;

  for (auto& bool_shape : state.bool_shapes) {
    if (!bool_shape.is_root) continue;
    auto& test_shape = new_state.bool_shapes.emplace_back();
    for (auto& border : bool_shape.border_points) {
      auto& polygon = test_shape.polygons.emplace_back();
      for (auto v : border) {
        auto id = state.control_points.at(v);
        polygon.points.push_back(id);
      }
    }
  }
  return new_state;
}

void compute_bool_operation(bool_state& state, const bool_operation& op) {
  auto& a = state.bool_shapes[op.shape_a];
  auto& b = state.bool_shapes[op.shape_b];

  // Convertiamo il vettore di interi in bool per semplificare le operazioni
  auto aa = vector<bool>(state.cells.size(), false);
  for (auto& c : a.cells) aa[c] = true;

  auto bb = vector<bool>(state.cells.size(), false);
  for (auto& c : b.cells) bb[c] = true;

  if (op.type == bool_operation::Type::op_union) {
    for (auto i = 0; i < aa.size(); i++) aa[i] = aa[i] || bb[i];
  } else if (op.type == bool_operation::Type::op_intersection) {
    for (auto i = 0; i < aa.size(); i++) aa[i] = aa[i] && bb[i];
  } else if (op.type == bool_operation::Type::op_difference) {
    for (auto i = 0; i < aa.size(); i++) aa[i] = aa[i] && !bb[i];
  } else if (op.type == bool_operation::Type::op_symmetrical_difference) {
    for (auto i = 0; i < aa.size(); i++) aa[i] = aa[i] != bb[i];
  }

  // Le shape 'a' e 'b' sono state usate nell'operazione,
  // quindi non sono root del csg tree
  a.is_root = false;
  b.is_root = false;

  // Creiamo una nuova shape risultato, settando come generatori le shape 'a'
  // e 'b' e riconvertendo il vettore di bool a interi
  auto  shape_id = state.bool_shapes.size();
  auto& c        = state.bool_shapes.emplace_back();
  c.generators   = {op.shape_a, op.shape_b};
  c.color        = state.bool_shapes[op.shape_a].color;
  auto sorting   = find_idx(state.shapes_sorting, op.shape_a);

  insert(state.shapes_sorting, sorting, (int)shape_id);

  for (auto i = 0; i < aa.size(); i++)
    if (aa[i]) c.cells.insert(i);
}

void compute_bool_operations(
    bool_state& state, const vector<bool_operation>& ops) {
  _PROFILE();
  for (auto& op : ops) {
    compute_bool_operation(state, op);
  }
}

mesh_point intersect_mesh(const bool_mesh& mesh, const shape_bvh& bvh,
    const scene_camera& camera, const vec2f& uv) {
  auto ray = camera_ray(
      camera.frame, camera.lens, camera.aspect, camera.film, uv);
  auto isec = intersect_triangles_bvh(bvh, mesh.triangles, mesh.positions, ray);
  return {isec.element, isec.uv};
}

vec3f get_cell_color(const bool_state& state, int cell_id, bool color_shapes) {
  if (state.bool_shapes.empty() && state.labels.empty()) return {1, 1, 1};
  if (color_shapes) {
    for (int s = (int)state.shapes_sorting.size() - 1; s >= 0; s--) {
      auto& bool_shape = state.bool_shapes[state.shapes_sorting[s]];
      if (bool_shape.cells.count(cell_id) && bool_shape.is_root) {
        return bool_shape.color;
      }
    }
    return {1, 1, 1};
  } else {
    auto color = vec3f{0, 0, 0};
    int  count = 0;

    for (int p = 1; p < state.labels[cell_id].size(); p++) {
      auto label = state.labels[cell_id][p];
      if (label > 0) {
        color += get_color(p);
        count += 1;
      }
    }

    if (count > 0) {
      color /= count;
    } else {
      color = {0.9, 0.9, 0.9};
    }

    return color;
  }
}

hash_map<int, vector<vec3i>>& debug_triangles() {
  static hash_map<int, vector<vec3i>> result = {};
  return result;
}

hash_map<int, vector<vec2i>>& debug_edges() {
  static hash_map<int, vector<vec2i>> result = {};
  return result;
}

hash_map<int, vector<vec2f>>& debug_nodes() {
  static hash_map<int, vector<vec2f>> result = {};
  return result;
}

hash_map<int, vector<int>>& debug_indices() {
  static hash_map<int, vector<int>> result = {};
  return result;
}

vector<int>& debug_result() {
  static vector<int> result = {};
  return result;
}
vector<bool>& debug_visited() {
  static vector<bool> result = {};
  return result;
}
vector<int>& debug_stack() {
  static vector<int> result = {};
  return result;
}
bool& debug_restart() {
  static bool result = {};
  return result;
}

}  // namespace yocto
