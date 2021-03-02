#include <yocto/yocto_commonio.h>
#include <yocto/yocto_trace.h>

#include "boolsurf.h"
#include "ext/json.hpp"

namespace yocto {

using json = nlohmann::json;
using std::array;

// support for json conversions
inline void to_json(json& j, const vec2f& value) {
  nlohmann::to_json(j, (const array<float, 2>&)value);
}

inline void from_json(const json& j, vec2f& value) {
  nlohmann::from_json(j, (array<float, 2>&)value);
}

inline void to_json(json& j, const frame3f& value) {
  nlohmann::to_json(j, (const array<float, 12>&)value);
}

inline void from_json(const json& j, frame3f& value) {
  nlohmann::from_json(j, (array<float, 12>&)value);
}

inline void to_json(json& js, const mesh_point& value) {
  js["face"] = value.face;
  js["uv"]   = value.uv;
}

inline void from_json(const json& js, mesh_point& value) {
  js.at("face").get_to(value.face);
  js.at("uv").get_to(value.uv);
}

inline void to_json(json& js, const bool_operation& op) {
  js["a"] = op.shape_a;
  js["b"] = op.shape_b;
  // js["type"] = bool_operation::type_names[(int)op_shape];
  js["type"] = (int)op.type;
}

inline void from_json(const json& js, bool_operation& op) {
  js.at("a").get_to(op.shape_a);
  js.at("b").get_to(op.shape_b);
  js.at("type").get_to(op.type);
}

inline void to_json(json& js, const trace_camera& camera) {
  js["frame"]        = camera.frame;
  js["orthographic"] = camera.orthographic;
  js["lens"]         = camera.lens;
  js["film"]         = camera.film;
  js["aspect"]       = camera.aspect;
  js["focus"]        = camera.focus;
  js["aperture"]     = camera.aperture;
}

inline void from_json(const json& js, trace_camera& camera) {
  js.at("frame").get_to(camera.frame);
  js.at("orthographic").get_to(camera.orthographic);
  js.at("lens").get_to(camera.lens);
  js.at("film").get_to(camera.film);
  js.at("aspect").get_to(camera.aspect);
  js.at("focus").get_to(camera.focus);
  js.at("aperture").get_to(camera.aperture);
}
}  // namespace yocto

inline bool load_json(const string& filename, json& js) {
  // error helpers
  auto error = ""s;
  auto text  = ""s;
  if (!load_text(filename, text, error)) {
    printf("[%s]: %s\n", __FUNCTION__, error.c_str());
    return false;
  }
  try {
    js = json::parse(text);
    return true;
  } catch (std::exception& e) {
    printf("[%s]: %s\n", __FUNCTION__, e.what());
    return false;
  }
}

struct bool_test {
  string              model;
  vector<mesh_point>  points;
  vector<vector<int>> polygons;

  vector<bool_operation> operations = {};
  trace_camera           camera     = {};
};

inline bool save_test(const bool_test& test, const string& filename) {
  auto js          = json{};
  js["points"]     = test.points;
  js["polygons"]   = test.polygons;
  js["model"]      = test.model;
  js["operations"] = test.operations;
  js["camera"]     = test.camera;

  auto error = ""s;
  if (!save_text(filename, js.dump(2), error)) {
    printf("[%s]: %s\n", __FUNCTION__, error.c_str());
    return false;
  }
  return true;
}

inline bool load_test(bool_test& test, const string& filename) {
  auto js = json{};
  if (!load_json(filename, js)) {
    return false;
  }

  try {
    test.points   = js["points"].get<vector<mesh_point>>();
    test.polygons = js["polygons"].get<vector<vector<int>>>();
    if (js.find("operations") != js.end()) {
      test.operations = js["operations"].get<vector<bool_operation>>();
    }
    if (js.find("camera") != js.end()) {
      test.camera = js["camera"].get<trace_camera>();
    }
    test.model = js["model"].get<string>();
  } catch (std::exception& e) {
    printf("[%s]: %s\n", __FUNCTION__, e.what());
    return false;
  }
  return true;
}

bool_state state_from_test(const bool_mesh& mesh, const bool_test& test) {
  auto state   = bool_state{};
  state.points = test.points;
  for (auto& polygon : test.polygons) {
    // Add new polygon to state.
    auto& mesh_polygon  = state.polygons.emplace_back();
    mesh_polygon.points = polygon;
  }

  for (auto& mesh_polygon : state.polygons) {
    for (int i = 0; i < mesh_polygon.points.size(); i++) {
      auto start = mesh_polygon.points[i];
      auto end   = mesh_polygon.points[(i + 1) % mesh_polygon.points.size()];
      auto path  = compute_geodesic_path(
          mesh, state.points[start], state.points[end]);
      auto segments = mesh_segments(
          mesh.triangles, path.strip, path.lerps, path.start, path.end);
      mesh_polygon.segments += segments;
    }
  }
  return state;
}

#include <yocto/yocto_color.h>

inline string tree_to_string(const vector<mesh_cell>& cells) {
  string result = "digraph {\n";
  result += "forcelabels=true\n";

  auto get_cell_color = [](const mesh_cell& cell, int cell_id) {
    auto color = vec3f{0, 0, 0};
    int  count = 0;
    for (int p = 0; p < cell.labels.size(); p++) {
      auto label = cell.labels[p];
      if (label > 0) {
        color += get_color(p);
        count += 1;
      }
    }
    if (count > 0) {
      color /= count;
      color += vec3f{1, 1, 1} * 0.1f * yocto::sin(cell_id);
    } else {
      color = {0.9, 0.9, 0.9};
    }
    return color;
  };

  for (int i = 0; i < cells.size(); i++) {
    auto& cell  = cells[i];
    auto  color = rgb_to_hsv(get_cell_color(cell, i));
    char  str[1024];
    sprintf(str, "%d [label=\"%d\" style=filled fillcolor=\"%f %f %f\"]\n", i,
        i, color.x, color.y, color.z);
    result += std::string(str);

    for (auto [neighbor, polygon] : cell.adjacency) {
      if (polygon < 0) continue;
      int  c     = neighbor;
      auto color = rgb_to_hsv(get_color(polygon));
      sprintf(str, "%d -> %d [ label=\"%d\" color=\"%f %f %f\"]\n", i, c,
          polygon, color.x, color.y, color.z);
      result += std::string(str);
    }
  }
  result += "}\n";
  return result;
}

inline void save_tree_png(
    const bool_state& state, string filename, const string& extra = "") {
  if (filename.empty()) filename = "data/tests/test.json";
  auto  graph = replace_extension(filename, extra + ".txt");
  FILE* file  = fopen(graph.c_str(), "w");
  fprintf(file, "%s", tree_to_string(state.cells).c_str());
  fclose(file);

  auto image = replace_extension(filename, extra + ".png");
  auto cmd   = "dot -Tpng "s + graph + " > " + image;
  printf("%s\n", cmd.c_str());
  system(cmd.c_str());
}