#include <yocto/yocto_commonio.h>

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

// support for json conversions
inline void to_json(json& js, const mesh_point& value) {
  js["face"] = value.face;
  js["uv"]   = value.uv;
}

inline void from_json(const json& js, mesh_point& value) {
  js.at("face").get_to(value.face);
  js.at("uv").get_to(value.uv);
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
};

inline bool save_test(const bool_test& test, const string& filename) {
  auto js        = json{};
  js["points"]   = test.points;
  js["polygons"] = test.polygons;
  js["model"]    = test.model;

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
    test.model    = js["model"].get<string>();
  } catch (std::exception& e) {
    printf("[%s]: %s\n", __FUNCTION__, e.what());
    return false;
  }
  return true;
}