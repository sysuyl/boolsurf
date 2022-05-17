#define NOMINMAX

// TODO(giacomo): review usless includes
#include <yocto/yocto_bvh.h>
//#include <yocto/yocto_common.h>
#include <boolsurf/boolsurf.h>
#include <boolsurf/boolsurf_io.h>
#include <cinolib/homotopy_basis.h>
#include <yocto/yocto_geometry.h>
#include <yocto/yocto_image.h>
#include <yocto/yocto_mesh.h>
#include <yocto/yocto_sceneio.h>
#include <yocto_gui/yocto_imgui.h>
#include <yocto_gui/yocto_window.h>

#include <unordered_map>
#include <unordered_set>

#include "render.h"

namespace yocto {

// Application state
struct app_state {
  // loading parameters
  string         model_filename       = "";
  string         test_filename        = "";
  string         output_test_filename = "data/tests/test.json";
  string         svg_filename         = "";
  int            svg_subdivs          = 2;
  float          drawing_size         = 0.01;
  bool           project_points       = false;
  int            num_sampled_polygons = 200;
  bool_test      test                 = {};
  bool_operation operation            = {};
  gui_window*    window               = nullptr;

  bool  thick_lines = false;
  float line_width  = 0.003f;

  bool color_shapes      = false;
  bool save_edges        = false;
  bool save_polygons     = false;
  bool save_generators   = false;
  bool smooth_generators = false;

  bool color_hashgrid  = false;
  bool show_polygons   = true;
  bool show_generators = true;
  // bool         show_arrows    = false;
  bool         use_projection = false;
  scene_camera camera         = {};

  // options
  shade_params drawgl_prms = {};

  // boolmesh info
  bool_mesh mesh          = bool_mesh{};
  bool_mesh mesh_original = bool_mesh{};

  bool_state state = {};

  // TODO (marzia): è ridicolo, cambialo!
  vector<bool_shape_shape> shape_shapes      = {};
  vector<shade_instance*>  generators_shapes = {};
  vector<shade_instance*>  cell_shapes       = {};
  // vector<shade_instance*> polygon_shapes = {};
  // vector<shade_instance*> arrow_shapes        = {};
  shade_instance*         hashgrid_shape      = nullptr;
  vector<shade_instance*> border_faces_shapes = {};

  // shade_instance*         inner_faces_shape = nullptr;
  // shade_instance*         outer_faces_shape = nullptr;

  vector<bool_state> history               = {};
  int                history_index         = -1;
  int                selected_point        = -1;
  float              selected_point_radius = 0;
  int                selected_cell         = -1;
  int                selected_shape        = 0;
  int                selected_polygon      = 0;
  int                selected_strategy     = 0;
  int                count                 = 0;
  // int                selected_polygon      = 0;
  hash_set<int> current_shape = {};

  // rendering state
  shade_scene*    glscene           = new shade_scene{};
  shade_camera*   glcamera          = nullptr;
  shade_instance* mesh_instance     = nullptr;
  shade_instance* edges_instance    = nullptr;
  shade_instance* vertices_instance = nullptr;
  shade_instance* temp_patch        = nullptr;

  shade_material* mesh_material   = nullptr;
  shade_material* edges_material  = nullptr;
  shade_material* points_material = nullptr;
  shade_material* paths_material  = nullptr;
  shade_material* isecs_material  = nullptr;

  struct {
    shade_material* red   = nullptr;
    shade_material* blue  = nullptr;
    shade_material* green = nullptr;
    shade_material* white = nullptr;
  } materials;

  vector<shade_material*> cell_materials       = {};
  vector<shade_instance*> instances            = {};
  shade_instance*         selected_point_shape = nullptr;

  //(marzia) Useful while debugging!
  vector<int> cell_patches   = {};
  int         current_patch  = 0;
  int         current_border = 0;
  int         current_cell   = 0;

  gui_widgets widgets = {};

  mesh_point last_clicked_point          = {};
  vec3f      last_clicked_point_position = {};
  mesh_point last_clicked_point_original = {};

  struct last_svg {
    vector<Svg_Shape> svg;
    int               last_shape;
  } last_svg          = {};
  bool_test temp_test = {};

  ~app_state() {
    if (glscene) delete glscene;
  }
};

shade_instance* get_polygon_shape(
    app_state* app, const vector<int>& polygon, int index) {
  auto polygon_shape    = add_shape(app->glscene);
  auto polygon_material = add_material(app->glscene);

  polygon_material->color = get_color(index);

  auto polygon_instance = add_instance(
      app->glscene, identity3x4f, polygon_shape, polygon_material);

  set_polygon_shape(
      polygon_instance, app->mesh, polygon, app->thick_lines, app->line_width);
  return polygon_instance;
}

shade_instance* get_polygon_shape(
    app_state* app, const mesh_polygon& polygon, int index) {
  auto polygon_shape    = add_shape(app->glscene);
  auto polygon_material = add_material(app->glscene);

  polygon_material->color = get_color(index);

  auto polygon_instance = add_instance(
      app->glscene, identity3x4f, polygon_shape, polygon_material);

  if (polygon.length > 0)
    set_polygon_shape(polygon_instance, app->mesh, polygon, app->thick_lines,
        app->line_width);

  // polygon_instance->depth_test = ogl_depth_test::always;
  // app->polygon_shapes += polygon_instance;

  // auto arrow_shape    = add_shape(app->glscene);
  // auto arrow_material = polygon_material;

  // auto arrow_instance = add_instance(
  //     app->glscene, identity3x4f, arrow_shape, arrow_material);
  // arrow_instance->depth_test = ogl_depth_test::always;
  // app->arrow_shapes += arrow_instance;

  return polygon_instance;
}

void set_polygon_shape(app_state* app, int shape_id, int polygon_id) {
  auto& mesh      = app->mesh;
  auto  positions = polygon_positions(
       app->state.bool_shapes[shape_id].polygons[polygon_id], mesh);
  auto normals = polygon_normals(
      app->state.bool_shapes[shape_id].polygons[polygon_id], mesh);

  auto polygon_shape = make_polygon_shape(
      mesh, positions, app->thick_lines, app->line_width);
  set_shape(app->shape_shapes[shape_id].polygons[polygon_id], polygon_shape);

  // auto arrow_shape = make_arrow_shape();
  // for (auto& p : arrow_shape.positions) p *= 0.005;
  // for (int i :
  //     {1, (int)positions.size() / 3, (int)(positions.size() * 2) / 3}) {
  //   arrow_shape.froms.push_back(positions[i]);
  //   auto left = normalize(positions[i + 1] - positions[i]);
  //   auto y    = cross(left, normals[i]);
  //   arrow_shape.tos.push_back(positions[i] - y);
  // }

  // set_shape(app->arrow_shapes[polygon_id], arrow_shape);
}

// void set_polygon_shape(app_state* app, int polygon_id) {
//   auto& mesh      = app->mesh;
//   auto  positions = polygon_positions(app->state.polygons[polygon_id], mesh);
//   auto  normals   = polygon_normals(app->state.polygons[polygon_id], mesh);

//   auto polygon_shape = make_polygon_shape(
//       mesh, positions, app->thick_lines, app->line_width);
//   set_shape(app->polygon_shapes[polygon_id], polygon_shape);

//   // auto arrow_shape = make_arrow_shape();
//   // for (auto& p : arrow_shape.positions) p *= 0.005;
//   // for (int i :
//   //     {1, (int)positions.size() / 3, (int)(positions.size() * 2) / 3}) {
//   //   arrow_shape.froms.push_back(positions[i]);
//   //   auto left = normalize(positions[i + 1] - positions[i]);
//   //   auto y    = cross(left, normals[i]);
//   //   arrow_shape.tos.push_back(positions[i] - y);
//   // }

//   // set_shape(app->arrow_shapes[polygon_id], arrow_shape);
// }

void update_polygon(
    app_state* app, int shape_id, int polygon_id, int point_id = 0) {
  auto& polygon       = app->state.bool_shapes[shape_id].polygons[polygon_id];
  auto& polygon_shape = app->shape_shapes[shape_id].polygons[polygon_id];

  recompute_polygon_segments(app->mesh, app->state, polygon, point_id);
  if (!app->show_polygons) return;

  if (polygon.length > 0) {
    set_polygon_shape(app, shape_id, polygon_id);
  } else if (polygon_shape->shape) {
    clear_shape(polygon_shape->shape);
  }
}

// void update_polygon(app_state* app, int polygon_id, int index = 0) {
//   auto& mesh_polygon = app->state.polygons[polygon_id];
//   // app->polygon_shapes.resize(app->state.polygons.size());
//   auto& polygon_shape = app->polygon_shapes[polygon_id];
//   // auto& arrow_shape   = app->arrow_shapes[polygon_id];

//   // TODO(giacomo): Solve this situation.
//   // if (!polygon_shape) {
//   //   polygon_shape = new shade_instance{};
//   // }

//   // if (!polygon_shape->shape) {
//   //   polygon_shape->shape = new shade_shape{};
//   // }

//   // Draw polygon.
//   recompute_polygon_segments(app->mesh, app->state, mesh_polygon, index);
//   if (!app->show_polygons) return;

//   if (mesh_polygon.length > 0) {
//     set_polygon_shape(app, polygon_id);
//   } else if (polygon_shape->shape) {
//     clear_shape(polygon_shape->shape);
//   }
// }

void update_polygons(app_state* app) {
  for (auto s = 0; s < app->state.bool_shapes.size(); s++) {
    for (int p = 0; p < app->state.bool_shapes[s].polygons.size(); p++) {
      update_polygon(app, s, p);
    }
  }

  for (auto s = app->state.bool_shapes.size(); s < app->shape_shapes.size();
       s++) {
    for (auto p = 0; p < app->shape_shapes[s].polygons.size(); p++) {
      clear_shape(app->shape_shapes[s].polygons[p]->shape);
    }
  }

  app->shape_shapes.resize(app->state.bool_shapes.size());
}

void commit_state(app_state* app) {
  PRINT_CALL();
  app->history_index += 1;
  app->history.resize(app->history_index + 1);
  app->history[app->history_index] = app->state;
}

inline void update_cell_shapes(app_state* app);
inline void update_cell_colors(app_state* app);

bool undo_state(app_state* app) {
  if (app->history_index <= 0) return false;
  app->history_index -= 1;
  app->state = app->history[app->history_index];
  update_polygons(app);
  update_cell_shapes(app);
  update_cell_colors(app);
  return true;
}
bool redo_state(app_state* app) {
  if (app->history_index >= app->history.size() - 1) return false;
  app->history_index += 1;
  app->state = app->history[app->history_index];
  update_polygons(app);
  update_cell_shapes(app);
  update_cell_colors(app);
  return true;
}

void load_shape(app_state* app, const string& filename) {
  app->model_filename = filename;

  auto error = ""s;
  if (!load_shape(filename, app->mesh, error)) {
    printf("%s\n", error.c_str());
    print_fatal("Error loading model " + filename);
  }

  auto modelname    = path_basename(app->model_filename);
  auto basefilename = "data/homotopy/"s + modelname + "_basis.json"s;
  init_mesh(app->mesh);

  if (path_exists(basefilename)) {
    printf("Loading: %s\n", basefilename.c_str());
    load_homotopy_basis(app->mesh, basefilename);

    init_mesh(app->mesh);
    app->mesh_original = app->mesh;

    app->mesh.homotopy_basis.smooth_basis = smooth_homotopy_basis(
        app->mesh.homotopy_basis, app->mesh, app->smooth_generators);

    for (auto b = 0; b < app->mesh.homotopy_basis.smooth_basis.size(); b++) {
      printf("Base: %d - size: %d\n", b,
          app->mesh.homotopy_basis.smooth_basis[b].length);
      auto& smooth_base = app->mesh.homotopy_basis.smooth_basis[b];
      auto  base_shape  = get_polygon_shape(app, smooth_base, 13);
      app->generators_shapes.push_back(base_shape);
    }
  }

  app->mesh_original = app->mesh;
}

void init_edges_and_vertices_shapes_and_points(
    app_state* app, bool thin = true) {
  auto  edges           = get_edges(app->mesh.triangles, {});
  float avg_edge_length = 0;
  for (auto& edge : edges) {
    avg_edge_length += length(
        app->mesh.positions[edge.x] - app->mesh.positions[edge.y]);
  }
  avg_edge_length /= edges.size();
  auto cylinder_radius = 0.01f * avg_edge_length;

  auto edges_shape = app->edges_instance->shape;
  if (thin) {
    set_quads(edges_shape, {});
    set_positions(edges_shape, app->mesh.positions);
    set_lines(edges_shape, edges);
    set_normals(edges_shape, {});
    set_texcoords(edges_shape, {});
    set_instances(edges_shape, {});
  } else {
    //    auto cylinder = make_uvcylinder({8, 1, 1}, {cylinder_radius, 1});
    //    for (auto& p : cylinder.positions) {
    //      p.z = p.z * 0.5 + 0.5;
    //    }
    //    set_quads(edges_shape, cylinder.quads);
    //    set_positions(edges_shape, cylinder.positions);
    //    set_normals(edges_shape, cylinder.normals);
    //    set_texcoords(edges_shape, cylinder.texcoords);
    //    set_instances(edges_shape, froms, tos);
  }

  auto vertices_shape  = app->vertices_instance->shape;
  auto vertices_radius = 3.0f * cylinder_radius;
  auto vertices        = make_sphere(3, vertices_radius);
  set_quads(vertices_shape, vertices.quads);
  set_positions(vertices_shape, vertices.positions);
  set_normals(vertices_shape, vertices.normals);
  set_texcoords(vertices_shape, vertices.texcoords);
  set_instances(vertices_shape, app->mesh.positions);
  // vertices_shape  = add_shape(glscene, {}, {}, {}, vertices.quads,
  //     vertices.positions, vertices.normals, vertices.texcoords, {});
  // set_instances(vertices_shape, app->mesh.positions);
}

void init_glscene(app_state* app, shade_scene* glscene, const bool_mesh& mesh) {
  // init scene
  init_scene(glscene, true);

  // TODO(fabio): move this function to math
  auto camera_frame = [](float lens, float aspect, float film = 0.036) {
    auto camera_dir  = normalize(vec3f{0, 0.5, 1});
    auto bbox_radius = 2.0f;
    auto camera_dist = bbox_radius * lens / (film / aspect);
    return lookat_frame(camera_dir * camera_dist, {0, 0, 0}, {0, 1, 0});
  };

  // camera
  //  if (progress_cb) progress_cb("convert camera", progress.x++, progress.y);
  app->glcamera = add_camera(glscene, camera_frame(0.050, 16.0f / 9.0f, 0.036),
      0.050, 16.0f / 9.0f, 0.036);
  app->glcamera->focus = length(app->glcamera->frame.o);

  app->mesh_material  = add_material(glscene, {0, 0, 0}, {1, 1, 1}, 1, 0, 0.4);
  app->edges_material = add_material(
      glscene, {0, 0, 0}, {0.4, 0.4, 0.4}, 1, 0, 0.4);
  app->points_material = add_material(glscene, {0, 0, 0}, {0, 0, 1}, 1, 0, 0.4);
  app->paths_material  = add_material(glscene, {0, 0, 0}, {1, 1, 1}, 1, 0, 0.4);
  app->isecs_material  = add_material(glscene, {0, 0, 0}, {0, 1, 0}, 1, 0, 0.4);

  set_unlit(app->edges_material, true);
  set_unlit(app->paths_material, true);

  auto colors = vector<vec3f>{
      {0.5, 0.5, 0.5},  // REMEMBER IT
      {1, 0, 0},
      {0, 0.5, 0},
      {0, 0, 1},
      {0, 0.5, 0.5},
      {1, 0.5, 0},
      {0.5, 0, 1},
      {0.5, 0, 0},
      {0, 0.5, 0},
      {0, 0, 0.5},
      {0, 0.5, 0.5},
      {0.5, 0.5, 0},
      {0.5, 0, 0.5},
  };
  app->cell_materials.resize(colors.size());
  for (int i = 0; i < colors.size(); i++) {
    app->cell_materials[i] = add_material(
        glscene, {0, 0, 0}, colors[i], 1, 0, 0.4);
  }

  app->materials.red   = add_material(glscene, {0, 0, 0}, {1, 0, 0}, 1, 0, 0.4);
  app->materials.green = add_material(glscene, {0, 0, 0}, {0, 1, 0}, 1, 0, 0.4);
  app->materials.blue  = add_material(glscene, {0, 0, 0}, {0, 0, 1}, 1, 0, 0.4);
  app->materials.white = add_material(glscene, {0, 0, 0}, {1, 1, 1}, 1, 0, 0.4);

  // shapes
  //  if (progress_cb) progress_cb("convert shape", progress.x++, progress.y);
  auto mesh_shape = add_shape(glscene, {}, {}, app->mesh.triangles, {},
      app->mesh.positions, app->mesh.normals, {}, {}, true);

  app->drawgl_prms.faceted = true;
  set_instances(mesh_shape, {}, {});
  app->drawgl_prms.background = {1, 1, 1, 1};

  auto edges_shape    = add_shape(glscene);
  auto vertices_shape = add_shape(glscene);
  app->mesh_instance  = add_instance(
       glscene, identity3x4f, mesh_shape, app->mesh_material);
  app->edges_instance = add_instance(
      glscene, identity3x4f, edges_shape, app->edges_material, true);
  app->vertices_instance = add_instance(
      glscene, identity3x4f, vertices_shape, app->points_material, true);
  init_edges_and_vertices_shapes_and_points(app);
}

// draw with shape
void draw_scene(const app_state* app, const gui_input& input) {
  draw_scene(app->glscene, app->glcamera, input.framebuffer_viewport,
      app->drawgl_prms);
}

void update_camera(app_state* app, const gui_input& input) {
  if (is_active(&app->widgets)) return;

  // handle mouse and keyboard for navigation
  if ((input.mouse_left || input.mouse_right) && !input.modifier_alt) {
    auto dolly  = 0.0f;
    auto pan    = zero2f;
    auto rotate = zero2f;
    if (input.mouse_left && !input.modifier_shift)
      rotate = (input.mouse_pos - input.mouse_last) / 100.0f;
    if (input.mouse_right)
      dolly = (input.mouse_pos.x - input.mouse_last.x) / 100.0f;
    if (input.mouse_left && input.modifier_shift)
      pan = (input.mouse_pos - input.mouse_last) / 100.0f;
    pan.x    = -pan.x;
    rotate.y = -rotate.y;

    std::tie(app->glcamera->frame, app->glcamera->focus) = camera_turntable(
        app->glcamera->frame, app->glcamera->focus, rotate, dolly, pan);
  }
  app->camera.frame    = app->glcamera->frame;
  app->camera.lens     = app->glcamera->lens;
  app->camera.film     = app->glcamera->film;
  app->camera.aspect   = app->glcamera->aspect;
  app->camera.aperture = app->glcamera->aperture;
  app->camera.focus    = app->glcamera->focus;
};

void drop(app_state* app, const gui_input& input) {
  if (input.dropped.size()) {
    app->model_filename = input.dropped[0];
    clear_scene(app->glscene);
    load_shape(app, app->model_filename);
    init_glscene(app, app->glscene, app->mesh);
    return;
  }
}

mesh_point intersect_mesh(const app_state* app, const gui_input& input) {
  auto uv = vec2f{input.mouse_pos.x / float(input.window_size.x),
      input.mouse_pos.y / float(input.window_size.y)};
  return intersect_mesh(app->mesh, app->camera, uv);
}

mesh_point intersect_mesh_original(
    const app_state* app, const gui_input& input) {
  auto uv = vec2f{input.mouse_pos.x / float(input.window_size.x),
      input.mouse_pos.y / float(input.window_size.y)};
  return intersect_mesh(app->mesh_original, app->camera, uv);
}

shade_instance* add_patch_shape(
    app_state* app, const vector<int>& faces, const vec3f& color) {
  auto instance       = add_instance(app->glscene);
  instance->shape     = add_shape(app->glscene, {}, {}, {}, {}, {}, {}, {}, {});
  instance->material  = add_material(app->glscene);
  *instance->material = *app->mesh_material;
  instance->material->color = color;
  set_patch_shape(instance->shape, app->mesh, faces);
  return instance;
}

// shade_instance* add_patch_shape(
//     app_state* app, const vector<int>& faces, shade_material* material) {
//   auto patch_shape = add_shape(app->glscene, {}, {}, {}, {}, {}, {}, {}, {});
//   set_patch_shape(patch_shape, app->mesh, faces);
//   return add_instance(app->glscene, identity3x4f, patch_shape, material);
// }
// void add_polygon_shape(app_state* app, const mesh_polygon& polygon, int
// index) {
//   add_polygon_shape(app, polygon, index);
// }

void add_shape_shape(app_state* app, int shape_id) {
  auto& shape        = app->state.bool_shapes[shape_id];
  auto& shape_shapes = app->shape_shapes.emplace_back();

  for (auto& polygon : shape.polygons) {
    auto p_shape = get_polygon_shape(app, polygon, shape_id);
    shape_shapes.polygons.push_back(p_shape);
  }
}

inline void update_cell_shapes(app_state* app) {
  for (int i = 0; i < app->state.cells.size(); i++) {
    auto& cell = app->state.cells[i];
    set_patch_shape(app->cell_shapes[i]->shape, app->mesh, cell.faces);
  }
}

inline void update_cell_colors(app_state* app) {
  auto& state = app->state;
  if (app->test.cell_colors.empty() || app->color_shapes) {
    app->test.cell_colors.resize(state.cells.size());
    for (int i = 0; i < state.cells.size(); i++) {
      app->test.cell_colors[i] = get_cell_color(state, i, app->color_shapes);
    }
  }

  for (int i = 0; i < state.cells.size(); i++) {
    app->cell_shapes[i]->material->color = app->test.cell_colors[i];
  }
}

void update_svg(app_state* app) {
  init_from_svg(app->state, app->mesh, app->last_clicked_point,
      app->last_svg.svg, app->drawing_size, app->svg_subdivs);

  if (app->state.bool_shapes.empty()) return;
  auto last_shape = app->state.bool_shapes.size() - 1;
  add_shape_shape(app, last_shape);

  update_polygons(app);
}

void save_test(
    app_state* app, const bool_state& state, const string& filename) {
  app->test.points     = state.points;
  app->test.polygons   = {{}};
  app->test.shapes     = {{0}};
  app->test.are_closed = {false};

  for (auto s = 0; s < state.bool_shapes.size(); s++) {
    auto test_shape = vector<int>();
    for (auto p = 0; p < state.bool_shapes[s].polygons.size(); p++) {
      auto& polygon = state.bool_shapes[s].polygons[p];

      if (polygon.points.size()) {
        test_shape.push_back((int)app->test.polygons.size());
        app->test.polygons.push_back(polygon.points);
        app->test.are_closed.push_back(polygon.is_closed);
      }
    }

    if (test_shape.size()) app->test.shapes.push_back(test_shape);
  }

  if (app->test.are_closed.size() == 1) app->test.are_closed = {};

  app->test.camera.frame    = app->glcamera->frame;
  app->test.camera.lens     = app->glcamera->lens;
  app->test.camera.aspect   = app->glcamera->aspect;
  app->test.camera.film     = app->glcamera->film;
  app->test.camera.aperture = app->glcamera->aperture;
  app->test.camera.focus    = app->glcamera->focus;

  save_test(app->test, filename);
}

void init_from_test(app_state* app) {
  app->state.polygons.clear();
  app->state.points = app->test.points;

  if (app->test.has_camera) {
    app->glcamera->frame    = app->test.camera.frame;
    app->glcamera->lens     = app->test.camera.lens;
    app->glcamera->aspect   = app->test.camera.aspect;
    app->glcamera->film     = app->test.camera.film;
    app->glcamera->aperture = app->test.camera.aperture;
    app->glcamera->focus    = app->test.camera.focus;
  }

  // Init bool_state
  // if (app->test.screenspace) {
  //   app->state = state_from_screenspace_test(
  //       app->mesh, app->test, app->drawing_size, app->use_projection);
  // } else {
  //   app->state = state_from_test(
  //       app->mesh, app->test, app->drawing_size, app->use_projection);
  // }
  app->state = state_from_test(
      app->mesh, app->test, app->drawing_size, app->use_projection);

  for (auto s = 0; s < app->state.bool_shapes.size(); s++) {
    add_shape_shape(app, s);
  }
}

}  // namespace yocto
