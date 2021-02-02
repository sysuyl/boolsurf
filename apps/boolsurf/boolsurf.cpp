#include "app.h"
using namespace yocto;

#include <deque>

#ifdef _WIN32
#undef near
#undef far
#endif

// draw with shading
void draw_widgets(app_state* app, const gui_input& input) {
  auto widgets = &app->widgets;
  begin_imgui(widgets, "boolsurf", {0, 0}, {320, 720});

  if (begin_header(widgets, "view")) {
    auto  glmaterial = app->mesh_material;
    auto& params     = app->drawgl_prms;
    draw_checkbox(widgets, "faceted", params.faceted);
    // continue_line(widgets);
    draw_checkbox(widgets, "lines", app->glscene->instances[1]->hidden, true);
    // continue_line(widgets);
    draw_checkbox(widgets, "points", app->glscene->instances[2]->hidden, true);
    draw_coloredit(widgets, "color", glmaterial->color);
    draw_slider(widgets, "resolution", params.resolution, 0, 4096);
    draw_combobox(
        widgets, "lighting", (int&)params.lighting, shade_lighting_names);
    draw_checkbox(widgets, "wireframe", params.wireframe);
    continue_line(widgets);
    draw_checkbox(widgets, "double sided", params.double_sided);
    // draw_slider(widgets, "exposure", params.exposure, -10, 10);
    // draw_slider(widgets, "gamma", params.gamma, 0.1f, 4);
    // draw_slider(widgets, "near", params.near, 0.01f, 1.0f);
    // draw_slider(widgets, "far", params.far, 1000.0f, 10000.0f);
    end_header(widgets);
  }
  if (begin_header(widgets, "inspect")) {
    draw_label(widgets, "shape", app->name);
    draw_label(widgets, "filename", app->filename);
    auto ioshape = app->ioshape;
    draw_label(widgets, "points", std::to_string(ioshape->points.size()));
    draw_label(widgets, "lines", std::to_string(ioshape->lines.size()));
    draw_label(widgets, "triangles", std::to_string(ioshape->triangles.size()));
    draw_label(widgets, "quads", std::to_string(ioshape->quads.size()));
    draw_label(widgets, "positions", std::to_string(ioshape->positions.size()));
    draw_label(widgets, "normals", std::to_string(ioshape->normals.size()));
    draw_label(widgets, "texcoords", std::to_string(ioshape->texcoords.size()));
    draw_label(widgets, "colors", std::to_string(ioshape->colors.size()));
    draw_label(widgets, "radius", std::to_string(ioshape->radius.size()));
    draw_label(widgets, "quads pos", std::to_string(ioshape->quadspos.size()));
    draw_label(
        widgets, "quads norm", std::to_string(ioshape->quadsnorm.size()));
    draw_label(widgets, "quads texcoord",
        std::to_string(ioshape->quadstexcoord.size()));
    end_header(widgets);
  }

  end_imgui(widgets);
}

geodesic_path compute_path(const mesh_polygon& polygon,
    const vector<mesh_point>& points, const bool_mesh& mesh) {
  auto size  = polygon.points.size();
  auto start = polygon.points[size - 2];
  auto end   = polygon.points[size - 1];
  auto path = compute_geodesic_path(mesh, points[start], points[end]);  // check
  return path;
}

void mouse_input(app_state* app, const gui_input& input) {
  if (is_active(&app->widgets)) return;

  if (input.modifier_alt &&
      input.mouse_left.state == gui_button::state::releasing) {
    auto isec = intersect_shape(app, input);
    if (isec.hit) {
      if (app->polygons.size() == 1) app->polygons.push_back(mesh_polygon{});
      if (is_closed(app->polygons.back()))
        app->polygons.push_back(mesh_polygon{});

      auto& polygon = app->polygons.back();
      auto  point   = mesh_point{isec.element, isec.uv};

      if (polygon.points.size()) {
        auto& last_point = app->points[polygon.points.back()];
        if ((point.face == last_point.face) && (point.uv == last_point.uv))
          return;
      }

      app->points.push_back(point);
      polygon.points.push_back(app->points.size() - 1);

      draw_mesh_point(
          app->glscene, app->mesh, app->paths_material, point, 0.0008f);

      if (polygon.points.size() > 1) {
        auto geo_path = compute_path(polygon, app->points, app->mesh);
        auto instance = draw_path(
            app->glscene, app->mesh, app->paths_material, geo_path, 0.0005f);
        app->instances.push_back(instance);

        auto segments = mesh_segments(app->mesh.triangles, geo_path.strip,
            geo_path.lerps, geo_path.start, geo_path.end);

        update_mesh_polygon(polygon, segments);
      }
    }
  }
}

void do_the_thing(app_state* app) {
  // Hashgrid from triangle idx to <polygon idx, edge_idx, segment idx,
  // segment start uv, segment end uv> to handle intersections and
  // self-intersections

  auto hashgrid          = unordered_map<int, vector<hashgrid_entry>>();
  auto intersections     = unordered_map<vec2i, vector<intersection>>();
  auto triangle_segments = unordered_map<int, vector<triangle_segment>>{};

  for (auto p = 0; p < app->polygons.size(); p++) {
    auto& polygon = app->polygons[p];
    for (auto s = 0; s < polygon.segments.size(); s++) {
      auto& segment = polygon.segments[s];
      hashgrid[segment.face].push_back({p, s, segment.start, segment.end});
    }
  }

  for (auto& [face, entries] : hashgrid) {
    for (auto i = 0; i < entries.size() - 1; i++) {
      auto& AB = entries[i];
      for (auto j = i + 1; j < entries.size(); j++) {
        auto& CD = entries[j];

        auto l = intersect_segments(AB.start, AB.end, CD.start, CD.end);
        if (l.x <= 0.0f || l.x >= 1.0f || l.y <= 0.0f || l.y >= 1.0f) continue;

        //        C
        //        |
        // A -- point -- B
        //        |
        //        D

        auto uv       = lerp(CD.start, CD.end, l.y);
        auto point    = mesh_point{face, uv};
        auto point_id = (int)app->points.size();
        app->points.push_back(point);

        auto idx = (int)app->mesh.positions.size();
        auto pos = eval_position(
            app->mesh.triangles, app->mesh.positions, point);
        app->mesh.positions.push_back(pos);

        intersections[{AB.polygon, AB.segment}].push_back({idx, l.x});
        intersections[{CD.polygon, CD.segment}].push_back({idx, l.y});
      }
    }
  }

  for (auto pid = 0; pid < app->polygons.size(); pid++) {
    auto& segments = app->polygons[pid].segments;
    auto  first_id = (int)app->mesh.positions.size();
    auto  id       = first_id;
    for (auto s = 0; s < segments.size(); s++) {
      auto& segment = segments[s];
      auto  uv      = segment.start;
      auto  pos     = eval_position(
          app->mesh.triangles, app->mesh.positions, {segment.face, uv});
      app->mesh.positions.push_back(pos);

      if (intersections.find({pid, s}) != intersections.end()) {
        auto& isecs = intersections[{pid, s}];
        sort(isecs.begin(), isecs.end(),
            [](auto& a, auto& b) { return a.lerp < b.lerp; });

        for (auto& [id1, l] : isecs) {
          auto uv1 = lerp(segment.start, segment.end, l);
          triangle_segments[segment.face].push_back({pid, id, id1, uv, uv1});
          uv = uv1;
          id = id1;
        };
      }

      auto id1 = (int)app->mesh.positions.size();
      if (s == segments.size() - 1) id1 = first_id;

      triangle_segments[segment.face].push_back(
          {pid, id, id1, uv, segment.end});
      id = id1;
    }
  }

  //(marzia) collapse the triangle_segments iterations
  // Not now, they're useful while debugging
  auto face_edgemap = unordered_map<vec2i, vec2i>{};
  for (auto& [face, segments] : triangle_segments) {
    auto [a, b, c] = app->mesh.triangles[face];
    auto nodes     = vector<vec2f>{{0, 0}, {1, 0}, {0, 1}};
    auto indices   = vector<int>{a, b, c};

    for (auto s = 0; s < segments.size(); s++) {
      auto& [p, id, id1, uv, uv1] = segments[s];
      if (find_idx(indices, id) == -1) {
        nodes.push_back(uv);
        indices.push_back(id);
      }

      if (find_idx(indices, id1) == -1) {
        nodes.push_back(uv1);
        indices.push_back(id1);
      }

      auto edge          = make_edge_key({id, id1});
      face_edgemap[edge] = {-1, -1};
    }

    auto triangles = triangulate(nodes);

    for (auto i = 0; i < triangles.size(); i++) {
      auto& [x, y, z] = triangles[i];
      auto i0         = indices[x];
      auto i1         = indices[y];
      auto i2         = indices[z];

      auto triangle_idx = app->mesh.triangles.size();
      app->mesh.triangles.push_back({i0, i1, i2});

      update_face_edgemap(face_edgemap, {i0, i1}, triangle_idx);
      update_face_edgemap(face_edgemap, {i1, i2}, triangle_idx);
      update_face_edgemap(face_edgemap, {i2, i0}, triangle_idx);
    }

    app->mesh.triangles[face] = {0, 0, 0};
  }

  for (auto& [face, segments] : triangle_segments) {
    for (auto i = 0; i < segments.size(); i++) {
      auto& [p, id, id1, uv, uv1] = segments[i];
      auto edge                   = vec2i{id, id1};
      auto edge_key               = make_edge_key(edge);

      auto faces      = face_edgemap[edge_key];
      auto& [a, b, c] = app->mesh.triangles[faces.x];
      if ((edge == vec2i{b, a}) || (edge == vec2i{c, b}) ||
          (edge == vec2i{a, c})) {
        swap(faces.x, faces.y);
      }

      if (faces.x != -1) app->polygons[p].inner_faces.push_back(faces.x);
      if (faces.y != -1) app->polygons[p].outer_faces.push_back(faces.y);
    }
  }

  // Removing face duplicates
  for (auto i = 1; i < app->polygons.size(); i++) {
    auto& inner = app->polygons[i].inner_faces;
    auto& outer = app->polygons[i].outer_faces;

    sort(inner.begin(), inner.end());
    inner.erase(unique(inner.begin(), inner.end()), inner.end());

    sort(outer.begin(), outer.end());
    outer.erase(unique(outer.begin(), outer.end()), outer.end());
  }

  //(marzia) Why do we need this?
  for (auto& ist : app->instances) ist->hidden = true;

  app->mesh.normals = compute_normals(app->mesh.triangles, app->mesh.positions);
  app->mesh.adjacencies = face_adjacencies(app->mesh.triangles);
  set_positions(app->mesh_shape, app->mesh.positions);
  set_triangles(app->mesh_shape, app->mesh.triangles);
  set_normals(app->mesh_shape, app->mesh.normals);
  init_edges_and_vertices_shapes_and_points(app);

  auto tags          = compute_face_tags(app->mesh, app->polygons);
  auto face_polygons = unordered_map<int, vector<int>>();

  auto cells      = vector<vector<int>>();
  auto cell_faces = unordered_map<int, vector<int>>();

  for (auto p = 1; p < app->polygons.size(); p++) {
    auto check = [&](int face, int polygon) {
      return find_in_vec(tags[face], polygon) == -1;
    };

    auto start_out   = app->polygons[p].outer_faces;
    auto visited_out = flood_fill(app->mesh, start_out, -p, check);
    for (auto o : visited_out) face_polygons[o].push_back(-p);

    auto start_in   = app->polygons[p].inner_faces;
    auto visited_in = flood_fill(app->mesh, start_in, p, check);
    for (auto i : visited_in) face_polygons[i].push_back(p);

    // auto color_out = app->cell_materials[(2 * p)]->color;
    // auto color_in  = app->cell_materials[(2 * p) - 1]->color;

    // app->patch_out[p].push_back(app->glscene->instances.size());
    // add_patch_shape(app, visited_out, color_out, 0.0002f * p);

    // app->patch_in[p].push_back(app->glscene->instances.size());
    // add_patch_shape(app, visited_in, color_in, 0.00025f * p);
  }

  // Inverting face_polygons map
  for (auto& [face, polygons] : face_polygons) {
    auto idx = find_idx(cells, polygons);
    if (idx == -1) {
      idx = (int)cells.size();
      cells.push_back(polygons);
    }

    cell_faces[idx].push_back(face);
  }

  for (auto i = 0; i < cells.size(); i++) {
    printf("Cell: %d -> ", i);
    for (auto c : cells[i]) printf("%d ", c);
    printf("\n\t Faces: %d \n", cell_faces[i].size());

    auto color = app->cell_materials[i + 1]->color;
    add_patch_shape(app, cell_faces[i], color, 0.00025f * (i + 1));
  }

  // Previous Implementation
  // auto graph = compute_graph(
  //     app->points.size(), edge_map, counterclockwise);
  // // print_graph(graph);

  // auto edge_info  = compute_edge_info(edge_map, app->polygons);
  // auto components = compute_connected_components(graph);

  // for (auto& component : components) {
  //   auto cells = compute_cells(
  //       component, app->points, app->mesh, app->polygons.size());
  //   // draw_arrangement(app->glscene, app->mesh, app->cell_materials,
  //   //     app->points, arrangement);

  //   print_graph(component);
  //   print_cells(cells);

  //   auto dual_graph = compute_dual_graph(cells, edge_info);
  //   // print_dual_graph(dual_graph);
  //   auto outer_face = compute_outer_face(dual_graph);

  //   visit_dual_graph(dual_graph, cells, outer_face);
  // }

  // Boolean operation example
  // auto ids = vector<int>(arrangement.size(), 0);
  // for (auto i = 0; i < arrangement.size(); i++)
  //   if (!arrangement[i].embedding[1]) ids[i] = 1;

  // polygon_and(arrangement, ids, 0);
  // // polygon_or(arrangement, ids, 2);

  // auto result = vector<cell_polygon>();
  // for (auto i = 0; i < ids.size(); i++) {
  //   if (ids[i]) {
  //     result.push_back(arrangement[i]);
  //     for (auto a : arrangement[i].points) printf("%d ", a);
  //   }
  //   printf("\n");
  // }

  // draw_arrangement(
  //     app->glscene, app->mesh, app->cell_materials, app->points,
  //     result);
}

void key_input(app_state* app, const gui_input& input) {
  for (auto idx = 0; idx < input.key_buttons.size(); idx++) {
    auto button = input.key_buttons[idx];
    if (button.state != gui_button::state::pressing) continue;

    switch (idx) {
      case (int)gui_key('I'): {
        do_the_thing(app);
      } break;
        // case (int)gui_key('L'): {
        //   for (auto& [_, patches] : app->patch_in)
        //     for (auto p : patches)
        //       app->glscene->instances[p]->hidden =
        //           !app->glscene->instances[p]->hidden;
        // } break;
        // case (int)gui_key('R'): {
        //   for (auto& [_, patches] : app->patch_out)
        //     for (auto p : patches)
        //       app->glscene->instances[p]->hidden =
        //           !app->glscene->instances[p]->hidden;
        // } break;
        // case (int)gui_key('N'): {
        //   for (auto& [key, patches] : app->patch_in)
        //     for (auto p : patches)
        //       app->glscene->instances[p]->hidden =
        //           (key == app->current_polygon) ? false : true;

        //   for (auto& [key, patches] : app->patch_out)
        //     for (auto p : patches)
        //       app->glscene->instances[p]->hidden =
        //           (key == app->current_polygon) ? false : true;

        //   app->current_polygon = (app->current_polygon + 1) %
        //                          app->polygons.size();
        // } break;

      case (int)gui_key('S'): {
        save_polygons("test.json", app->points, {{0, 1}, {1, 2, 0}});
      }
      case (int)gui_key('C'): {
        auto old_camera = app->glcamera;
        app->points.clear();
        app->polygons.clear();
        app->polygons.push_back(mesh_polygon{});
        load_shape(app, app->filename);
        clear_scene(app->glscene);
        init_glscene(app, app->glscene, app->mesh, {});
        app->glcamera = old_camera;
      } break;

      case (int)gui_key::enter: {
        auto& polygon = app->polygons.back();
        if (polygon.points.size() < 3 || is_closed(polygon)) return;

        auto point = polygon.points.front();
        polygon.points.push_back(point);

        auto geo_path = compute_path(polygon, app->points, app->mesh);
        draw_path(
            app->glscene, app->mesh, app->paths_material, geo_path, 0.0005f);

        auto segments = mesh_segments(app->mesh.triangles, geo_path.strip,
            geo_path.lerps, geo_path.start, geo_path.end);

        update_mesh_polygon(polygon, segments);

        break;
      }
    }
  }
}

void update_app(const gui_input& input, void* data) {
  auto app = (app_state*)data;

  update_camera(app, input);
  mouse_input(app, input);
  key_input(app, input);

  drop(app, input);

  draw_scene(app, input);
  draw_widgets(app, input);
}

int main(int argc, const char* argv[]) {
  auto app         = new app_state{};
  auto filename    = "tests/_data/shapes/bunny.obj"s;
  auto camera_name = ""s;

  auto window = new gui_window{};

  // parse command line
  auto cli = make_cli("yboolsurf", "views shapes inteactively");
  add_option(cli, "--camera", camera_name, "Camera name.");
  add_option(
      cli, "--resolution,-r", app->drawgl_prms.resolution, "Image resolution.");
  add_option(cli, "--lighting", app->drawgl_prms.lighting, "Lighting type.",
      shade_lighting_names);
  add_option(cli, "shape", filename, "Shape filename", true);
  add_option(cli, "--msaa", window->msaa, "Multisample anti-aliasing.");
  parse_cli(cli, argc, argv);

  init_window(window, {1280 + 320, 720}, "boolsurf", true);
  window->user_data = app;

  app->filename = filename;
  load_shape(app, filename);

  init_glscene(app, app->glscene, app->mesh, {});
  if (window->msaa > 1) set_ogl_msaa();
  set_ogl_blending(true);

  app->widgets = create_imgui(window);

  run_ui(window, update_app);

  // clear
  clear_scene(app->glscene);

  // done
  return 0;
}
