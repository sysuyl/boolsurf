#include <boolsurf/boolsurf.h>
#include <boolsurf/boolsurf_io.h>
#include <yocto/yocto_cli.h>
#include <yocto/yocto_sampling.h>
#include <yocto/yocto_sceneio.h>
#include <yocto/yocto_shape.h>
#include <yocto/yocto_trace.h>

#include "serialize/serialize.h"
using namespace yocto;

struct test_stats {
  string model_filename   = ""s;
  int    model_triangles  = 0;
  int    genus            = 0;
  double triangulation_ns = 0.0;
  double flood_fill_ns    = 0.0;
  double labelling_ns     = 0.0;
  double boolean_ns       = 0.0;
  double total_ns         = 0.0;
  int    polygons         = 0;
  int    control_points   = 0;
  int    added_points     = 0;
  int    sliced_triangles = 0;
  int    added_triangles  = 0;
  int    graph_nodes      = 0;
  int    graph_edges      = 0;
  int    graph_cycles     = 0;
};

void save_image(const string& output_image_filename, const bool_mesh& mesh,
    const bool_state& state, const scene_camera& camera, bool color_shapes,
    int spp) {
  if (output_image_filename == "no-output") return;
  auto scene = make_scene(mesh, state, camera, color_shapes);

  auto params    = trace_params{};
  auto error     = string{};
  params.sampler = trace_sampler_type::eyelight;
  params.samples = spp;
  auto image     = trace_image(scene, params);
  save_image(output_image_filename, image, error);
}

void save_image(
    const string& output_image_filename, const scene_model& scene, int spp) {
  if (output_image_filename == "no-output") return;
  auto params    = trace_params{};
  auto error     = string{};
  params.sampler = trace_sampler_type::eyelight;
  params.samples = spp;
  auto image     = trace_image(scene, params);
  save_image(output_image_filename, image, error);
}

vector<mesh_cell> make_cell_graph_fruit() {
  auto graph = vector<mesh_cell>();
  auto cell0 = mesh_cell();
  cell0.adjacency.insert({1, 2});
  cell0.adjacency.insert({3, 1});
  graph.push_back(cell0);

  auto cell1 = mesh_cell();
  cell1.adjacency.insert({2, 1});
  cell1.adjacency.insert({0, -2});
  graph.push_back(cell1);

  auto cell2 = mesh_cell();
  cell2.adjacency.insert({1, -1});
  cell2.adjacency.insert({3, -2});
  graph.push_back(cell2);

  auto cell3 = mesh_cell();
  cell3.adjacency.insert({2, 2});
  cell3.adjacency.insert({0, -1});
  graph.push_back(cell3);
  return graph;
}

bool check_ambient_cell(const bool_state& state) {
  if (state.cells.size() == 1) return false;

  // Getting max faces in ambient cell
  auto zero               = vector<int>(state.labels[0].size(), 0);
  auto ambient_cell_faces = 0;
  for (int c = 0; c < state.cells.size(); c++) {
    auto& cell = state.cells[c];
    if (state.labels[c] != zero) continue;
    ambient_cell_faces = max(ambient_cell_faces, (int)cell.faces.size());
  }

  printf("ambient_cell_faces: %d\n", ambient_cell_faces);
  for (auto& cell : state.cells) {
    if (cell.faces.size() > ambient_cell_faces) return false;
  }

  return true;
}

bool check_cell_adjacency(const vector<mesh_cell>& cells) {
  auto correct_cells = make_cell_graph_fruit();
  if (cells.size() != correct_cells.size()) return false;

  auto cells_degree         = hash_map<int, int>();
  auto correct_cells_degree = hash_map<int, int>();

  for (auto c = 0; c < cells.size(); c++) {
    auto cell_degree = cells[c].adjacency.size();
    cells_degree[cell_degree] += 1;

    auto correct_cell_degree = correct_cells[c].adjacency.size();
    correct_cells_degree[correct_cell_degree] += 1;
  }

  for (auto& [key, value] : cells_degree) {
    if (!contains(correct_cells_degree, key)) return false;
    if (cells_degree[key] != correct_cells_degree[key]) return false;
  }

  // (Marzia): Aggiungi altri check sul grafo
  return true;
}

// TODO(marzia): move to boolsurf_io?
bool_state state_from_screenspace_test(
    bool_mesh& mesh, bool_test& test, float drawing_size, bool use_projection) {
  int  seed          = 0;
  auto stop          = false;
  auto rng           = make_rng(seed);
  auto state         = bool_state{};
  auto mesh_original = mesh;

  while (!stop) {
    state    = {};
    auto cam = scene_camera{};
    auto eye = sample_sphere(rand2f(rng)) * 3.5;

    auto position = vec3f{0, 0, 0};
    cam.frame     = lookat_frame(eye, position, {0, 1, 0});
    cam.focus     = length(eye - position);

    auto center = intersect_mesh(mesh, cam, vec2f{0.5, 0.5});
    test.camera = make_camera(mesh, seed);

    add_polygons(state, mesh, test.camera, test, center, drawing_size, false);
    test.camera = cam;

    try {
      auto timer = print_timed("[compute_cells]");
      compute_cells(mesh, state);
      compute_shapes(state);

      stop = check_ambient_cell(state);
      stop = stop && check_cell_adjacency(state.cells);
    } catch (const std::exception&) {
      stop = true;
    }

    mesh = mesh_original;
    if (stop) break;
  }
  return state;
}

int main(int num_args, const char* args[]) {
  auto test_filename         = ""s;
  auto output_image_filename = "data/render.png"s;
  auto output_scene_filename = ""s;
  auto spp                   = 4;
  auto model_filename        = ""s;
  auto svg_filename          = ""s;
  auto svg_subdivs           = 2;
  auto drawing_size          = 0.005f;
  auto color_shapes          = false;
  auto stats_filename        = ""s;
  auto append_stats          = false;
  auto stats                 = test_stats{};

  // parse command line
  auto cli = make_cli("test", "test boolsurf algorithms");
  add_argument(
      cli, "input", test_filename, "Input test filename (.json).", {}, false);
  add_option(cli, "output_image", output_image_filename,
      "Output image filename (.png).");
  add_option(cli, "output_scene", output_scene_filename, "");
  add_option(cli, "spp", spp, "Samples per pixel.");
  add_option(cli, "model", model_filename, "Input model filename.");
  add_option(cli, "svg", svg_filename, "Input svg filename.");
  add_option(cli, "svg-subdivs", svg_subdivs, "Svg subdivisions.");
  add_option(cli, "drawing-size", drawing_size, "Size of mapped drawing.");
  add_option(cli, "color-shapes", color_shapes, "Color shapes.");
  add_option(cli, "stats", stats_filename, "output stats");
  add_option(cli, "append-stats", append_stats, "append statistics");
  parse_cli(cli, num_args, args);

  if (!test_filename.size()) print_fatal("No input filename");

  auto test      = bool_test{};
  auto extension = path_extension(test_filename);
  if (extension == ".svg") {
    auto script_path = normalize_path("scripts/svg_parser.py"s);
    auto output      = normalize_path("data/tests/tmp.json"s);
    auto cmd         = "python3 "s + script_path + " "s + test_filename + " "s +
               output + " "s + to_string(svg_subdivs);
    auto ret_value = system(cmd.c_str());
    if (ret_value != 0) print_fatal("Svg conversion failed " + test_filename);

    test_filename = output;
  }

  if (!load_test(test, test_filename)) {
    print_fatal("Error loading test " + test_filename);
  }

  if (model_filename.size()) test.model = model_filename;
  stats.model_filename = test.model;

  if (stats_filename.size() && !append_stats) {
    auto stats_file = fopen(stats_filename.c_str(), "w");
    fprintf(stats_file,
        "model, model_triangles, genus, triangulation_ns, flood_fill_ns, labelling_ns, boolean_ns, total_ns, polygons, control_points, added_points, sliced_triangles, added_triangles, graph_nodes, graph_edges, graph_cycles\n");
    fclose(stats_file);
  }

  // Init mesh.
  auto error = string{};
  auto mesh  = bool_mesh{};
  {
    if (!load_shape(test.model, mesh, error)) {
      printf("%s\n", error.c_str());
      print_fatal("Error loading model " + test_filename);
    }

    init_mesh(mesh);
    stats.model_triangles = (int)mesh.triangles.size();
    printf("triangles: %d\n", (int)mesh.triangles.size());
    printf("positions: %d\n\n", (int)mesh.positions.size());
  }

  auto bvh = make_triangles_bvh(mesh.triangles, mesh.positions, {});

  // Init bool_state
  auto state = bool_state{};
  if (test.screenspace) {
    state = state_from_screenspace_test(mesh, test, drawing_size, false);
  } else {
    state = state_from_test(mesh, test, drawing_size, false);
  }

  stats.polygons = (int)state.polygons.size();
  for (auto& polygon : state.polygons) stats.added_points += polygon.length;
  stats.control_points = state.points.size();
  stats.control_points += (int)state.isecs_generators.size();

  // Execute triangulation
  auto triangulation_timer = simple_timer{};
  slice_mesh(mesh, state);

  stats.triangulation_ns = elapsed_nanoseconds(triangulation_timer);
  stats.total_ns += stats.triangulation_ns;
  stats.sliced_triangles = (int)mesh.triangulated_faces.size();
  for (auto& [face, triangles] : mesh.triangulated_faces)
    stats.added_triangles += triangles.size();

  // Flood-fill for graph creation
  auto flood_fill_timer = simple_timer{};
  state.cells           = make_mesh_cells(mesh.adjacencies, mesh.borders);
  update_virtual_adjacencies(state.cells, mesh.borders);

  stats.flood_fill_ns = elapsed_nanoseconds(flood_fill_timer);
  stats.total_ns += stats.flood_fill_ns;
  stats.graph_nodes = state.cells.size();
  for (auto& cell : state.cells) stats.graph_edges += cell.adjacency.size();
  stats.graph_edges /= 2;

  // Label propagation
  auto labelling_timer = simple_timer{};
  compute_cell_labels(state);

  stats.labelling_ns = elapsed_nanoseconds(labelling_timer);
  stats.total_ns += stats.labelling_ns;

  compute_shapes(state);

  // Saving output scene
  auto scene = make_scene(mesh, state, test.camera, color_shapes);
  if (output_scene_filename.size()) {
    // for (auto& shape : scene.shapes) {
    //   save_shape(shape, output_scene_filename + "/")
    // }
    save_scene(output_scene_filename, scene, error);
  }

  // Saving render and cell adjacency graph
  save_image(
      output_image_filename, mesh, state, test.camera, color_shapes, spp);
  auto graph_dir      = path_dirname(output_image_filename);
  auto graph_filename = path_basename(output_image_filename) +
                        string("_graph.png");
  auto graph_outfile = path_join(graph_dir, graph_filename);
  save_tree_png(state, graph_outfile.c_str(), "", color_shapes);

  if (color_shapes) {
    auto booleans_timer = simple_timer{};
    for (auto& operation : test.operations) {
      compute_bool_operation(state, operation);
    }
    stats.boolean_ns = elapsed_nanoseconds(booleans_timer);
    stats.total_ns += stats.boolean_ns;
  }

  // output timings and stats:
  // model, model_triangles, genus,
  // triangulation_secs, flood_fill_secs, labelling_secs, boolean_secs,
  // total_secs, polygons, control_points, added_points, sliced_triangles,
  // added_triangles, graph_nodes, graph_edges, graph_cycles
  if (stats_filename.size()) {
    auto stats_file = fopen(stats_filename.c_str(), "a");
    fprintf(stats_file,
        "%s, %d, %d, %f, %f, %f, %f, %f, %d, %d, %d, %d, %d, %d, %d, %d\n",
        stats.model_filename.c_str(), stats.model_triangles, stats.genus,
        stats.triangulation_ns, stats.flood_fill_ns, stats.labelling_ns,
        stats.boolean_ns, stats.total_ns, stats.polygons, stats.control_points,
        stats.added_points, stats.sliced_triangles, stats.added_triangles,
        stats.graph_nodes, stats.graph_edges, stats.graph_cycles);
    fclose(stats_file);
  }
}
