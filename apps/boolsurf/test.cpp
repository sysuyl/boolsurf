#include <yocto/yocto_commonio.h>
#include <yocto/yocto_sampling.h>

#include "boolsurf.h"
#include "boolsurf_io.h"
using namespace yocto;

vector<mesh_point> sample_points(const bool_mesh& mesh, const shape_bvh& bvh,
    const vec3f& camera_from, const bbox3f& bbox, const vec3f& camera_to,
    float camera_lens, float camera_aspect, uint64_t trial, int num_points = 4,
    int ray_trials = 10000) {
  // init data
  auto points  = vector<mesh_point>{};
  auto rng_ray = make_rng(9867198237913, trial * 2 + 1);
  // try to pick in the camera
  auto  ray_trial = 0;
  float aspect    = size(bbox).x / size(bbox).y;
  auto  uvs       = vector<vec2f>{};
  while (points.size() < num_points) {
    if (ray_trial++ >= ray_trials) break;
    auto uv = rand2f(rng_ray);
    if (points.size()) {
      uv = 2 * uv - vec2f{1, 1};
      uv *= min(aspect, 1 / aspect);
      uv += 2 * uvs[0] - vec2f{1, 1};
      uv = uv * 0.5 + vec2f{0.5, 0.5};
    }

    auto ray  = camera_ray(lookat_frame(camera_from, camera_to, {0, 1, 0}),
        camera_lens, camera_aspect, 0.036f, uv);
    auto isec = intersect_triangles_bvh(
        bvh, mesh.triangles, mesh.positions, ray);
    if (!isec.hit) continue;
    if (isec.element < 0 || isec.element > mesh.triangles.size()) continue;
    if (!(isec.uv.x >= 0 && isec.uv.x <= 1)) continue;
    if (!(isec.uv.y >= 0 && isec.uv.y <= 1)) continue;
    points.push_back({isec.element, isec.uv});
    uvs.push_back(uv);
  }
  // pick based on area
  auto rng_area = make_rng(9867198237913);
  auto cdf      = sample_triangles_cdf(mesh.triangles, mesh.positions);
  while (points.size() < num_points) {
    auto [triangle, uv] = sample_triangles(
        cdf, rand1f(rng_area), rand2f(rng_area));
    points.push_back({mesh_point{triangle, uv}});
  }
  return points;
}

int main(int num_args, const char* args[]) {
  auto test_filename = "data/tests/test.json"s;

  // parse command line
  auto cli = make_cli("test", "test boolsurf algorithms");
  add_option(cli, "input", test_filename, "Input test filename (.json).", true);
  parse_cli(cli, num_args, args);

  auto test = bool_test{};
  if (!load_test(test, test_filename)) {
    print_fatal("Error loading test " + test_filename);
  }

  auto          mesh = bool_mesh{};
  vector<vec2f> texcoords;
  vector<vec3f> colors;
  string        error;

  if (!load_mesh(test.model, mesh.triangles, mesh.positions, mesh.normals,
          texcoords, colors, error)) {
    printf("%s\n", error.c_str());
    print_fatal("Error loading model " + test.model);
  }
  init_mesh(mesh);
  printf("triangles: %d\n", (int)mesh.triangles.size());
  printf("adjacencies: %d\n", (int)mesh.adjacencies.size());
  printf("positions: %d\n", (int)mesh.positions.size());

  auto state = state_from_test(mesh, test);
  compute_cells(mesh, state);

  for (auto& operation : test.operations) {
    compute_bool_operation(state, operation);
  }

  // Welcome
  printf("Cells: %d\n", (int)state.cells.size());
}