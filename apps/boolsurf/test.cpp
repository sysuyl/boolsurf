#include <yocto/yocto_commonio.h>

#include "boolsurf_utils.h"
#include "io.h"
using namespace yocto;

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
  printf("triangles: %d\n", (int)mesh.triangles.size());
  printf("positions: %d\n", (int)mesh.positions.size());

  // Welcome
  printf("Hello, %s!\n", test_filename.c_str());
}