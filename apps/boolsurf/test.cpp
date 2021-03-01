#include <yocto/yocto_commonio.h>
using namespace yocto;

int main(int num_args, const char *args[]) {
  auto test_filename = "data/tests/test.json"s;

  // parse command line
  auto cli = make_cli("test", "test boolsurf algorithms");
  add_option(cli, "input", test_filename, "Input test filename (.json).", true);
  parse_cli(cli, num_args, args);

  // Welcome
  printf("Hello, %s!\n", test_filename.c_str());
}