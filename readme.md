# Yocto/GL: Tiny C++ Libraries for Data-Oriented Physically-based Graphics

![windows build badge](https://github.com/xelatihy/yocto-gl/workflows/windows-build/badge.svg)
![macos build badge](https://github.com/xelatihy/yocto-gl/workflows/macos-build/badge.svg)
![ubuntu build badge](https://github.com/xelatihy/yocto-gl/workflows/ubuntu-build/badge.svg)

Yocto/GL is a collection of small C++17 libraries for building
physically-based graphics algorithms released under the MIT license.
Yocto/GL is written in a deliberately data-oriented style for ease of
development and use.
Yocto/GL is split into small libraries to make code navigation easier.
See each header file for documentation.

- `yocto/yocto_math.{h}`: fixed-size vectors, matrices, rigid frames, rays,
  bounding boxes, transforms
- `yocto/yocto_color.{h}`: color conversion, color adjustment, tone mapping
  functions, Perlin noise, shading and integration utilities
- `yocto/yocto_geometry.{h}`: geometry functions, ray-primitive intersection,
  point-primitive overlap
- `yocto/yocto_noise.{h}`: Perlin noise
- `yocto/yocto_sampling.{h}`: random number generation, generation of points
  and directions, Monte Carlo utilities
- `yocto/yocto_shading.{h}`: evaluation and sampling of fresnel functions,
  bsdf lobes, transmittance lobes, phase functions
- `yocto/yocto_shape.{h,cpp}`: various utilities for manipulating
  triangle meshes, quads meshes and line sets, computation of normals and
  tangents, linear and Catmull-Clark subdivision, mesh loading and saving,
  procedural shapes generation, ray intersection and closest point queries
- `yocto/yocto_mesh.{h,cpp}`: computational geometry utilities for triangle
  meshes, mesh geodesic, mesh cutting, mesh loading and saving
- `yocto/yocto_bvh.{h,cpp}`: ray intersection and closest point queries of
  triangle meshes, quads meshes, line sets and instances scenes using a
  two-level bounding volume hierarchy
- `yocto/yocto_image.{h,cpp}`: simple image data type, image resizing,
  tone mapping, color correction, image loading and saving,
  procedural images, procedural sun-sky, advanced color conversion utilities
- `yocto/yocto_sceneio.{h,cpp}`: simple scene and scene loading and
  saving of Ply/Obj/Pbrt/glTF and a custom and scalable Json format
- `yocto/yocto_trace.{h,cpp}`: path tracing of surfaces and hairs supporting
  area and environment illumination, microfacet GGX and subsurface scattering,
  multiple importance sampling
- `yocto/yocto_modelio.{h,cpp}`: parsing and writing for Ply/Obj/Pbrt formats
- `yocto/yocto_commonio.h`: printing utilities, file io utilities,
  command line parsing
- `yocto/yocto_json.h`: JSON data type, json io utilities,
  command line parsing
- `yocto/yocto_common.h`: container and iterator utilities
- `yocto/yocto_parallel.h`: concurrency utilities

You can see Yocto/GL in action in the following applications written to
test the library:

- `apps/yscenetrace.cpp`: command-line path-tracer
- `apps/ysceneitrace.cpp`: interactive path-tracer
- `apps/ysceneitraces.cpp`: simpler version of `apps/ysceneitrace.cpp` for demos
- `apps/ysceneproc.cpp`: command-line scene manipulation and conversion
- `apps/yshapeproc.cpp`: command-line mesh manipulation and conversion
- `apps/yimageview.cpp`: Hdr/Ldr image viewer with tonemapping and color grading
- `apps/yimageviews.cpp`: simpler version of `apps/yimageview.cpp` for demos
- `apps/yimageproc.cpp`: command-line image manipulation
- `apps/yimagedenoise.cpp`: command-line image denoiser that uses Intel Open Image
  Denoise
- `apps/ysceneview.cpp`: simple OpenGL viewer

Here are some test images rendered with the path tracer. More images are
included in the [project site](https://xelatihy.github.io/yocto-gl/).

![Example materials: matte, plastic, metal, glass, subsurface, normal mapping](images/features1.png)

![Example shapes: procedural shapes, Catmull-Clark subdivision, hairs, displacement mapping](images/features2.png)

![Image rendered with Yocto/GL path tracer. Model by Disney Animation Studios.](images/island.png)

## Design Considerations

Yocto/GL follows a "data-oriented programming model" that makes data explicit.
Data is stored in simple structs and accessed with free functions or directly.
All data is public, so we make no attempt at encapsulation.
Most objects is Yocto/GL have value semantic, while large data structures
use reference semnatic with strict ownership. This means that everything
can be trivially serialized and there is no need for memory management.

We do this since this makes Yocto/GL easier to extend and quicker to learn,
with a more explicit data flow that is easier when writing parallel code.
Since Yocto/GL is mainly used for research and teaching,
explicit data is both more hackable and easier to understand.

In terms of code style we prefer a functional approach rather than an
object oriented one, favoring free functions to class methods. All functions
and data are defined in sibling namespaces contained in the `yocto` namespace
so libraries can call all others, but have to do so explicitly.

The use of templates in Yocto was the reason for many refactoring, going
from no template to heavy template use. At this point, Yocto uses some templates
for readability. In the future, we will increase the use of templates in math
code, while keeping many APIs explicitly typed.

We do not use exception for error reporting, but only to report "programmers"
errors. For example, IO operations use boolean flags and error strings for
human readable errors, while exceptions are used when preconditions or
post conditions are violated in functions.

The current version of the library (2.x) is a major refactoring of the previous
library versions (1.x) in three main aspects. First, we now allow the use of
reference semantic via pointers and adopt it for all large objects, while
keeping value semantic for all others. We did this to avoid erroneous copies
that cannot detected and avoided at compile time. Second, we had trouble
interacting with C libraries that mostly use reference semantic. Third, we
reduce the use of exceptions, again for better integration with external code.

## Credits

Main contributors:

- Fabio Pellacini (lead developer): [web](http://pellacini.di.uniroma1.it), [github](https://github.com/xelatihy)
- Edoardo Carra: [github](https://github.com/edoardocarra)
- Giacomo Nazzaro: [github](https://github.com/giacomonazzaro)

This library includes code from the [PCG random number generator](http://www.pcg-random.org),
boost `hash_combine`, and public domain code from `github.com/sgorsten/linalg`,
`gist.github.com/badboy/6267743` and `github.com/nothings/stb_perlin.h`.
Other external libraries are included with their own license.

## Compilation

This library requires a C++17 compiler and is know to compiled on
OsX (Xcode >= 11), Windows (MSVC 2019) and Linux (gcc >= 9, clang >= 9).

You can build the example applications using CMake with
`mkdir build; cd build; cmake ..; cmake --build .`

Yocto/GL depends on `stb_image.h`, `stb_image_write.h`, `stb_image_resize.h` and
`tinyexr.h` for image loading, saving and resizing, `cgltf.h` and `json.hpp`
for glTF and JSON support, and `filesystem.hpp` to support C++17 filesystem API
when missing. All dependencies are included in the distribution.

Yocto/GL optionally supports building OpenGL demos, which are handled by including
glad, GLFW, ImGui as dependencies in apps. OpenGL support might eventually
become part of the Yocto/GL libraries. OpenGL support is enabled by defining
the cmake option `YOCTO_OPENGL` and contained in the `yocto_gui` library.

Yocto/GL optionally supports the use of Intel's Embree for ray casting.
At this point, we rely on prebuilt binaries distributed by Intel.
See the main CMake file for how to link to it. Embree support is enabled by
defining the cmake option `YOCTO_EMBREE`.

Yocto/GL optionally supports the use of Intel's Open Image Denoise for denoising
renders. At this point, we rely on prebuilt binaries distributed by Intel.
See the main CMake file for how to link to it. Open Image Denoise support is enabled by
defining the cmake option `YOCTO_DENOISE`. See `apps/yimagedenoise` for
a demonstration.
