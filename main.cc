#include <deal.II/base/point.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <iostream>

using namespace dealii;

void first_grid()
{
  Triangulation<2> square;
  GridGenerator::hyper_rectangle(square, Point<2> (0,0), Point<2>(7,2));

  Triangulation<2> circle;
  GridGenerator::hyper_ball(circle);

  Triangulation<2> combined;
  GridGenerator::merge_triangulations(square, circle, combined);

  // combined.refine_global(4);
  circle.refine_global(2);
  std::ofstream out("shape.svg");
  GridOut grid_out;
  grid_out.write_svg(combined, out);
  // grid_out.write_svg(circle, out);
  std::cout << "Grid written to shape.svg" << std::endl;
}
int main()
{
    first_grid();
    return 0;
}
