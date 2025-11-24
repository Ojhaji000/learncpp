#include <deal.ii/base/logstream.h>
#include <deal.ii/lac/vector.h>

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <iomanip>

using namespace dealii;

// Problem parameters (see Python file comments)
const double b   = 10.0; // m
const double h   = 0.5;  // m
const double k   = 1000.0; // N/m
const double EA0 = 5.0e6;  // N
const double DF  = 50.0; // N (External load increment)

// Solver parameters
const int    N       = 30;
const double tol     = 1.0e-6;
const int    iterMax = 5;

// =================================================================
// 1. Governing Equations (Scalar Functions)
// =================================================================

/**
 * Calculates the length of the truss member, l(v).
 */
double truss_length(const double v)
{
  return std::sqrt(b * b + (h - v) * (h - v));
}

/**
 * Calculates the total internal force F(v) (Truss + Spring).
 * This is the *negative* of the residual when set equal to Fext.
 * F(v) = -EA0 * (h-v)/l(v) * (l(v)-l(0))/l(0) + k * v
 */
double internal_force(const double v)
{
  const double l_v  = truss_length(v);
  const double l_0  = truss_length(0.0);
  const double delta_l = l_v - l_0;

  // Axial force component (Truss)
  const double F_axial = -EA0 * (h - v) / l_v * (delta_l / l_0);

  // Spring force component
  const double F_spring = k * v;

  return F_axial + F_spring;
}

/**
 * Calculates the tangent stiffness dF/dv.
 */
double tangent_stiffness(const double v)
{
  const double l_v = truss_length(v);
  const double l_0 = truss_length(0.0);
  const double d_l_d_v = (v - h) / l_v; // dl/dv

  // Term 1: Derivative of -EA0 * (h-v)/l(v) * (l(v)-l(0))/l(0)
  // This uses the product rule and chain rule (omitted intermediate steps)
  const double K_truss = (EA0 / l_v) * ((h - v) / l_v) * ((h - v) / l_v)
                       + (EA0 / l_v) * (truss_length(v) - l_0) / l_0
                       + (EA0 / l_v) * (h - v) / l_0 * d_l_d_v;

  // Term 2: Derivative of k * v
  const double K_spring = k;

  return K_truss + K_spring;
}

// =================================================================
// 2. Main Solver Class (mimicking a deal.ii driver structure)
// =================================================================

class ShallowTrussSolver
{
public:
  ShallowTrussSolver() : displacement(1), output_data(1, {0.0, 0.0}) {}
  void solve();

private:
  void post_process() const;

  // Displacement vector (size 1 for single DoF)
  Vector<double> displacement;

  // Store pairs of {displacement, force} for plotting
  std::vector<std::array<double, 2>> output_data;
};

void ShallowTrussSolver::solve()
{
  // Initialize displacement v, increment Dv, and force
  displacement[0] = 0.0; // v
  Vector<double> displacement_increment(1);
  displacement_increment[0] = 0.0; // Dv

  double Fext = 0.0;

  std::cout << std::fixed << std::setprecision(10);

  // Outer loop: Load step iterator
  for (int i = 0; i < N; ++i)
  {
    std::cout << "\n=================================" << std::endl;
    std::cout << " Load step " << i << std::endl;
    std::cout << "=================================" << std::endl;
    std::cout << " NR iter : |Fext - F(v)|" << std::endl;

    // Compute the new external force
    Fext += DF;

    double current_v = displacement[0];
    displacement_increment[0] = 0.0; // Reset Dv for new load step

    // Inner loop: Newton-Raphson iterations
    double error = 1.0;
    int iiter = 0;

    while (error > tol)
    {
      // 1. Current displacement guess
      const double v_guess = current_v + displacement_increment[0];

      // 2. Compute Residual R = F_int(v_guess) - F_ext
      // Note: The Python code solves for F_ext - F_int(v), so we use that sign
      const double residual = Fext - internal_force(v_guess);

      // 3. Compute Tangent Stiffness K = dF_int/dv
      const double tangent_K = tangent_stiffness(v_guess);

      // 4. Solve for dv (in 1D: dv = R / K)
      const double dv = residual / tangent_K;

      // 5. Update Delta v
      displacement_increment[0] += dv;

      // 6. Convergence check
      error = std::abs(Fext - internal_force(current_v + displacement_increment[0]));
      
      iiter++;
      std::cout << "  Iter " << iiter << " : " << error << std::endl;

      if (iiter == iterMax)
      {
        throw std::runtime_error("Newton-Raphson iterations did not converge!");
      }
    }

    // Update the displacement
    displacement[0] += displacement_increment[0];

    // Store the output: { v, F(v) }
    output_data.push_back({ displacement[0], internal_force(displacement[0]) });

    std::cout << "=================================" << std::endl;
  }

  // Post-processing: Write data to file
  post_process();
}

/**
 * Writes the load-displacement curve to a file for external plotting.
 */
void ShallowTrussSolver::post_process() const
{
  std::ofstream output_file("load_displacement.dat");
  output_file << "# v [m]  F [N]" << std::endl;
  output_file << std::scientific << std::setprecision(15);

  for (const auto& point : output_data)
  {
    output_file << point[0] << " " << point[1] << std::endl;
  }
  std::cout << "\nPost-processing data written to load_displacement.dat" << std::endl;
}

// =================================================================
// 3. Main Function
// =================================================================

int main()
{
  try
  {
    // deal.ii logs are disabled by default; this sets up the log stream
    deallog.depth_console(0);

    ShallowTrussSolver solver;
    solver.solve();
  }
  catch (std::exception& exc)
  {
    std::cerr << std::endl
              << "----------------------------------------------------" << std::endl;
    std::cerr << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------" << std::endl;
    return 1;
  }
  catch (...)
  {
    std::cerr << std::endl
              << "----------------------------------------------------" << std::endl;
    std::cerr << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------" << std::endl;
    return 1;
  }
  return 0;
}