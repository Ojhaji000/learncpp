#include <string>
#include <iostream>
#include <cmath>
#include <list>
#include <deal.II/numerics/data_out.h>
// #include <vector>
using namespace std;
using namespace dealii;
int main()
{
    // problem dimensions
    double b = 9; // in meters
    double h = 0.5; // in meters
    
    // Material parameters
    double k = 1000; // spring constant in Newton/meter
    double EA0 = 5e6; // (modulus of elasticity) x (cross sectional area) in Newton per meter squared

    // External load increment
    double DF = 50; // in Newtons
    
    // Solver parameters
    int N = 30;
    double tol = 1e-6;
    int iterMax = 5;

    // Some useful functions
    auto l = [=](double v)
    { return sqrt(pow(b,2) + pow(h-v,2)); };

    auto F = [=](double v)
    { return (-EA0*(h-v)/l(v) * (l(v)-l(0))/l(0) + k* v); };

    auto dFdv = [=](double v)
    { return (EA0/l(v)) * pow((h-v)/l(v),2) + k + (EA0/l(v))*(l(v)-l(0))/l(0); };
    
    /////////////////////////////
    // newton raphson iteraton //
    /////////////////////////////

    //Initialize
    double v = 0;
    double Dv = 0;
    double Fext = 0;

    double DFext = DF;

    list<list<double>> output;
    // = {{0, 0}};
    // vector<int> rangeVector {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,};
    for (int i = 0; i < N;i++)
    {
        cout << "=================" << endl;
        cout << " load step %i % i" << endl;
        cout << "=================" << endl;
        cout << " NR iter : |Fext-F(v)|" << endl;

        ////////////////////////////////////
        // compute the new external force //
        ////////////////////////////////////

        Fext = Fext + DFext;

        // Initialize Newton-Raphson iteration parameters
        double error = 1;
        double iiter = 0;

        /////////////////////////////////////////
        // compute the derivative of F w.r.t v //
        /////////////////////////////////////////

        while(error > tol)
        {
            //////////////////
            // solve for dv //
            //////////////////

            double dv = (1 / dFdv(v + Dv)) * (Fext - F(v + Dv));

            ////////////////////
            // Update Delta v //
            ////////////////////

            Dv+=dv;


            /////////////////////
            //COnvergence check//
            /////////////////////
            error = abs(Fext - F(v + Dv));

            // Increment the Newton-Raphson iteration counter
            iiter += 1;

            cout << "Iter"<< iiter<< ":"<< error << endl;

            if(iiter == iterMax)
            {
                throw runtime_error("RUNTIME ERROR: Newton-Raphson iterations did not converge!");
            }
        }
        
        // update the displacement
        v += Dv;
        Dv = 0;

        // store the output
        output = {{v, F(v)}};

        cout << "======================" << endl;
    }

    /////////////////////
    // post processing //
    /////////////////////

    dealii::DataOut<dim> data_out;
}