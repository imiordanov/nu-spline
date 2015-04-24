
/************************************************************
 *  Header file for the nu-spline interpolation algorithm.  *
 *                                                          *
 *  First edition: April 2015                               *
 ************************************************************/

#include <iostream>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/multi/geometries/multi_point.hpp>

using namespace boost::geometry;
using namespace std;

/*
 *  Typedef for spherical equatorial points.
 */
typedef boost::geometry::model::point
    <
      double, 2, boost::geometry::cs::spherical_equatorial<boost::geometry::radian>
    > spherical_point;

/*
 *  Typedef for cartesian points.
 */
typedef boost::geometry::model::point
    <
      double, 3, boost::geometry::cs::cartesian
    > cartesian_point;

/*
 *  Typedef for multi point container (used only with spherical equatorial points)
 */
typedef boost::geometry::model::multi_point<spherical_point> mp_type;

/*
 *  Class definition.
 */
class ComputeNuSpline {
  
  private:
    mp_type q;              // Holds the input points (quaternions)
    mp_type  d;             // Holds the computed control points
    mp_type sp;             // Interpolating spline. Computed on request.
    int N, MAX_ITERS;       // Number of input points, maximum iterations
    double tol;             // Tolerance (equivalent to desired accuracy)
    double last_error;      // Holds the last error calculated (for informative purposes)
    vector<double> t;             // Holds the knot values (denoted in the paper by t_i, hence t)
    vector<double> n;             // Holds the tension values (denoted by nu_i, hence n)
    vector<double> h;             // Holds the spacing, i.e. h[i] = t[i] - t[i-1], i = 1, ..., N
    

    double alpha(int i);        // Compute the coefficient alpha_i (eq. 6 in the paper)
    double beta(int i);         // Compute the coefficient beta_i (eq. 8 in the paper)
    spherical_point F(int idx); // Compute next approximation for the d_i (eq. 5 in the paper)
    double lambda(int i);       // Compute lambda_i (used in the above equations)
    double mu(int i);           // Compute mu_i (used  in the above equations)
    double delta(int i);        // Compute delta_i (used in the above equations)
    double gamma(int i);        // Compute gamma_i (used in the above equations)

    spherical_point L(int i);            // Function used to compute the interpolating spline.
    spherical_point R(int i);            // Again, used to compute the interpolating spline.

    // Compute relative difference (relative change in approximation from the previous step).
    double reldiff(mp_type prev);

    // Due to inability to directly apply scalar-point multiplication and point-point addition,
    // these two functions are employed. Operator overloading perhaps?
    spherical_point prod(double c, spherical_point p);
    spherical_point add(spherical_point q, spherical_point p);
    

    // Find the point lying on the arc between p and q for the parameter 0 <= t <= 1
    // (when t = 0 it returns p, when t = 1 returns q and when t = 0.5 returns the midpoint of the arc)
    spherical_point G(spherical_point p, spherical_point q, double t);

  public:
    // Default constructor (currently unused)
    ComputeNuSpline();

    // Constructor (currently used)
    ComputeNuSpline(mp_type quats, vector<double> nvals, double tol, int maxiters, bool uniform);
    
    // Default destructor (currently not doing anything)
    ~ComputeNuSpline();

    // Set-methods (used to manually set input values)
    void SetQuats(mp_type quats)   {   q = quats;    }
    void SetTvals(vector<double> tvals)  {   t = tvals;    }
    void SetNvals(vector<double> nvals)  {   n = nvals;    }
    
    // Get-methods (used to retrieve input or output values)
    mp_type GetQads()     {   return q;   }
    vector<double> GetTvals()   {   return t;   }
    vector<double> GetNvals()   {   return n;   }
    mp_type GetDvals()    {   return d;   }
    double GetLastError() {   return last_error; }
    mp_type GetNuSpline(int N);

    // Main mathod (invokes approximation)
    bool Execute();

};