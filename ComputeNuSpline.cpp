
/***************************************************************
 *  Method implementations for the corresponding header file.  *
 *                                                             *
 *  First edition: April 2015                                  *
 ***************************************************************/

#include <iostream>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/multi/geometries/multi_point.hpp>
#include "ComputeNuSpline.hpp"

#define _USE_MATH_DEFINES
#include <math.h>

#define pi M_PI


template <typename SphericalPoint>
boost::geometry::model::point<double, 3, boost::geometry::cs::cartesian>
to_Cartesian(SphericalPoint const& spherical_point)
{
    typedef boost::geometry::model::point
        <
            double, 3, boost::geometry::cs::cartesian
        > cartesian_point_type;

    double lon = boost::geometry::get_as_radian<0>(spherical_point);
    double lat = boost::geometry::get_as_radian<1>(spherical_point);

    double x = cos(lat)*cos(lon);
    double y = cos(lat)*sin(lon);
    double z = sin(lat);

    return cartesian_point_type(x, y, z);
}


template <typename CartesianPoint>
boost::geometry::model::point
<
          double, 2, boost::geometry::cs::spherical_equatorial
          <
              boost::geometry::radian
          >
>
to_spherical_equatorial(CartesianPoint const& cartesian_point)
{
    typedef boost::geometry::model::point
        <
            double, 2, boost::geometry::cs::spherical_equatorial
            <
                boost::geometry::radian
            >
        > spherical_point_type;

    double x = boost::geometry::get<0>(cartesian_point);
    double y = boost::geometry::get<1>(cartesian_point);
    double z = boost::geometry::get<2>(cartesian_point);

    double lat = atan2(z, sqrt(x*x+y*y));
    double lon = atan2(y, x);

    return spherical_point_type(lon, lat);
}

/*
 * Default constructor. Needs perhaps initialization, although
 * currently only the argument-based constructor is used.
 */
ComputeNuSpline::ComputeNuSpline() {
  // TODO: initialization?
}


/*
 * Argument-based constructor. Initializes private variables to given data.
 */
ComputeNuSpline::ComputeNuSpline(mp_type quats, vector<double> nvals, double tolerance, int maxit, bool uniform) {
  this->q = quats;
  this->n = nvals;
  this->tol = tolerance;
  this->MAX_ITERS = maxit;

  cout << "tol = " << this->tol << endl;
  cout << "iters = " << this->MAX_ITERS << endl;

  // Initialize the approximations d_i for the control points
  N = 0;
  for (mp_type::iterator it = q.begin(); it != q.end(); it++) {
    append(d, *it);    
    N++;
  } 

  cout << "N is " << N << endl;

  // Compute knot values and distances.
  t.push_back(0.0);
  h.push_back(0);
  if (!uniform) {
    for (mp_type::iterator i = q.begin() + 1; i != q.end(); i++) {
      h.push_back(boost::geometry::distance(*i, *(i-1)));
    }
  } else {
    double hstep = 1.0/(N-1.0);
    for (mp_type::iterator i = q.begin() + 1; i != q.end(); i++) {
      h.push_back(hstep);
    }
  }
  h.push_back(0);

  for (int i = 1; i < h.size()-1; i++) {
    t.push_back(t[i-1] + h[i]);
  }

}


/*
 * Default desctructor. This would need to free the space taken up by the local variables.
 */
ComputeNuSpline::~ComputeNuSpline() {
  // TODO: delete objects
}


/* 
 *  Product of a scalar and a point (i.e. 2D-vector)
 */
spherical_point ComputeNuSpline::prod(double c, spherical_point p) {
  p.set<0>(p.get<0>() * c);
  p.set<1>(p.get<1>() * c);
  return p;
}


/*
 *  Addition of two points (i.e. 2D-vectors)
 */
spherical_point ComputeNuSpline::add(spherical_point q, spherical_point p) {
  p.set<0>(p.get<0>() + q.get<0>());
  p.set<1>(p.get<1>() + q.get<1>());
  return p;
}


/*
 *  Compute the coefficient alpha_i
 */
double ComputeNuSpline::alpha(int i) {
    assert(i > 0 && i < N);
  return boost::geometry::distance(d[i-1], d[i]);
}


/*
 *  Compute the coefficient beta_i
 */
double ComputeNuSpline::beta(int i) {
  assert(i > 0 && i < N-1);
  //cout << "In beta: i = " << i << ", lambda(i) = " << lambda(i) << ", mu(i) = " << mu(i) << endl;
  return boost::geometry::distance( G(d[i-1], d[i], lambda(i)), G(d[i], d[i+1], mu(i)) );
}


/*
 *  Compute the point lying on the spherical arc between p and q for the given
 *  value of the parameter t where 0 <= t <= 1.
 */
spherical_point ComputeNuSpline::G(spherical_point p, spherical_point q, double t) {
  bool b1 = (t > 0.0);
  bool b2 = (t == 0.0);
  bool b3 = (t < 1.0);
  bool b4 = (t == 1.0);
  assert(!(t < 0.0 || t > 1.0));
  double theta = boost::geometry::distance(p, q);
  if (std::abs(theta) <= std::numeric_limits<double>::epsilon()) {
      return p;
  }
  double den = sin(theta);  

  double lon1 = get<0>(p);
  double lat1 = get<1>(p);
  double lon2 = get<0>(q);
  double lat2 = get<1>(q);
  double A = sin((1.0-t) * theta) / den;
  double B = sin(t * theta) / den;
  double x = A*cos(lat1)*cos(lon1) +  B*cos(lat2)*cos(lon2);
  double y = A*cos(lat1)*sin(lon1) +  B*cos(lat2)*sin(lon2);
  double z = A*sin(lat1)           +  B*sin(lat2);
  double len = x*x+y*y+z*z;
  //  assert( std::abs(len-1) <= 1e-6);
  double lat=atan2(z,sqrt(x*x+y*y));
  double lon=atan2(y,x);
  return spherical_point(lon, lat);

}


/*
 *  Compute the coefficient lambda_i
 */
double ComputeNuSpline::lambda(int i) {
  assert(i > 0 && i < N-1);
  double num, den;
  num = gamma(i-1)*h[i-1] + h[i];
  den = gamma(i-1)*h[i-1] + h[i] + gamma(i)*h[i+1];
  //cout << "In lambda: i = " << i << " and I am returning " << num/den << endl;
  return num / den;
}


/*
 *  Compute the coefficient mu_i
 */
double ComputeNuSpline::mu(int i) {
    assert(i >= 0 && i < N);
  double num, den;
  num = gamma(i)*h[i];
  den = gamma(i)*h[i] + h[i+1] + gamma(i+1)*h[i+2];
  //cout << "In mu, we have gamma(i) = " << gamma(i) << ", gamma(i+1) = " << gamma(i+1) << ", h[i] = " << h[i] <<
  //        ", h[i+1] = " << h[i+1] << " and h[i+2] = " << h[i+2] << endl;
  //cout << "The result for i = " << i << " is therefore " << num/den << endl << endl;
  return num / den;
}


/*
 *  Compute the coefficient delta_i
 */
double ComputeNuSpline::delta(int i) {
    assert(i > 0 && i < N-1);
  double num, den;
  num = h[i];
  den = h[i] + h[i+1];
  return num / den;
}


/*
 *  Compute the coefficient gamma_i
 */
double ComputeNuSpline::gamma(int i) {
    assert(i >= 0 && i < N);
  double num, den;
  num = 2.0*(h[i] + h[i+1]);
  den = n[i]*h[i]*h[i+1] + num;
  return num / den;
}


/*
 *  Compute the next approximation for the control point d_i (eq. 5 in the paper)
 */
spherical_point ComputeNuSpline::F(int i)
{
    namespace bg = boost::geometry;
  
  double den, num1, num2, num3;
  double ai, bi, di, li, mi, ai1;
  
  ai = alpha(i);
  ai1 = alpha(i+1);
  di = delta(i);
  bi = beta(i);
  li = lambda(i);
  mi = mu(i);

  num1 = sin(bi);
  num2 = ( sin( (1.0 - di) * bi ) * sin( (1.0 - li) * ai ) ) / sin( ai );
  num3 = ( sin( di * bi ) * sin( mi * ai1 ) ) / sin( ai1 );
  den = ( sin( (1.0 - di) * bi ) * sin( li * ai ) ) / sin( ai ) + ( sin( di * bi ) * sin( (1.0 - mi) * ai1 ) ) / sin( ai1 );

  typedef bg::model::point<double, 3, bg::cs::cartesian> cartesian_point_type;

  cartesian_point_type qi = to_Cartesian(q[i]);
  cartesian_point_type di_prev = to_Cartesian(d[i-1]);
  cartesian_point_type di_next = to_Cartesian(d[i+1]);

  bg::multiply_value(qi, num1 / den);
  bg::multiply_value(di_prev, -num2 / den);
  bg::multiply_value(di_next, -num3 / den);

  bg::add_point(qi, di_prev);
  bg::add_point(qi, di_next);

  return to_spherical_equatorial(qi);

#if 0
  spherical_point quat1 = prod(num1 / den, q[i]);
  spherical_point quat2 = prod(-1.0 * num2 / den, d[i-1]);
  spherical_point quat3 = prod(-1.0 * num3 / den, d[i+1]);

  spherical_point interp = add(add(quat1, quat2), quat3);
  return interp;
#endif
}


/*
 *  Compute the relative change in approximation (step 4 in the algorithm of the paper)
 */
double ComputeNuSpline::reldiff(mp_type prev)
{
    namespace bg = boost::geometry;
    typedef bg::model::point<double, 3, bg::cs::cartesian> cartesian_point_type;

  double num, den;
  num = den = 0.0;
  for (int i = 0; i < N; i++) {


      cartesian_point_type cur = to_Cartesian(d[i]);
      cartesian_point_type pre = to_Cartesian(prev[i]);

      bg::multiply_value(pre, -1.0);
      bg::add_point(pre, cur); // prev becomes the difference

      num += bg::get<0>(pre) * bg::get<0>(pre) 
          + bg::get<1>(pre) * bg::get<1>(pre) 
          + bg::get<2>(pre) * bg::get<2>(pre);

      den += bg::get<0>(cur) * bg::get<0>(cur) 
          + bg::get<1>(cur) * bg::get<1>(cur) 
          + bg::get<2>(cur) * bg::get<2>(cur);

#if 0
    spherical_point p = add(d[i], prod(-1.0, prev[i]));
    num += (p.get<0>()*p.get<0>() + p.get<1>()*p.get<1>());
    den += (d[i].get<0>()*d[i].get<0>() + d[i].get<1>()*d[i].get<1>());
#endif
  }

  return sqrt(num / den);
}


/* 
 *  Run the approximation algorithm.
 */
bool ComputeNuSpline::Execute() {
  
  int iteration = 0;
  double norm = 10.0;
  mp_type prev = d;

  std::vector<spherical_point> dnew(N);
  dnew[0] = d[0];
  dnew[N-1] = d[N-1];

  while (norm > tol && iteration < MAX_ITERS) {
    for (int i = 1; i < N - 1; i++) { // Note that the iteration involves only the inner nodes!
        //        d[i] = F(i);
        dnew[i] = F(i);
    }
    for (int i = 1; i < N - 1; i++) {
        d[i] = dnew[i];
    }
    iteration++;
    norm = reldiff(prev);
    last_error = norm;
    prev = d;
    std::cout << "At iteration " << iteration << " the error is " << norm << std::endl;
  }
  
  if (norm <= tol)
    return true;
  else
    return false;
}

mp_type ComputeNuSpline::GetNuSpline(int tsteps) {
  double htime = 1.0 / tsteps;

  for (int i = 1; i < N; i++) {
    spherical_point Ri = R(i-1);
    spherical_point Li = L(i); 
    double tt = 0.0;
    for (int j = 0; j <= tsteps; j++) {
      spherical_point p1 = G(q[i-1], Ri, tt);
      spherical_point p2 = G(Ri, Li, tt);
      spherical_point p3 = G(Li, q[i], tt);
      spherical_point z1 = G(p1, p2, tt);
      spherical_point z2 = G(p2, p3, tt);
      spherical_point w = G(z1, z2, tt);

      if (j < tsteps-1) 
        tt += htime;
      else
        tt = 1.0;

      append(sp, w);
    }
  }

  return sp;
}

spherical_point ComputeNuSpline::L(int i) {
    assert(i > 0 && i < N);
  double t = ( gamma(i-1)*h[i-1] + h[i] ) / ( gamma(i-1)*h[i-1] + h[i] + gamma(i)*h[i+1] );
  return G(d[i-1], d[i], t); 
}


spherical_point ComputeNuSpline::R(int i) {
    assert(i >= 0 && i < N);
  double t = ( gamma(i)*h[i] ) / ( gamma(i)*h[i] + h[i+1] + gamma(i+1)*h[i+2] );
  return G(d[i], d[i+1], t); 
}
