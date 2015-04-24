/******************************************************
 * Driver file for the test nu-spline interpolation   *
 * algorithm, based on the work of G. Nielson.        *
 *                                                    *
 * First edition of the code: April 2015              *
 ******************************************************/

#include <iostream>
#include <fstream>
#include <string> 
#include <vector>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/multi/geometries/multi_point.hpp>
#include "ComputeNuSpline.cpp"

using namespace std;

void load_options(char* filename, map<string, string> * options);
void load_points(string filename, mp_type *pts);
void load_tensions(string filename, vector<double> *tensions);

/*
 * The program expects exactly 1 parameter indicating the file containing IO configuration.
 */
int main(int argc, char** argv) {

  /*
   * Check for enough input arguments.
   */
  if (argc < 2) {
    cout << "Input data needs to be provided from external files!" << endl;
    cout << "Please use as indicated: " << argv[0] << " config_file" << endl;
    return -1;
  }

  cout << "Reading data..." << endl;
  map<string,string> options;
  load_options(argv[1], &options);

  cout << "Points file: " << options["points_file"] << endl;
  cout << "Tensions file: " << options["tensions_file"] << endl;
  cout << "Tolerance: " << options["tolerance"] << endl;
  cout << "Iterations: " << options["iterations"] << endl;

  ofstream ofile;
  mp_type points;
  vector<double> tensions;
  load_points(options["points_file"], &points);
  load_tensions(options["tensions_file"], &tensions);

  /* 
   *  Display data (for debugging purposes). In the loop below the points read in spherical cartesian form
   *  are written in cartesian form into the file "points" so they can be later plotted by gnuplot.
   */
  cout << "Data read successfully!" << endl << endl;
  cout << "Points: " << endl;
  ofile.open("points");
  for (mp_type::iterator it = points.begin(); it != points.end(); it++) {
    cartesian_point p;
    transform(*it, p);
    cout << "   " << dsv(*it) << endl;
    ofile << p.get<0>() << " " << p.get<1>() << " " << p.get<2>() << endl;
  }
  ofile.close();

  /* 
   * Display tension values.
   */
  cout << endl << "Tensions:" << endl;
  for (vector<double>::iterator it = tensions.begin(); it != tensions.end(); it++) {
    cout << *it << "  ";
  }
  cout << endl << endl;

  /*
   *  Create an instance here (initialization allows to retrieve knot values)
   */
  ComputeNuSpline cmp = ComputeNuSpline(points, tensions, atof(options["tolerance"].c_str()), atoi(options["iterations"].c_str()), atoi(options["uniform"].c_str()) == 1);

  /*
   *  Display knot values
   */
  vector<double> t = cmp.GetTvals();
  cout << endl << "Knots:" << endl;
  for (vector<double>::iterator it = t.begin(); it != t.end(); it++) {
    cout << *it << "  ";
  }
  cout << endl;
  

  cout << "Approximation in course..." << endl;

  /* 
   *  Execute the main method of the class.
   *  The value returned indicates whether the approximation has converged in the steps allowed.
   */
  if( cmp.Execute() ) {
    cout << endl << "Approximation converged with requested accuracy!" << endl;
  } else {
    cout << endl << "Approximation reached maximum number of iterations!" << endl;
    cout << "Last error is " << cmp.GetLastError() << endl;
  }

  /*
   * Output the control nodes. The output is also written in the file "approx"
   * in order to allow plotting by gnuplot.
   */
  cout << endl << "The approximated nodes (in cartesian coordinates) are:" << endl;
  ofile.open("knots");
  mp_type d = cmp.GetDvals();
  for (mp_type::iterator it = d.begin(); it != d.end(); it++) {
    cartesian_point p;
    transform(*it, p);
    cout << dsv(p) << endl;
    ofile << p.get<0>() << " " << p.get<1>() << " " << p.get<2>() << endl;
  }
  ofile.close();

  mp_type sp = cmp.GetNuSpline(atoi(options["spline_points"].c_str()));
  ofile.open("spline");
  for (mp_type::iterator it = sp.begin(); it != sp.end(); it++) {
    // cout << dsv(*it) << endl;
    cartesian_point p;
    transform(*it, p);
    ofile << p.get<0>() << " " << p.get<1>() << " " << p.get<2>() << endl;
  }
  ofile.close();

  cout << endl << "Plotting data... ";
  if (atoi(options["plot_knots"].c_str()) == 1)
    system("gnuplot -persist -e \"sphere='sphere'; points='points'; knots='knots'; spline='spline'\" plotall_knots.gnu");
  else
    system("gnuplot -persist -e \"sphere='sphere'; points='points'; spline='spline'\" plotall.gnu");
  cout << " DONE! Exiting." << endl;

  return 0;
}




void load_options(char* filename, map<string, string> * options) {
  ifstream f;
  string key, value;
  f.open(filename);
  while (f >> key >> value) {
    (*options)[key] = value;
  }
  f.close();

  return;
}

void load_points(string filename, mp_type *pts) {
  ifstream f;
  f.open(filename.c_str());
  double phi, theta;
  while (f >> phi >> theta) {
    append(*pts, spherical_point(phi, theta));
  }
  f.close();

  return;
}


void load_tensions(string filename, vector<double> *tensions) {
  ifstream f;
  f.open(filename.c_str());
  double t;
  while(f >> t) {
    tensions->push_back(t);
  }
  f.close();

  return;
}