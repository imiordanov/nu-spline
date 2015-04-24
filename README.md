# nu-spline

Spherical nu-spline interpolation.

Implementation based on the paper "ν-Quaternion Splines for the Smooth
Interpolation of Orientations" by Gregory M. Nielson, IEEE TRANSACTIONS ON VISUALIZATION AND COMPUTER GRAPHICS, VOL. 10, NO. 2, MARCH/APRIL 2004.

The code is based on the Boost framework, in particular the Geometry package.

Compile with:
```
  g++ driver.cpp
```
and execute with:
```
  ./a.out config
```
Different running configurations are to be implemented via configuration files, refer to the example present.

In subsequent versions the structure of the branch will hopefully improve.




Configuration file guidelines
--------------------------------

The configuration is quite simplistic and allows for some parameterization of the program without requiring recompilation of the source. The switches supported are below:

| Parameter             | Example                 | Description                                                                 |
|-----------------------|-------------------------|-----------------------------------------------------------------------------|
|points_file            |tests/helix15/points     | File containing the points to interpolate (Spherical Equatorial 2D coordinates)
|tensions_file          |tests/helix15/tensions   | File containing the tension values (ν-values in the paper)
|tolerance              |1e-10                    | Error below which the solution is acceptable
|iterations             |10000                    | Max iterations before interrupting exwcution (if convergence is not achieved)
|spline_points          |43                       | Now many points to produce in each interval (more = smoother curve)
|uniform                |0                        | Uniform knot spacing (set to 1) or non-uniform (set to 0), refer to the paper
|plot_knots             |0                        | Whether to plot the control points (set to 1) or not (set to 0) 
|plot_sphere_file       |plot/sphere              | File containing the points of the reference sphere 
|plot_points_file       |plot/points              | Where to write the interpolation points in Cartesian 3D coordinates
|plot_spline_file       |plot/spline              | Where to write the points of the interpolating spline
|plot_control_pts_file  |plot/control_points      | Where to write the control points