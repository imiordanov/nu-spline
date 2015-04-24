lim = 1;
set xrange [-lim : lim]
set yrange [-lim : lim]
set zrange [-lim : lim]

set view equal
set size ratio -1

if (!exists("knots"))
  splot sphere with dots,   \
      points with points, \
      spline with lines