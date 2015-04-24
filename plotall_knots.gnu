lim = 1;
set xrange [-lim : lim]
set yrange [-lim : lim]
set zrange [-lim : lim]

set view equal
set size ratio -1

splot sphere with dots,   \
      points with points, \
      spline with lines, \
      knots with linespoints