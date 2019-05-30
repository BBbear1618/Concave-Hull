# Concave-Hull

Require:
  numpy, scipy, plotly

plotly is only for visualization of Delaunay, Alpha shape and Concave hull, and can be commented if you want.
Before using plotly, you have to register first. Go to https://plot.ly/python/getting-started/ for more info.

Usage:
  import concave_hull
  ch = concave_hull(pts, alpha)
  
Input:
  pts: a set of 2D points in an array
  alpha: threshold for the generation of Alpha shape

Output:
  a list of rings, each ring is a list of point indices in pts. The ring is closed, i.e. the last element of a ring is the same as the first element of it.
