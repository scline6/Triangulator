# Triangulator
A header-only library for Constrained Delaunay Triangulation


## Features

1. Triangulate a single contour using ear-clipping method and satifsying the [Fast Industrial-Strength Triangulation of Polygons (FIST)](http://www.cosy.sbg.ac.at/~held/projects/triang/triang.html) robustness guarantees.
2. Triangulate a nested contour by performing nesting & bridging, then calling Triangulator::triangulate1.
3. Optional Constrained Delaunay Triangulation (CDT) via Delaunay-Lawson flips.
4. Optional Mesh Refinement via Chew's 2nd Algorithm.


## Distiguishing Characteristics

1. Triangulator uses Hilbert spatial hashing to prune distant vertices from the point-in-triangle tests and accelerates the ear-clipping algorithm to O(n*log(n)) performance.
2. Triangulator's Ear-clipping algorithm uses a priority queue to greedily select the ear with the best triangle quality.  This minimizes the number of Delaunay-Lawson flips and generally avoids the worst case of O(n^2) for the CDT portion.
3. Both Ear-clipping and Delaunay flips are carried out with in-place algorithms in contiguous memory.  This approach gives optimal cache performance.
4. Triangulator handles either winding order (CCW or CW) for single contours and always returns triangles in the same winding order.  For nested contours, this is almost true, but Triangulator will adjust interior contours as necessary.


## Usage

1. Coming soon.


## Roadmap

1. Test Triangulator::triangulate1 on challenging polygons with touching points, self-intersections, duplicate points, collinear points, etc.
2. Debug issues with Triangulator::performDelaunayFlips.
3. Implement and test Triangulator::triangulate, which will handle nester polygon contours.
4. Implement and test Triangulator::performChewRefinement.


## Personal Goals

1. Build a Tensorflow deep-learning model to perform CDT by training it against Triangulator.  (This is the original reason I started the project.)
