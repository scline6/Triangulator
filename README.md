# Triangulator
A header-only library for Constrained Delaunay Triangulation


## Features

1. Triangulate a single contour using ear-clipping method and satifsying the [Fast Industrial-Strength Triangulation of Polygons (FIST)](http://www.cosy.sbg.ac.at/~held/projects/triang/triang.html) robustness guarantees.
2. Triangulate a nested contour by performing nesting & bridging, then calling Triangulator::triangulate1.
3. Optional Constrained Delaunay Triangulation (CDT) via Delaunay flips.
4. Optional Mesh Refinement via Chew's 2nd Algorithm.


## Distiguishing Characteristics

1. Triangulator uses Hilbert spatial hashing to prune distant vertices from the point-in-triangle tests.
2. Triangulator's uses a priority queue to greedily select the ear with the best triangle quality, which reduces the number of Delaunay flip operations.
3. Both ear-clipping (circular list on vector) and Delaunay flips (quad-edge) are carried out with in-place algorithms in contiguous memory for optimal cache performance.
4. Triangulator handles either winding order (CCW or CW) for single contours and returns triangles in winding order that is consistent with the input.


## Usage

1. Coming soon.


## Roadmap

1. Test Triangulator::triangulate1 on challenging polygons with touching points, self-intersections, duplicate points, collinear points, etc.
2. Implement and test Triangulator::triangulate, which will handle nested polygon contours.
3. Implement and test Triangulator::performChewRefinement.


## Personal Goals

1. Build a Tensorflow deep-learning model to perform CDT by training it against Triangulator.  (This is the original reason I started the project.)
