# Triangulator
A header-only library for Constrained Delaunay Triangulation


## Features

1. Triangulate a single contour using ear-clipping method and satifsying the [Fast Industrial-Strength Triangulation of Polygons (FIST)](http://www.cosy.sbg.ac.at/~held/projects/triang/triang.html) robustness guarantees.
2. Triangulate a contour with hold by bridging the holes and making a single contour.
3. Provide an option to retain or discard degenerate triangles.  Geometry devs will retain and graphics devs will discard.
3. Provide an option for Constrained Delaunay Triangulation (CDT) via Delaunay flips.
4. Provide an option to refine the triangle mesh using Chew's 2nd Algorithm.  This is for finite element folks.


## Distiguishing Characteristics

1. Triangulator uses Hilbert spatial hashing to prune distant vertices from the point-in-triangle tests.
2. Triangulator's priority queue greedily selects the best triangle quality, which reduces Delaunay flips.
3. Both ear-clipping (circular list as vector) and Delaunay flips (quad-edge) use contiguous memory for good cache performance.
4. Triangulator handles CW or CCW and resulting triangles are consistent with input (for single contour).


## Usage

1. Coming soon.


## Roadmap

1. Test Triangulator::triangulate1 on challenging polygons with touching points, self-intersections, duplicate points, collinear points, etc.
2. Implement and test Triangulator::triangulate, which will handle nested polygon contours.
3. Implement and test Triangulator::performChewRefinement.


## Personal Goals

1. Build a Tensorflow deep-learning model to perform CDT by training it against Triangulator.  This is the original reason for the project.
