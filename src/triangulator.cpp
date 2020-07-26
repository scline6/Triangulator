#include "triangulator.h"




void Triangulator::performDelaunayFlips(const std::vector<Triangulator::Vec2>& contour,
                                        std::vector<std::array<std::size_t,3> >& triangles)
{
    // Make a jagged array of all edges with adjacent triangle, vertex, and edge information
    std::vector<std::vector<Triangulator::Edge> > allEdges(contour.size());
    std::size_t totalEdgeCount = 0;
    for (std::size_t i = 0; i < triangles.size(); i++)
    {
        for (std::size_t j = 0; j < 3; j++)
        {
            const std::size_t& v0 = triangles[i][j];
            const std::size_t& v1 = triangles[i][(j+1)%3];
            const std::size_t& v2 = triangles[i][(j+2)%3];
            std::size_t v12min, v12max, m;
            if (v1 < v2)
            {
                v12min = v1;
                v12max = v2;
                m = 0;
            }
            else {
                v12min = v2;
                v12max = v1;
                m = 1;
            }
            bool found = false;
            std::size_t k;
            for (k = 0; k < allEdges[v12min].size(); k++)
            {
                if (allEdges[v12min][k].vertex[1] == v12max)
                {
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                totalEdgeCount++;
                k = allEdges[v12min].size();
                allEdges[v12min].push_back(Edge());
                allEdges[v12min][k].vertex = {v12min, v12max};
                allEdges[v12min][k].count = 0;
            }
            allEdges[v12min][k].adjTriangle[m] = i;
            allEdges[v12min][k].adjVertex[m] = v0;
            allEdges[v12min][k].count++;
        }
    }

    // Transfer jagged array into flat non-boundary-only array for in-place edge swaps (and better cache-locality too)
    std::vector<Triangulator::Edge> nbEdges;
    nbEdges.reserve(totalEdgeCount);
    std::vector<std::size_t> nbLookup(allEdges.size()+1);
    for (std::size_t i = 0; i < allEdges.size(); i++)
    {
        nbLookup[i] = nbEdges.size();
        for (std::size_t j = 0; j < allEdges[i].size(); j++)
        {
            if (allEdges[i][j].count != 2) continue;    // Skip boundary edges
            nbEdges.push_back(allEdges[i][j]);
        }
    }
    nbLookup.back() = nbEdges.size();

    // Pre-compute triangle circumcircles
    std::vector<Triangulator::Circumcircle> circumcircles(triangles.size());
    for (std::size_t i = 0; i < triangles.size(); i++)
    {
        const std::size_t& v0 = triangles[i][0];
        const std::size_t& v1 = triangles[i][1];
        const std::size_t& v2 = triangles[i][2];
        Triangulator::triangleCircumcircle(contour[v0], contour[v1], contour[v2], circumcircles[i].center, circumcircles[i].radSq);
    }

    // Figure out the four adjacent edges for each edge and compute the diametral circle of each edge
    for (std::size_t i = 0; i < nbEdges.size(); i++)
    {
        Triangulator::Edge& edge = nbEdges[i];
        std::array<std::size_t,4> v = {edge.vertex[0], edge.adjVertex[1], edge.vertex[1], edge.adjVertex[0]};    //v0, vR, v1, vL
        for (std::size_t j = 0; j < 4; j++)
        {
            edge.adjEdge[j] = TRI_UNDEFINED_INDEX;
            const std::size_t& va = v[j];
            const std::size_t& vb = v[(j+1)%4];
            const std::size_t vabmin = std::min(va, vb);
            const std::size_t vabmax = std::max(va, vb);
            std::size_t k;
            for (k = nbLookup[vabmin]; k < nbLookup[vabmin+1]; k++)
            {
                if (nbEdges[k].vertex[1] == vabmax)
                {
                    edge.adjEdge[j] = k;
                    break;
                }
            }
        }
        const Triangulator::Vec2& a = contour[edge.vertex[0]];
        const Triangulator::Vec2& b = contour[edge.vertex[1]];
    }

    // Initialize priority queue for Delaunay edge flips
    std::deque<std::size_t> priorityQueue;
    std::unordered_set<std::size_t> prioritySet;
    for (std::size_t i = 0; i < nbEdges.size(); i++)
    {
        priorityQueue.push_back(i);
        prioritySet.insert(i);
    }

    // Start Delaunay edge flipping loop
    while (priorityQueue.size() > 0)
    {
        // Pop the highest priority edge from the queue and check if it needs flipping
        const std::size_t edgeIndex = priorityQueue.front();
        priorityQueue.pop_front();
        prioritySet.erase(edgeIndex);
        Triangulator::Edge& edge = nbEdges[edgeIndex];
        const std::size_t v0 = edge.vertex[0];
        const std::size_t v1 = edge.vertex[1];
        const std::size_t vL = edge.adjVertex[0];
        const std::size_t vR = edge.adjVertex[1];
        const std::size_t& t0 = edge.adjTriangle[0];
        const std::size_t& t1 = edge.adjTriangle[1];
        const bool needsFlip = (distanceSquared(circumcircles[t0].center, contour[vR]) < circumcircles[t0].radSq - TRI_EPSILON) ||
                               (distanceSquared(circumcircles[t1].center, contour[vL]) < circumcircles[t1].radSq - TRI_EPSILON);
        if (!needsFlip) continue;    // If neither adjacent vertex the opposite triangle's circumcircle, then it is already good, so skip it

        // Perform edge flip
        triangles[t0] = {v0, vR, vL};
        triangles[t1] = {v1, vL, vR};
        edge.vertex[0] = vR;
        edge.vertex[1] = vL;
        edge.adjVertex[0] = v0;
        edge.adjVertex[1] = v1;

        // Update triangle circumcircles
        Triangulator::triangleCircumcircle(contour[v0], contour[vR], contour[vL], circumcircles[t0].center, circumcircles[t0].radSq);
        Triangulator::triangleCircumcircle(contour[v0], contour[vL], contour[vR], circumcircles[t1].center, circumcircles[t1].radSq);

        // Update connectivity in adjacent edges (a.k.a. quad-edges)
        if (edge.adjEdge[0] < nbEdges.size())
        {
            Triangulator::Edge& edge0R = nbEdges[edge.adjEdge[0]];
            if (v0 == edge0R.vertex[0])
            {
                edge0R.adjTriangle[0] = t0;
                edge0R.adjVertex[0] = vL;
                edge0R.adjEdge[2] = edgeIndex;
                edge0R.adjEdge[3] = edge.adjEdge[3];
            }
            else {
                edge0R.adjTriangle[1] = t0;
                edge0R.adjVertex[1] = vL;
                edge0R.adjEdge[0] = edgeIndex;
                edge0R.adjEdge[1] = edge.adjEdge[3];
            }
        }
        if (edge.adjEdge[1] < nbEdges.size())
        {
            Triangulator::Edge& edgeR1 = nbEdges[edge.adjEdge[1]];
            if (vR == edgeR1.vertex[0])
            {
                edgeR1.adjVertex[0] = vL;
                edgeR1.adjEdge[2] = edge.adjEdge[2];
                edgeR1.adjEdge[3] = edgeIndex;
            }
            else {
                edgeR1.adjVertex[1] = vL;
                edgeR1.adjEdge[0] = edge.adjEdge[2];
                edgeR1.adjEdge[1] = edgeIndex;
            }
        }
        if (edge.adjEdge[2] < nbEdges.size())
        {
            Triangulator::Edge& edge1L = nbEdges[edge.adjEdge[2]];
            if (v1 == edge1L.vertex[0])
            {
                edge1L.adjTriangle[0] = t1;
                edge1L.adjVertex[0] = vR;
                edge1L.adjEdge[2] = edgeIndex;
                edge1L.adjEdge[3] = edge.adjEdge[1];
            }
            else {
                edge1L.adjTriangle[1] = t1;
                edge1L.adjVertex[1] = vR;
                edge1L.adjEdge[0] = edgeIndex;
                edge1L.adjEdge[1] = edge.adjEdge[1];
            }
        }
        if (edge.adjEdge[3] < nbEdges.size())
        {
            Triangulator::Edge& edgeL0 = nbEdges[edge.adjEdge[3]];
            if (vL == edgeL0.vertex[0])
            {
                edgeL0.adjVertex[0] = vR;
                edgeL0.adjEdge[2] = edge.adjEdge[0];
                edgeL0.adjEdge[3] = edgeIndex;
            }
            else {
                edgeL0.adjVertex[1] = vR;
                edgeL0.adjEdge[0] = edge.adjEdge[0];
                edgeL0.adjEdge[1] = edgeIndex;
            }
        }

        // Rotate the edges by 1/4 since the edge was swapped
        edge.adjEdge = {edge.adjEdge[1], edge.adjEdge[2], edge.adjEdge[3], edge.adjEdge[0]};

        // Push the quad-edge neighbors into the queue
        for (std::size_t i = 0; i < 4; i++)
        {
            if (prioritySet.count(edge.adjEdge[i]) > 0)
            {
                priorityQueue.push_back(edge.adjEdge[i]);
                prioritySet.insert(edge.adjEdge[i]);
            }
        }
    }
}




std::vector<std::array<std::size_t,3> >
Triangulator::triangulate1(const std::vector<Triangulator::Vec2>& contour,
                           Triangulator::Diagnostics& diagnostics,
                           const bool& delaunay,
                           const bool& includeDegen)
{
    // Compute polygon area (which flags CW vs CCW) and the bounding box (which is for the spatial hashing)
    const Triangulator::Vec2& first = contour.front();
    const Triangulator::Vec2& last  = contour.back();
    const bool lastIsFirst = (std::abs(first[0] - last[0]) <= TRI_EPSILON && std::abs(first[1] - last[1]) <= TRI_EPSILON);
    const std::size_t N = lastIsFirst ? contour.size()-1 : contour.size();
    diagnostics.polygonArea = 0.0;
    Triangulator::Vec2 lowerCorner = contour[0];
    Triangulator::Vec2 upperCorner = lowerCorner;
    for (std::size_t i = 0; i < N; i++)
    {
        const Triangulator::Vec2& q = contour[i];
        const Triangulator::Vec2& r = contour[(i + 1) % N];
        diagnostics.polygonArea += (q[0] * r[1] - q[1] * r[0]);
        lowerCorner[0] = std::min(lowerCorner[0], q[0]);
        lowerCorner[1] = std::min(lowerCorner[1], q[1]);
        upperCorner[0] = std::max(upperCorner[0], q[0]);
        upperCorner[1] = std::max(upperCorner[1], q[1]);
    }
    diagnostics.polygonArea *= 0.5;
    lowerCorner[0] -= TRI_EPSILON;
    lowerCorner[1] -= TRI_EPSILON;
    upperCorner[0] += TRI_EPSILON;
    upperCorner[1] += TRI_EPSILON;
    const double inverseTotalWidth = 1.0 / std::max(upperCorner[0] - lowerCorner[0], upperCorner[1] - lowerCorner[1]);

    // Initialize the dynamic polygon as a doubly-linked list sorted in hash-order (good for spatial hashing and cache efficiency)
    bool spatialHashingOn = (N > TRI_SPATIAL_HASH_CUTOFF);
    std::vector<Triangulator::Ear> polygon(N);
    for (std::size_t i = 0; i < N; i++)
    {
        polygon[i].pos = contour[i];
        polygon[i].originalIndex = i;
        polygon[i].hash = spatialHashingOn ? Triangulator::hilbertHashReal(polygon[i].pos, lowerCorner, inverseTotalWidth) : i;
    }
    if (spatialHashingOn) std::sort(polygon.begin(), polygon.end());
    std::vector<std::size_t> windingOrderToHashOrder(N);
    for (std::size_t i = 0; i < N; i++)
    {
        windingOrderToHashOrder[polygon[i].originalIndex] = i;
    }
    for (std::size_t i = 0; i < N; i++)
    {
        if (diagnostics.polygonArea >= 0.0)
        {
            polygon[i].prevEar = windingOrderToHashOrder[(polygon[i].originalIndex + N - 1) % N];
            polygon[i].nextEar = windingOrderToHashOrder[(polygon[i].originalIndex + 1) % N];
        }
        else {
            polygon[i].prevEar = windingOrderToHashOrder[(polygon[i].originalIndex + 1) % N];
            polygon[i].nextEar = windingOrderToHashOrder[(polygon[i].originalIndex + N - 1) % N];
        }

        polygon[i].prevHash = (i == 0   ? TRI_UNDEFINED_INDEX : i - 1);
        polygon[i].nextHash = (i == N-1 ? TRI_UNDEFINED_INDEX : i + 1);
    }
    windingOrderToHashOrder.clear();    // Not needed anymore

    // Compute area and triangle quality measure
    for (std::size_t i = 0; i < N; i++)
    {
        const Vec2& p = contour[polygon[i].prevEar];
        const Vec2& q = contour[i];
        const Vec2& r = contour[polygon[i].nextEar];
        //polygon[i].signedArea = Triangulator::triangleSignedArea(p, q, r);
        polygon[i].reflex = Triangulator::isReflexVertexFast(i, polygon, polygon[i].signedArea);
        polygon[i].triQuality = Triangulator::triangleRadiusRatio(p, q, r);
    }

    // Initialize priority queue with all non-negative area (non-reflex vertices)
    Triangulator::PriorityQueue<Priority, std::size_t> earQueue;
    for (std::size_t i = 0; i < N; i++)
    {
        if (polygon[i].signedArea < -TRI_EPSILON) continue;
        Triangulator::Priority priority(Triangulator::PointInPolygon::NOT_TESTED_YET,
                                        polygon[i].signedArea,
                                        polygon[i].reflex,
                                        polygon[i].triQuality);
        earQueue.insert(i, priority);
    }

    // Start the ear-clipping loop
    std::vector<std::array<std::size_t,3> > triangles;
    diagnostics.totalTriangleArea = 0.0;
    diagnostics.numTrianglesWithPointOnBoundary = 0;
    diagnostics.numTrianglesWithPointInside = 0;
    while (triangles.size() < N - 2 && earQueue.size() > 0)
    {
        // Pick the ear with best quality from the queue
        const auto earPriorityPair = earQueue.pop_front();
        const std::size_t& ear = earPriorityPair.first;
        const Triangulator::PointInPolygon& oldIntersection = earPriorityPair.second.intersection;
        const std::size_t prevEar = polygon[ear].prevEar;
        const std::size_t nextEar = polygon[ear].nextEar;
        const Triangulator::Vec2& o = polygon[polygon[prevEar].prevEar].pos;
        const Triangulator::Vec2& p = polygon[prevEar].pos;
        const Triangulator::Vec2& q = polygon[ear].pos;
        const Triangulator::Vec2& r = polygon[nextEar].pos;
        const Triangulator::Vec2& s = polygon[polygon[nextEar].nextEar].pos;
        const double& A = polygon[ear].signedArea;
        const double& Q = polygon[ear].triQuality;
        const std::size_t prevHash = polygon[ear].prevHash;
        const std::size_t nextHash = polygon[ear].nextHash;

        // Perform intersection tests to determine if the candidate really is an ear
        Triangulator::PointInPolygon newIntersection = oldIntersection;    // initialize value
        if (oldIntersection == Triangulator::PointInPolygon::NOT_TESTED_YET)
        {
            spatialHashingOn = (N - triangles.size() > TRI_SPATIAL_HASH_PHASE_OUT);
            newIntersection = Triangulator::PointInPolygon::EXTERIOR;
            std::int32_t hashLower, hashUpper;
            if (spatialHashingOn)
            {
                const Vec2 lower = {std::min(std::min(p[0], q[0]), r[0]), std::min(std::min(p[1], q[1]), r[1])};
                const Vec2 upper = {std::max(std::max(p[0], q[0]), r[0]), std::max(std::max(p[1], q[1]), r[1])};
                hashLower = Triangulator::hilbertHashReal(lower, lowerCorner, inverseTotalWidth);
                hashUpper = Triangulator::hilbertHashReal(upper, lowerCorner, inverseTotalWidth);
            }
            else {
                hashLower = std::numeric_limits<std::int32_t>::min();
                hashUpper = std::numeric_limits<std::int32_t>::max();
            }
            int vertex = polygon[ear].nextHash;
            int sweepDirection = 1;
            while (true)
            {
                // Test if vertex  is in interior, on boundary, or in exterior of the ear
                if (vertex != prevEar && vertex != nextEar)
                {
                    Triangulator::PointInPolygon result = Triangulator::PointInTriangleTest2D(polygon[vertex].pos, p, q, r);
                    if (result == Triangulator::PointInPolygon::INTERIOR)
                    {
                        newIntersection = result;
                        break;
                    }
                    else if (result == Triangulator::PointInPolygon::BOUNDARY)
                    {
                        newIntersection = result;
                    }
                }

                // Traverse other vertices by hash order in two sweeps (upward then downward)
                vertex = sweepDirection > 0 ? polygon[vertex].nextHash : polygon[vertex].prevHash;
                if (sweepDirection > 0 && (polygon[vertex].hash > hashUpper || vertex == TRI_UNDEFINED_INDEX) )
                {
                    vertex = polygon[ear].prevHash;    // Go back and descend from ear vertex
                    sweepDirection = -1;
                }
                else if (sweepDirection < 0 && (polygon[vertex].hash < hashLower || vertex == TRI_UNDEFINED_INDEX) )
                {
                    break;
                }
            }

            // Push intersecting and touching cases back into the queue, use them as ears only as a last resort
            if (newIntersection != Triangulator::PointInPolygon::EXTERIOR)
            {
                Triangulator::Priority priority(newIntersection,
                                                polygon[ear].signedArea,
                                                polygon[ear].reflex,
                                                polygon[ear].triQuality);
                earQueue.insert(ear, priority);
                continue;
            }
        }

        // It is an ear, so clip it by updating connectivity in the data structures
        diagnostics.totalTriangleArea += std::abs(polygon[ear].signedArea);
        if (newIntersection == Triangulator::PointInPolygon::BOUNDARY)
        {
            diagnostics.numTrianglesWithPointOnBoundary++;
        }
        else if (newIntersection == Triangulator::PointInPolygon::INTERIOR)
        {
            diagnostics.numTrianglesWithPointInside++;
        }
        triangles.push_back({polygon[prevEar].originalIndex, polygon[ear].originalIndex, polygon[nextEar].originalIndex});
        polygon[prevEar].nextEar = nextEar;
        polygon[nextEar].prevEar = prevEar;
        polygon[ear].prevEar = TRI_UNDEFINED_INDEX;
        polygon[ear].nextEar = TRI_UNDEFINED_INDEX;
        if (prevHash != TRI_UNDEFINED_INDEX) polygon[prevHash].nextHash = nextHash;
        if (nextHash != TRI_UNDEFINED_INDEX) polygon[nextHash].prevHash = prevHash;
        polygon[ear].prevHash = TRI_UNDEFINED_INDEX;
        polygon[ear].nextHash = TRI_UNDEFINED_INDEX;
        polygon[prevEar].reflex     = Triangulator::isReflexVertexFast(prevEar, polygon, polygon[prevEar].signedArea);
        polygon[prevEar].triQuality = Triangulator::triangleRadiusRatio(o, p, r);
        polygon[nextEar].reflex     = Triangulator::isReflexVertexFast(nextEar, polygon, polygon[nextEar].signedArea);
        polygon[nextEar].triQuality = Triangulator::triangleRadiusRatio(p, r, s);
        if (polygon[prevEar].signedArea >= -TRI_EPSILON)
        {
            Triangulator::Triangulator::Priority priority(PointInPolygon::NOT_TESTED_YET,
                                                          polygon[prevEar].signedArea,
                                                          polygon[prevEar].reflex,
                                                          polygon[prevEar].triQuality);
            earQueue.update(prevEar, priority);
        }
        if (polygon[nextEar].signedArea >= -TRI_EPSILON)
        {
            Triangulator::Triangulator::Priority priority(PointInPolygon::NOT_TESTED_YET,
                                                          polygon[nextEar].signedArea,
                                                          polygon[nextEar].reflex,
                                                          polygon[nextEar].triQuality);
            earQueue.update(nextEar, priority);
        }
    }

    // If delaunay option is on, then perform Delaunay flips until constrained Delaunay triangulation (CDT) is achieved
    if (delaunay) this->performDelaunayFlips(contour, triangles);

    // If removeDegen option is on, then remove degenerate triangles from the result, this is good for graphical applications
    if (!includeDegen)
    {
        std::vector<std::array<std::size_t,3> > trianglesReduced;
        trianglesReduced.reserve(triangles.size());
        for (std::size_t i = 0; i < triangles.size(); i++)
        {
            const Vec2& p = contour[triangles[i][0]];
            const Vec2& q = contour[triangles[i][1]];
            const Vec2& r = contour[triangles[i][2]];
            const double area = Triangulator::triangleSignedArea(p, q, r);
            if (abs(area) < TRI_EPSILON) continue;
            trianglesReduced.push_back(triangles[i]);
        }
        triangles = std::move(trianglesReduced);
    }

    // Return resulting triangulation
    diagnostics.areaDiff = diagnostics.totalTriangleArea - abs(diagnostics.polygonArea);
    diagnostics.numTriangles = int(triangles.size());
    diagnostics.expectedNumTriangles = int(N) - 2;
    diagnostics.triDiff = diagnostics.numTriangles - diagnostics.expectedNumTriangles;
    return std::move(triangles);
}




std::vector<std::array<std::array<std::size_t,2>,3> >
Triangulator::triangulate(const std::vector<std::vector<Vec2> >& contours,
                          Triangulator::Diagnostics& diagnostics,
                          const bool& delaunay,
                          const bool& includeDegen)
{    
    std::vector<std::array<std::array<std::size_t,2>,3> > triangles;
    if (contours.size() == 0) return triangles;

    // Calculate polygon area of each contour to verify winding order
    std::vector<double> polygonAreas(contours.size());
    for (std::size_t i = 0; i < contours.size(); i++)
    {
        polygonAreas[i] = 0.0;
        const std::size_t& N = contours[i].size();
        for (std::size_t j = 0; j < N; j++)
        {
            const Triangulator::Vec2& q = contours[i][j];
            const Triangulator::Vec2& r = contours[i][(j+1)%N];
            polygonAreas[i] += 0.5 * (q[0] * r[1] - q[1] * r[0]);
        }
    }

    // Initialize bridged contour graph - it starts as a multiple circular lists, and after bridging it is one big circular list
    std::size_t numVerticesInBridgedContour = 2 * contours.size() - 2;
    for (std::size_t i = 0; i < contours.size(); i++)
    {
        numVerticesInBridgedContour += contours[i].size();
    }
    std::vector<Triangulator::Vertex> bridgedContour;
    bridgedContour.reserve(numVerticesInBridgedContour);
    for (std::size_t i = 0; i < contours.size(); i++)
    {
        for (std::size_t j = 0; j < contours[i].size(); j++)
        {
            Vertex vertex;
            vertex.pos = contours[i][j];
            vertex.contourIndex = i;
            vertex.vertexIndex = j;
            bridgedContour.push_back(vertex);
        }
    }
    std::sort(bridgedContour.begin(), bridgedContour.end());    // Sort in order of increasing x
    std::vector<std::vector<std::size_t> > originalOrderToXOrder(contours.size());
    for (std::size_t i = 0; i < contours.size(); i++)
    {
        originalOrderToXOrder[i].resize(contours[i].size());
    }
    for (std::size_t k = 0; k < bridgedContour.size(); k++)
    {
        const std::size_t& i = bridgedContour[k].contourIndex;
        const std::size_t& j = bridgedContour[k].vertexIndex;
        originalOrderToXOrder[i][j] = k;
    }
    for (std::size_t k = 0; k < bridgedContour.size(); k++)
    {
        const std::size_t& i = bridgedContour[k].contourIndex;
        const std::size_t& j = bridgedContour[k].vertexIndex;
        const std::size_t& N = contours[i].size();
        std::size_t prev, next;
        if ( (i == 0 && polygonAreas[i] >= 0.0) || (i > 0 && polygonAreas[i] <= 0.0) )
        {
            bridgedContour[k].prevVertex = originalOrderToXOrder[i][(j+N-1)%N];
            bridgedContour[k].nextVertex = originalOrderToXOrder[i][(j+1)%N];
        }
        else {
            bridgedContour[k].prevVertex = originalOrderToXOrder[i][(j+1)%N];
            bridgedContour[k].nextVertex = originalOrderToXOrder[i][(j+N-1)%N];
        }
        bridgedContour[k].prevX = (k == 0 ? TRI_UNDEFINED_INDEX : k-1);
        bridgedContour[k].nextX = (k == bridgedContour.size()-1 ? TRI_UNDEFINED_INDEX : k+1);
    }

    // Find the leftmost point of each contour
    std::vector<std::size_t> leftmost(contours.size());
    std::size_t i = 0;
    for (std::size_t k = 0; k < bridgedContour.size(); k++)
    {
        if (bridgedContour[k].contourIndex == i)
        {
            leftmost[i] = k;
            i++;
            if (i >= leftmost.size()) break;    // Found the last one, so we're done
        }
    }

    // I use a novel hole-bridging method.  David Eberly's method is the standard.
    // Eberly's method:
    //     1. Look for a polygon edge due left of the hole's leftmost point.
    //     2. Form a triangle with the edge (sort of) and leftmost point.  Look for other vertices inside.
    //     3. The closest vertex (either inside or on the edge) forms the bridge with the leftmost point.
    // The problem is Step #1.  It is an O(N) search.  If you build a tree or grid, it is O(log N) or O(1).
    // But building the tree or grid is O(N) and there usually are not enough holes to justify it.
    // So if we accept O(N) methods, then there is an alternative simpler approach.
    // My method:
    //    1. Look for the closest edge that is wholly or partially to the left of the hole's leftmost point.
    //        a. If the edge is wholly to the left, use point-to-edge distance.
    //        b. If the edge is partially to the left, use point-to-point distance with the left point.
    //    2. Once you have the closest edge, find the closer of its two vertices.  Bridge it with the leftmost point.

    // Create bridges until the contour graph is one big circular list
    for (std::size_t m = 1; m < leftmost.size(); m++)
    {
        // Find a the closest edge in the bridgedContour that is to the left of the hole's leftmost point
        double distSqToClosestEdge = std::numeric_limits<double>::max();
        std::array<std::size_t,2> closestEdge;
        const Triangulator::Vec2& a = bridgedContour[leftmost[m]].pos;
        std::size_t vertex = bridgedContour[leftmost[m]].prevX;
        while(vertex != TRI_UNDEFINED_INDEX)
        {
            const Triangulator::Vec2& b = bridgedContour[vertex].pos;
            //if (pow(a[0] - b[0],2) > distSqToClosestEdge) break;    // Once ax-bx>closest, the points are too far, so stop - Can you do this?
            for (std::size_t n = 0; n < 2; n++)    // Two connected edges on each vertex
            {
                const std::size_t& neighbor = (n == 0 ? bridgedContour[vertex].prevVertex : bridgedContour[vertex].nextVertex);
                const Triangulator::Vec2& c = bridgedContour[neighbor].pos;
                if ( !Triangulator::lessThan(b, c) ) continue;    // If c is to the right of b, then skip it to avoid computing twice
                if ( Triangulator::lessThan(a, c))    // If c is to the right of a, then forget edge bc, just compute distance to vertex b
                {
                    const double distSq = distanceSquared(a, b);
                    if (distSq < distSqToClosestEdge)
                    {
                        distSqToClosestEdge = distSq;
                        closestEdge = {vertex, TRI_UNDEFINED_INDEX};
                    }
                }
                else {
                    const double distSq = pointToEdgeDistanceSquared2D(a, b, c, TRI_EPSILON);
                    if (distSq < distSqToClosestEdge)
                    {
                        distSqToClosestEdge = distSq;
                        closestEdge = {vertex, neighbor};
                    }
                }
            }
            vertex = bridgedContour[vertex].prevX;
        }

        // Find a vertex in the bridgeContour that is to the left of the hole's leftmost point and is visible
        std::size_t unobstructedVertex = closestEdge[0];
        if (closestEdge[1] != TRI_UNDEFINED_INDEX)
        {
            const double dist0 = distanceSquared(a, bridgedContour[closestEdge[0]].pos);
            const double dist1 = distanceSquared(a, bridgedContour[closestEdge[1]].pos);
            if (dist1 < dist0) unobstructedVertex = closestEdge[1];
        }

        // Insert two bridge vertices and make the proper reconnections
        Triangulator::Vertex L = bridgedContour[unobstructedVertex];    // Copy vertex on the big loop
        Triangulator::Vertex R = bridgedContour[leftmost[m]];           // Copy vertex on the hole we are trying to bridge
        bridgedContour.push_back(R);
        bridgedContour.push_back(L);
        const std::size_t n = bridgedContour.size();
        bridgedContour[n-2].prevX = leftmost[m];
        bridgedContour[leftmost[m]].nextX = n-2;
        bridgedContour[n-1].prevX = unobstructedVertex;
        bridgedContour[unobstructedVertex].nextX = n-1;
        bridgedContour[n-2].nextVertex = n-1;
        bridgedContour[n-1].prevVertex = n-2;
        bridgedContour[bridgedContour[unobstructedVertex].nextVertex].prevVertex = n-1;
        bridgedContour[unobstructedVertex].nextVertex = leftmost[m];
        bridgedContour[bridgedContour[leftmost[m]].prevVertex].nextVertex = n-2;
        bridgedContour[leftmost[m]].prevVertex = unobstructedVertex;
    }

    // Convert bridgedContour into a vertex vector and call the triangulate1 routine
    std::vector<Triangulator::Vec2> singleContour;
    std::vector<std::size_t> singleToBridgedIndex;
    singleContour.reserve(bridgedContour.size());
    singleToBridgedIndex.reserve(bridgedContour.size());
    const std::size_t start = 0;
    std::size_t v = start;
    do {
        singleContour.push_back(bridgedContour[v].pos);
        singleToBridgedIndex.push_back(v);
        v = bridgedContour[v].nextVertex;
    } while (v != start);
    const std::vector<std::array<std::size_t,3> > triangles1 = this->triangulate1(singleContour, diagnostics, delaunay, includeDegen);

    // Convert triangle vertices from singleContour index to a original contour/vertex pairs
    triangles.reserve(triangles1.size());
    for (std::size_t i = 0; i < triangles1.size(); i++)
    {
        std::array<std::array<std::size_t,2>,3> tri;
        for (std::size_t j = 0; j < 3; j++)
        {
            const std::size_t& k = singleToBridgedIndex[triangles1[i][j]];
            tri[j] = {bridgedContour[k].contourIndex, bridgedContour[k].vertexIndex};
        }
        triangles.push_back(tri);
    }
    return triangles;
}









void Triangulator::performChewRefinement(const std::vector<Vec2>& contour,
                                         std::vector<std::array<std::size_t,3> >& triangles)
{
    // IMPLEMENT LATER //
}




int Triangulator::DEBUG__WRITE_IMAGE(const std::vector<Vec2>& contour,
                                      std::vector<std::array<std::size_t,3> >& triangles,
                                      const QString& imageFilePath)
{
    // Set some image parameters
    const int IMAGE_WIDTH      = 1000;
    const int LINE_WIDTH       = 3;
    const double GROWTH_FACTOR = 1.01;

    // Compute the bounding box for the polygon and triangulation
    Triangulator::Vec2 lower = contour[0];
    Triangulator::Vec2 upper = contour[0];
    for (std::size_t i = 0; i < contour.size(); i++)
    {
        lower[0] = std::min(lower[0], contour[i][0]);
        lower[1] = std::min(lower[1], contour[i][1]);
        upper[0] = std::max(upper[0], contour[i][0]);
        upper[1] = std::max(upper[1], contour[i][1]);
    }
    const Triangulator::Vec2 center = {0.5 * (lower[0] + upper[0]), 0.5 * (lower[1] + upper[1])};
    const double paddedWidth = GROWTH_FACTOR * std::max( upper[0] - lower[0], upper[1] - lower[1] );
    const double scaleFactor = double(IMAGE_WIDTH) / paddedWidth;
    const double imageHalfWidth = 0.5 * double(IMAGE_WIDTH);

    // Initialize the image and painter
    QImage image(IMAGE_WIDTH, IMAGE_WIDTH, QImage::Format::Format_RGB32);
    image.fill(0);
    QPainter painter;
    painter.begin(&image);
    painter.setBrush( QBrush(Qt::lightGray, Qt::BrushStyle::SolidPattern) );
    painter.setPen( QPen(Qt::white, LINE_WIDTH, Qt::SolidLine) );
    std::array<QColor,12> palette = {Qt::red, Qt::green, Qt::blue,
                                     Qt::cyan, Qt::magenta, Qt::yellow,
                                     Qt::darkRed, Qt::darkGreen, Qt::darkBlue,
                                     Qt::darkCyan, Qt::darkMagenta, Qt::darkYellow};

    // Draw triangle edges
    for (std::size_t i = 0; i < triangles.size(); i++)
    {
        const QColor color = palette[i % 12];
        painter.setBrush( QBrush(color, Qt::BrushStyle::SolidPattern) );
        QPointF points[3];
        for (std::size_t j = 0; j < 3; j++)
        {
            points[j].setX( imageHalfWidth + scaleFactor * (contour[triangles[i][j]][0] - center[0]) );
            points[j].setY( imageHalfWidth - scaleFactor * (contour[triangles[i][j]][1] - center[1]) );
        }
        painter.drawPolygon(points, 3, Qt::WindingFill);
    }

    // Finalize painter and save image
    painter.end();
    int result = image.save(imageFilePath);
    return result;
}




int Triangulator::DEBUG__WRITE_IMAGE(const std::vector<std::vector<Triangulator::Vec2> >& contours,
                                     std::vector<std::array<std::array<std::size_t,2>,3> >& triangles,
                                     const QString& imageFilePath)
{
    std::size_t numVerticesTotal = 0;
    for (std::size_t i = 0; i < contours.size(); i++)
    {
        numVerticesTotal += contours[i].size();
    }

    std::vector<Triangulator::Vec2> singleContour(numVerticesTotal);
    std::vector<std::vector<std::size_t> > pairToFlatIndex(contours.size());
    std::size_t k = 0;
    for (std::size_t i = 0; i < contours.size(); i++)
    {
        pairToFlatIndex[i].resize(contours[i].size());
        for (std::size_t j = 0; j < contours[i].size(); j++)
        {
            singleContour[k] = contours[i][j];
            pairToFlatIndex[i][j] = k;
            k++;
        }
    }

    std::vector<std::array<std::size_t,3> > trianglesFlat(triangles.size());
    for (std::size_t i = 0; i < triangles.size(); i++)
    {
        for (std::size_t j = 0; j < 3; j++)
        {
            trianglesFlat[i][j] = pairToFlatIndex[triangles[i][j][0]][triangles[i][j][1]];
        }
    }

    return this->DEBUG__WRITE_IMAGE(singleContour, trianglesFlat, imageFilePath);
}




std::vector<std::array<std::size_t,2> >
Triangulator::DEBUG__TEST_DELAUNAY_PROPERTY(const std::vector<Vec2>& contour,
                                            const std::vector<std::array<std::size_t,3> >& triangles)
{
    std::vector<std::array<std::size_t,2> > encroachedList;
    for (std::size_t i = 0; i < triangles.size(); i++)
    {
        const Triangulator::Vec2& p = contour[triangles[i][0]];
        const Triangulator::Vec2& q = contour[triangles[i][1]];
        const Triangulator::Vec2& r = contour[triangles[i][2]];
        Triangulator::Vec2 center;
        double radSq;
        Triangulator::triangleCircumcircle(p, q, r, center, radSq);
        for (std::size_t j = 0; j < contour.size(); j++)
        {
            if (j == triangles[i][0] || j == triangles[i][1] || j == triangles[i][2]) continue;
            const double distSq = distanceSquared(center, contour[j]);
            if (distSq < radSq) encroachedList.push_back({i, j});
        }
    }
    return encroachedList;
}



