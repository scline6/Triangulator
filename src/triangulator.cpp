#include "triangulator.h"




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
        polygon[prevHash].nextHash = nextHash;
        polygon[nextHash].prevHash = prevHash;
        polygon[ear].prevHash = TRI_UNDEFINED_INDEX;
        polygon[ear].nextHash = TRI_UNDEFINED_INDEX;
        //polygon[prevEar].signedArea = Triangulator::triangleSignedArea(o, p, r);
        polygon[prevEar].reflex     = Triangulator::isReflexVertexFast(prevEar, polygon, polygon[prevEar].signedArea);
        polygon[prevEar].triQuality = Triangulator::triangleRadiusRatio(o, p, r);
        //polygon[nextEar].signedArea = Triangulator::triangleSignedArea(p, r, s);
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




std::vector<std::array<std::size_t,3> >
Triangulator::triangulate(const std::vector<std::vector<Vec2> >& contours,
                          Triangulator::Diagnostics& diagnostics,
                          const bool& delaunay,
                          const bool& includeDegen)
{
    // PERFORM NESTING //
    // this->triangulate1(bridgedContour, diagnostics, delaunay, includeDegen); //
    // PERFORM INDEX MAPPING
}




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
    std::vector<std::size_t> nbLookup(allEdges.size());
    for (std::size_t i = 0; i < allEdges.size(); i++)
    {
        nbLookup[i] = nbEdges.size();
        for (std::size_t j = 0; j < allEdges[i].size(); j++)
        {
            if (allEdges[i][j].count != 2) continue;    // Skip boundary edges
            nbEdges.push_back(allEdges[i][j]);
        }
    }

    // Figure out the four adjacent edges for each edge and compute the diametral circle of each edge
    for (std::size_t i = 0; i < nbEdges.size(); i++)
    {
        Triangulator::Edge& edge = nbEdges[i];
        std::array<std::size_t,4> v = {edge.vertex[0], edge.adjVertex[1], edge.vertex[1], edge.adjVertex[0]};    //v0, vR, v1, vL
        for (std::size_t j = 0; j < 4; j++)
        {
            const std::size_t& va = v[j];
            const std::size_t& vb = v[(j+1)%4];
            const std::size_t vabmin = std::min(va, vb);
            const std::size_t vabmax = std::max(va, vb);
            std::size_t k;
            for (k = nbLookup[vabmin]; k < nbEdges.size(); k++)
            {
                if (nbEdges[k].vertex[1] == vabmax) break;
            }
            edge.adjEdge[j] = k;
        }
        const Triangulator::Vec2& a = contour[edge.vertex[0]];
        const Triangulator::Vec2& b = contour[edge.vertex[1]];
        edge.midPt = {0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1])};
        edge.radSq = distanceSquared(edge.midPt, a);
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
        const bool needsFlip = (distanceSquared(edge.midPt, contour[vL]) < edge.radSq - TRI_EPSILON) ||
                               (distanceSquared(edge.midPt, contour[vR]) < edge.radSq - TRI_EPSILON);
        if (!needsFlip) continue;    // If neither adjacent vertex is inside the edge's diametral circle, then don't flip it

        // Perform edge flip
        const std::size_t& t0 = edge.adjTriangle[0];
        const std::size_t& t1 = edge.adjTriangle[1];
        triangles[t0] = {v0, vR, vL};
        triangles[t1] = {v1, vL, vR};
        edge.vertex[0] = vR;
        edge.vertex[1] = vL;
        edge.adjVertex[0] = v0;
        edge.adjVertex[1] = v1;
        const Vec2& a = contour[edge.vertex[0]];
        const Vec2& b = contour[edge.vertex[1]];
        edge.midPt = {0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1])};
        edge.radSq = distanceSquared(edge.midPt, a);

        // Update connectivity in adjacent edges
        Triangulator::Edge& edge0R = nbEdges[edge.adjEdge[0]];
        Triangulator::Edge& edgeR1 = nbEdges[edge.adjEdge[1]];
        Triangulator::Edge& edge1L = nbEdges[edge.adjEdge[2]];
        Triangulator::Edge& edgeL0 = nbEdges[edge.adjEdge[3]];
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
        edge.adjEdge = {edge.adjEdge[1], edge.adjEdge[2], edge.adjEdge[3], edge.adjEdge[0]};
        for (std::size_t i = 0; i < 4; i++)
        {
            if (prioritySet.count(edge.adjEdge[i]) > 0) priorityQueue.push_back(edge.adjEdge[i]);
        }
    }
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

