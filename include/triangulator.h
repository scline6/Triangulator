#ifndef TRIANGULATOR_H
#define TRIANGULATOR_H




#include <array>
#include <vector>
#include <deque>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <cmath>
#include <algorithm>
#include <limits>



// FOR DEBUGGING, REMOVE LATER
#include <iostream>
#include <QImage>
#include <QPainter>
#include <QFile>
#include <QFileInfo>




const double      TRI_EPSILON                = 1e-12;
const std::size_t TRI_UNDEFINED_INDEX        = std::numeric_limits<std::size_t>::max();
const std::size_t TRI_SPATIAL_HASH_CUTOFF    = 25;
const std::size_t TRI_SPATIAL_HASH_PHASE_OUT = 10;




class Triangulator
{

public:    // public API functions


    typedef  std::array<double,2>  Vec2;


    typedef  std::array<std::array<double,2>,2>  Mat2x2;


    struct Diagnostics
    {
        double polygonArea;
        double totalTriangleArea;
        double areaDiff;
        int numTriangles;
        int expectedNumTriangles;
        int triDiff;
        int numTrianglesWithPointOnBoundary;
        int numTrianglesWithPointInside;
    };


    std::vector<std::array<std::size_t,3> >
    triangulate1(const std::vector<Vec2>& contour,
                 Diagnostics& diagnostics,
                 const bool& delaunay = true,
                 const bool& includeDegen = true);


    std::vector<std::array<std::array<std::size_t,2>,3> >
    triangulate(const std::vector<std::vector<Vec2> >& contours,
                Diagnostics& diagnostics,
                const bool& delaunay = true,
                const bool& includeDegen = true);


    void performDelaunayFlips(const std::vector<Vec2>& contour,
                              std::vector<std::array<std::size_t,3> >& triangles);


    void performChewRefinement(const std::vector<Vec2>& contour,
                               std::vector<std::array<std::size_t,3> >& triangles);


    int DEBUG__WRITE_IMAGE(const std::vector<Vec2>& contour,
                           std::vector<std::array<std::size_t,3> >& triangles,
                           const QString& imageFilePath);


    int DEBUG__WRITE_IMAGE(const std::vector<std::vector<Vec2> >& contours,
                           std::vector<std::array<std::array<std::size_t,2>,3> >& triangles,
                           const QString& imageFilePath);


    std::vector<std::array<std::size_t,2> >
    DEBUG__TEST_DELAUNAY_PROPERTY(const std::vector<Vec2>& contour,
                                  const std::vector<std::array<std::size_t,3> >& triangles);




protected:


    enum IntersectionResult
    {
        EXTERIOR                = 0,
        NOT_TESTED_YET          = 1,
        BOUNDARY                = 2,
        BOUNDARY_BUT_EDGE_CROSS = 3,
        INTERIOR                = 4,
    };


    // Ear struct is used in circular list for ear-clipping
    struct Ear
    {
        Vec2 pos;
        std::int32_t hash;
        std::size_t originalIndex;
        std::size_t prevEar;
        std::size_t nextEar;
        std::size_t prevHash;
        std::size_t nextHash;
        double signedArea;
        double triQuality;
        IntersectionResult reflex;
        bool operator<(const Ear& other) { return this->hash < other.hash; }
    };


    static bool lessThan(const Vec2& a, const Vec2& b)
    {
        if (a[0] < b[0]) return true;
        if (a[0] > b[0]) return false;
        if (a[1] < b[1]) return true;    // Decide tie-break with y-coordinate
        return false;
    }


    // Similar to Ear struct, but used only for bridging holes
    struct Vertex
    {
        Vec2 pos;
        std::size_t contourIndex;
        std::size_t vertexIndex;
        std::size_t prevVertex;
        std::size_t nextVertex;
        std::size_t prevX;
        std::size_t nextX;
        bool operator<(const Vertex& other) { return lessThan(this->pos, other.pos); }
    };


//    Delicate priority logic - lower numbers are higher priority
//
//                             earArea      Positive Degen/Non-Ref Degenerate Degen/Reflex Negative
//    intersection                             (0)        (1)          (2)         (3)        (4)
//
//    EXTERIOR                 (0)              0          2            4          15         20
//    NOT_TESTED_YET           (1)              1          3            5          16         21
//    BOUNDARY                 (2)              8          6            7          17         22
//    BOUNDARY_BUT_EDGE_CROSS  (3)             11         10            9          18         23
//    INTERIOR                 (4)             14         13           12          19         24


    struct Priority
    {
        IntersectionResult intersection;
        double signedArea;
        IntersectionResult reflex;
        double triQuality;
        Priority(const IntersectionResult& pit, const double& sa, const IntersectionResult& r, const double& tq) : intersection(pit), signedArea(sa), reflex(r), triQuality(tq) {}
        bool operator<(const Priority& other) const
        {
            int thisAreaScore  = (std::abs(this->signedArea) <= TRI_EPSILON ? 2 : (this->signedArea > 0.0 ? 0 : 4));    //0 positive, 2 degen, 4 negative
            int otherAreaScore = (std::abs(other.signedArea) <= TRI_EPSILON ? 2 : (other.signedArea > 0.0 ? 0 : 4));    //0 positive, 2 degen, 4 negative
            thisAreaScore  += (this->reflex == IntersectionResult::EXTERIOR ? 1 : (this->reflex == IntersectionResult::INTERIOR ? -1 : 0));
            otherAreaScore += (other.reflex == IntersectionResult::EXTERIOR ? 1 : (other.reflex == IntersectionResult::INTERIOR ? -1 : 0));
            const std::array<std::array<std::size_t,5>,5> score = {std::array<std::size_t,5>{  0,  2,  4, 15, 20 },
                                                                   std::array<std::size_t,5>{  1,  3,  5, 16, 21 },
                                                                   std::array<std::size_t,5>{  8,  6,  7, 17, 22 },
                                                                   std::array<std::size_t,5>{ 11, 10,  9, 18, 23 },
                                                                   std::array<std::size_t,5>{ 14, 13, 12, 19, 24 }};
            const std::size_t thisCompositeScore  = score[int(this->intersection)][thisAreaScore];
            const std::size_t otherCompositeScore = score[int(other.intersection)][otherAreaScore];
            if (thisCompositeScore < otherCompositeScore) return true;
            if (this->triQuality > other.triQuality) return true;    // All other things being equal, pick the best triangle quality
            return false;
        }
    };


    template <typename TPriority, typename TItem>
    struct PriorityQueue
    {
        std::multimap<TPriority,TItem> order;
        std::unordered_map<TItem,TPriority> lookup;
        void insert(const TItem& item, const TPriority& priority)
        {
            this->order.insert( std::make_pair(priority, item) );
            this->lookup.insert( std::make_pair(item, priority) );
        }
        void remove(const TItem& item)
        {
            auto it = this->lookup.find(item);
            if (it == this->lookup.end()) return;
            const TPriority& priority = it->second;
            auto eqrng = this->order.equal_range(priority);
            for (auto it2 = eqrng.first; it2 != eqrng.second; it2++)
            {
                if (it2->second == item)
                {
                    this->order.erase(it2);
                    this->lookup.erase(it);
                    return;
                }
            }
        }
        void update(const TItem& item, const TPriority& priority)
        {
            this->remove(item);
            this->insert(item, priority);
        }
        std::pair<TItem,TPriority> pop_front()    // Must not be empty
        {
            auto it = this->order.begin();
            const TPriority& priority = it->first;
            const TItem& item = it->second;
            this->order.erase(it);
            this->lookup.erase(item);
            return std::make_pair(item, priority);
        }
        std::size_t size() { return this->order.size(); }
    };


    struct Edge
    {
        std::array<std::size_t,2> vertex;
        std::size_t count;
        std::array<std::size_t,2> adjTriangle;
        std::array<std::size_t,2> adjVertex;
        std::array<std::size_t,4> adjEdge;
        Vec2 cCenterL;
        Vec2 cCenterR;
        double cRadSqL;
        double cRadSqR;
    };


    struct Circumcircle
    {
        Vec2 center;
        double radSq;
    };




protected:    // static functions


    static double distanceSquared(const Vec2& a, const Vec2& b)
    {
        const double dx = b[0] - a[0];
        const double dy = b[1] - a[1];
        return (dx * dx + dy * dy);
    }


    static double pointToEdgeDistanceSquared2D(const Vec2& p, const Vec2& a, const Vec2& b, const double& EPS)
    {
        const Vec2 beta = {b[0] - a[0], b[1] - a[1]};
        const double denom = beta[0] * beta[0] + beta[1] * beta[1];
        if (abs(denom) <= EPS) return distanceSquared(p, a);    // beta=0 means a=b
        const double t = ((p[0] - a[0]) * beta[0] + (p[1] - a[1]) * beta[1]) / denom;
        return distanceSquared(p, {a[0] + beta[0] * t, a[1] + beta[1] * t});
    }


    static IntersectionResult pointInTriangleTest2D(const Vec2& a, const Vec2& p, const Vec2& q, const Vec2& r)
    {
        std::array<double,3> test;
        test[0] = (a[0] - r[0]) * (q[1] - r[1]) - (q[0] - r[0]) * (a[1] - r[1]);
        test[1] = (a[0] - p[0]) * (r[1] - p[1]) - (r[0] - p[0]) * (a[1] - p[1]);
        test[2] = (a[0] - q[0]) * (p[1] - q[1]) - (p[0] - q[0]) * (a[1] - q[1]);
        std::array<int,3> intTest;
        intTest[0] = (test[0] > TRI_EPSILON ? 1 : (test[0] < -TRI_EPSILON ? -1 : 0));
        intTest[1] = (test[1] > TRI_EPSILON ? 1 : (test[1] < -TRI_EPSILON ? -1 : 0));
        intTest[2] = (test[2] > TRI_EPSILON ? 1 : (test[2] < -TRI_EPSILON ? -1 : 0));
        int sum  = intTest[0] + intTest[1] + intTest[2];
        int prod = intTest[0] * intTest[1] * intTest[2];
        if (abs(sum) == 3) return IntersectionResult::INTERIOR;
        if (abs(sum) == 2) return IntersectionResult::BOUNDARY;    // On edge
        if (abs(sum) == 1 && prod == 0) return IntersectionResult::BOUNDARY;    // On vertex
        return IntersectionResult::EXTERIOR;
    }


    static Vec2 solve2x2(const Mat2x2& A, const Vec2& B, bool& singular, const double& EPS)
    {
        const double det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
        singular = (det < EPS);
        if (singular) return Vec2();
        return {(A[1][1] * B[0] - A[0][1] * B[1]) / det, (A[0][0] * B[1] - A[1][0] * B[0]) / det};
    }


    static IntersectionResult edgeEdgeIntersectionTest2D(const Vec2& a, const Vec2& b, const Vec2& p, const Vec2& q, const double& EPS)
    {
        Mat2x2 A = {Vec2{b[0]-a[0], p[0]-q[0]}, Vec2{b[1]-a[1], p[1]-q[1]}};
        Vec2 B = {p[0]-a[0], p[1]-a[1]};
        bool singular;
        Vec2 uv = solve2x2(A, B, singular, EPS);
        if (abs(uv[0]) <= EPS || abs(uv[1]) <= EPS) return IntersectionResult::BOUNDARY;
        if (uv[0] < 0.0 || uv[0] > 1.0 || uv[1] < 0.0 || uv[1] > 1.0) return IntersectionResult::EXTERIOR;
        return IntersectionResult::INTERIOR;
    }


    static IntersectionResult hybridIntersectionTest(const std::size_t& ear, const std::size_t& vertex, const std::vector<Ear>& polygon, const double& EPS)
    {
        const Vec2& p = polygon[polygon[ear].prevEar].pos;
        const Vec2& q = polygon[ear].pos;
        const Vec2& r = polygon[polygon[ear].nextEar].pos;
        const Vec2& b = polygon[vertex].pos;
        const IntersectionResult result = pointInTriangleTest2D(b, p, q, r);
        if (result != IntersectionResult::BOUNDARY) return result;
        const Vec2& a = polygon[polygon[vertex].prevEar].pos;
        const Vec2& c = polygon[polygon[vertex].nextEar].pos;
        const IntersectionResult e0 = edgeEdgeIntersectionTest2D(p, r, b, a, EPS);
        if (e0 != IntersectionResult::EXTERIOR) return IntersectionResult::BOUNDARY_BUT_EDGE_CROSS;
        const IntersectionResult e1 = edgeEdgeIntersectionTest2D(p, r, b, c, EPS);
        if (e1 != IntersectionResult::EXTERIOR) return IntersectionResult::BOUNDARY_BUT_EDGE_CROSS;
        return IntersectionResult::BOUNDARY;
    }


    // Compute signed area (+ is a good non-reflex vertex, - is a reflex vertex, near-0 is hard to say)
    static double triangleSignedArea(const Vec2& p, const Vec2& q, const Vec2& r)
    {
        return 0.5 * ((p[0] * q[1] - p[1] * q[0]) + (q[0] * r[1] - q[1] * r[0]) + (r[0] * p[1] - r[1] * p[0]));
    }


    // Use winding area to determine is a vertex is a reflex vertex or not
    static IntersectionResult isReflexVertex(const std::size_t& ear, const std::vector<Ear>& polygon)
    {
        int winding = 0;
        const Vec2& a = polygon[ear].pos;
        const Vec2& q = polygon[polygon[ear].prevEar].pos;
        const Vec2& r = polygon[polygon[ear].nextEar].pos;
        double inwardness = (r[0] - q[0]) * (a[1] - q[1]) - (r[1] - q[1]) * (a[0] - q[0]);
        if (abs(inwardness) < TRI_EPSILON)
        {
            const double u = abs(r[0] - q[0]) > abs(r[1] - q[1]) ? (a[0] - q[0]) / (r[0] - q[0]) : (a[1] - q[1]) / (r[1] - q[1]);
            if (-TRI_EPSILON < u && u < 1.0 + TRI_EPSILON) return IntersectionResult::BOUNDARY;
        }
        if (q[1] <= a[1] && a[1] < r[1] && inwardness > 0.0) winding++;
        if (r[1] <= a[1] && a[1] < q[1] && inwardness < 0.0) winding--;
        std::size_t vertex = polygon[ear].nextEar;
        while (vertex != ear)
        {
            const Vec2& q = polygon[vertex].pos;
            const Vec2& r = polygon[polygon[vertex].nextEar].pos;
            double inwardness = (r[0] - q[0]) * (a[1] - q[1]) - (r[1] - q[1]) * (a[0] - q[0]);
            if (abs(inwardness) < TRI_EPSILON)
            {
                const double u = abs(r[0] - q[0]) > abs(r[1] - q[1]) ? (a[0] - q[0]) / (r[0] - q[0]) : (a[1] - q[1]) / (r[1] - q[1]);
                if (-TRI_EPSILON < u && u < 1.0 + TRI_EPSILON) return IntersectionResult::BOUNDARY;
                if (q[1] <= a[1] && a[1] < r[1] && inwardness > 0.0) winding++;
                if (r[1] <= a[1] && a[1] < q[1] && inwardness < 0.0) winding--;
            }
            vertex = polygon[vertex].nextEar;
        }
        return (winding == 0 ? IntersectionResult::EXTERIOR : IntersectionResult::INTERIOR);
    }


    // Test for reflex vertex using area primarily, but if abs(area)<EPS, then use winding order test
    static IntersectionResult isReflexVertexFast(const std::size_t& ear, const std::vector<Ear>& polygon, double& signedArea)
    {
        const Vec2& p = polygon[polygon[ear].prevEar].pos;
        const Vec2& q = polygon[ear].pos;
        const Vec2& r = polygon[polygon[ear].nextEar].pos;
        signedArea = triangleSignedArea(p, q, r);
        if (abs(signedArea) <= TRI_EPSILON) return isReflexVertex(ear, polygon);
        if (signedArea > 0.0) return IntersectionResult::EXTERIOR;
        return IntersectionResult::INTERIOR;
    }


    static double triangleRadiusRatio(const Vec2& p, const Vec2& q, const Vec2& r)
    {
        const double area = std::abs(triangleSignedArea(p, q, r));
        const Vec2 pq = {q[0]-p[0], q[1]-p[1]};
        const Vec2 qr = {r[0]-q[0], r[1]-q[1]};
        const Vec2 rp = {p[0]-r[0], p[1]-r[1]};
        const double lenPQ = sqrt(pq[0] * pq[0] + pq[1] * pq[1]);
        const double lenQR = sqrt(qr[0] * qr[0] + qr[1] * qr[1]);
        const double lenRP = sqrt(rp[0] * rp[0] + rp[1] * rp[1]);
        const double prod = lenPQ * lenQR * lenRP;
        return (prod < TRI_EPSILON ? 0.0 : 2.0 * area * area * (lenPQ + lenQR + lenRP) / (prod * prod * prod));
    }


    static void triangleCircumcircle(const Vec2& p, const Vec2& q, const Vec2& r, Vec2& center, double& radSq)
    {
        const std::array<double,3> e = {distanceSquared(q, r), distanceSquared(r, p), distanceSquared(p, q)};
        std::array<double,3> w = {e[0] * (e[1] + e[2] - e[0]), e[1] * (e[2] + e[0] - e[1]), e[2] * (e[0] + e[1] - e[2])};
        const double denom = w[0] + w[1] + w[2];
        w[0] /= denom; w[1] /= denom; w[2] /= denom;
        center = {w[0] * p[0] +  w[1] * q[0] + w[2] * r[0], w[0] * p[1] +  w[1] * q[1] + w[2] * r[1]};
        radSq = distanceSquared(center, p);
    }


    static void hilbertRotate(const std::int32_t& n, const std::int32_t& rx, const std::int32_t& ry, std::int32_t& x, std::int32_t& y) {
        if (ry != 0) return;
        if (rx == 1)
        {
            x = n - 1 - x;
            y = n - 1 - y;
        }
        const int tmp = x;
        x = y;
        y = tmp;
    }


    static std::int32_t hilbertHashInt(const std::int32_t& n, std::int32_t& x, std::int32_t& y) {
        std::int32_t d = 0;
        for (std::int32_t s = n/2; s > 0; s /= 2)
        {
            const std::int32_t rx = (x & s) > 0;
            const std::int32_t ry = (y & s) > 0;
            d += s * s * ((3 * rx) ^ ry);
            hilbertRotate(n, rx, ry, x, y);
        }
        return d;
    }


    static int hilbertHashReal(const Vec2& p, const Vec2& lower, const double& inverseTotalWidth)
    {
        const double inverseCellWidth = 32767.0 * inverseTotalWidth;
        std::int32_t xInt = static_cast<int32_t>(inverseCellWidth * (p[0] - lower[0]));
        std::int32_t yInt = static_cast<int32_t>(inverseCellWidth * (p[1] - lower[1]));
        return hilbertHashInt(32768, xInt, yInt);   //32768 = 2^15, could go up to 2^30 or 2^31, but that is just extra looping
    }

};




#endif // TRIANGULATOR_H



