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


    std::vector<std::array<std::size_t,3> >
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




protected:


    enum PointInPolygon
    {
        EXTERIOR       = 0,
        NOT_TESTED_YET = 1,
        BOUNDARY       = 2,
        INTERIOR       = 3,
    };


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
        PointInPolygon reflex;
        bool operator<(const Ear& other) { return this->hash < other.hash; }
    };


//    Delicate priority logic - lower numbers are higher priority
//
//                    earArea      Positive Degen/Non-Ref Degenerate Degen/Reflex Negative
//    intersection                    (0)        (1)          (2)         (3)        (4)
//
//    EXTERIOR        (0)              0          2            4          12         16
//    NOT_TESTED_YET  (1)              1          3            5          13         17
//    BOUNDARY        (2)              8          6            7          14         18
//    INTERIOR        (3)             11         10            9          15         19


    struct Priority
    {
        PointInPolygon intersection;
        double signedArea;
        PointInPolygon reflex;
        double triQuality;
        Priority(const PointInPolygon& pit, const double& sa, const PointInPolygon& r, const double& tq) : intersection(pit), signedArea(sa), reflex(r), triQuality(tq) {}
        bool operator<(const Priority& other) const
        {
            int thisAreaScore  = (std::abs(this->signedArea) <= TRI_EPSILON ? 2 : (this->signedArea > 0.0 ? 0 : 4));    //0 positive, 2 degen, 4 negative
            int otherAreaScore = (std::abs(other.signedArea) <= TRI_EPSILON ? 2 : (other.signedArea > 0.0 ? 0 : 4));    //0 positive, 2 degen, 4 negative
            thisAreaScore  += (this->reflex == PointInPolygon::EXTERIOR ? 1 : (this->reflex == PointInPolygon::INTERIOR ? -1 : 0));
            otherAreaScore += (other.reflex == PointInPolygon::EXTERIOR ? 1 : (other.reflex == PointInPolygon::INTERIOR ? -1 : 0));
            const std::array<std::array<std::size_t,5>,4> score = {std::array<std::size_t,5>{ 0,   2,  4, 12, 16},
                                                                   std::array<std::size_t,5>{ 1,   3,  5, 13, 17},
                                                                   std::array<std::size_t,5>{ 8,   6,  7, 14, 18},
                                                                   std::array<std::size_t,5>{ 11, 10,  9, 15, 19}};
            const std::size_t thisCompositeScore  = score[int(this->intersection)][thisAreaScore];
            const std::size_t otherCompositeScore = score[int(other.intersection)][otherAreaScore];
            if (thisCompositeScore < otherCompositeScore) return true;
            if (this->triQuality > other.triQuality) return true;
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
        Vec2 midPt;      // Midpoint of edge, which is center of diametral circle
        double radSq;    // Diametral Circle squared-radius
    };




protected:    // static functions


    static double distanceSquared(const Vec2& a, const Vec2& b)
    {
        const double dx = b[0] - a[0];
        const double dy = b[1] - a[1];
        return (dx * dx + dy * dy);
    }


    static PointInPolygon PointInTriangleTest2D(const Vec2& a, const Vec2& p, const Vec2& q, const Vec2& r)
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
        if (abs(sum) == 3) return PointInPolygon::INTERIOR;
        if (abs(sum) == 2) return PointInPolygon::BOUNDARY;    // On edge
        if (abs(sum) == 1 && prod == 0) return PointInPolygon::BOUNDARY;    // On vertex
        return PointInPolygon::EXTERIOR;
    }


    // Compute signed area (+ is a good non-reflex vertex, - is a reflex vertex, near-0 is hard to say)
    static double triangleSignedArea(const Vec2& p, const Vec2& q, const Vec2& r)
    {
        return 0.5 * ((p[0] * q[1] - p[1] * q[0]) + (q[0] * r[1] - q[1] * r[0]) + (r[0] * p[1] - r[1] * p[0]));
    }


    // Use winding area to determine is a vertex is a reflex vertex or not
    static PointInPolygon isReflexVertex(const std::size_t& ear, const std::vector<Ear>& polygon)
    {
        int winding = 0;
        const Vec2& a = polygon[ear].pos;
        const Vec2& q = polygon[polygon[ear].prevEar].pos;
        const Vec2& r = polygon[polygon[ear].nextEar].pos;
        double inwardness = (r[0] - q[0]) * (a[1] - q[1]) - (r[1] - q[1]) * (a[0] - q[0]);
        if (abs(inwardness) < TRI_EPSILON)
        {
            const double u = abs(r[0] - q[0]) > abs(r[1] - q[1]) ? (a[0] - q[0]) / (r[0] - q[0]) : (a[1] - q[1]) / (r[1] - q[1]);
            if (-TRI_EPSILON < u && u < 1.0 + TRI_EPSILON) return PointInPolygon::BOUNDARY;
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
                if (-TRI_EPSILON < u && u < 1.0 + TRI_EPSILON) return PointInPolygon::BOUNDARY;
                if (q[1] <= a[1] && a[1] < r[1] && inwardness > 0.0) winding++;
                if (r[1] <= a[1] && a[1] < q[1] && inwardness < 0.0) winding--;
            }
            vertex = polygon[vertex].nextEar;
        }
        return (winding == 0 ? PointInPolygon::EXTERIOR : PointInPolygon::INTERIOR);
    }


    // Test for reflex vertex using area primarily, but if abs(area)<EPS, then use winding order test
    static PointInPolygon isReflexVertexFast(const std::size_t& ear, const std::vector<Ear>& polygon, double& signedArea)
    {
        const Vec2& p = polygon[polygon[ear].prevEar].pos;
        const Vec2& q = polygon[ear].pos;
        const Vec2& r = polygon[polygon[ear].nextEar].pos;
        signedArea = triangleSignedArea(p, q, r);
        if (abs(signedArea) <= TRI_EPSILON) return isReflexVertex(ear, polygon);
        if (signedArea > 0.0) return PointInPolygon::EXTERIOR;
        return PointInPolygon::INTERIOR;
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



