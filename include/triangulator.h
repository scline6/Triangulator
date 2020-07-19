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
#include <iostream>




const double      TRI_EPSILON                = 1e-12;
const std::size_t TRI_UNDEFINED_INDEX        = std::numeric_limits<std::size_t>::max();
const std::size_t TRI_SPATIAL_HASH_CUTOFF    = 20;
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




protected:


    enum PointInTriangle
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
        bool operator<(const Ear& other) { return this->hash < other.hash; }
    };


//    Delicate priority logic - lower numbers are higher priority
//
//                    earArea      Positive Degenerate Negative
//    intersection                    (0)       (1)       (2)
//
//    EXTERIOR        (0)              0         3         8
//    NOT_TESTED_YET  (1)              1         4         9
//    BOUNDARY        (2)              2         5        10
//    INTERIOR        (3)              7         6        11


    struct Priority
    {
        PointInTriangle intersection;
        double signedArea;
        double triQuality;
        Priority(const PointInTriangle& pit, const double& sa, const double& tq) : intersection(pit), signedArea(sa), triQuality(tq) {}
        bool operator<(const Priority& other) const
        {
            int thisAreaScore  = (std::abs(this->signedArea) <= TRI_EPSILON ? 1 : (this->signedArea > 0.0 ? 0 : 2));    //0 positive, 1 degen, 2 negative
            int otherAreaScore = (std::abs(other.signedArea) <= TRI_EPSILON ? 1 : (other.signedArea > 0.0 ? 0 : 2));    //0 positive, 1 degen, 2 negative
            const std::array<std::array<std::size_t,3>,4> score = {std::array<std::size_t,3>{ 0,  3,  8},
                                                                   std::array<std::size_t,3>{ 1,  4,  9},
                                                                   std::array<std::size_t,3>{ 2,  5, 10},
                                                                   std::array<std::size_t,3>{ 7,  6, 11}};
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
            //if (it == this->order.end()) return;
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


    static PointInTriangle pointInTriangleTest2D(const Vec2& a, const Vec2& p, const Vec2& q, const Vec2& r)
    {
        const double qrTest = (a[0] - r[0]) * (q[1] - r[1]) - (q[0] - r[0]) * (a[1] - r[1]);
        const double rpTest = (a[0] - p[0]) * (r[1] - p[1]) - (r[0] - p[0]) * (a[1] - p[1]);
        const double pqTest = (a[0] - q[0]) * (p[1] - q[1]) - (p[0] - q[0]) * (a[1] - q[1]);
        const bool allNegative = (qrTest < -TRI_EPSILON) && (rpTest < -TRI_EPSILON) && (pqTest < -TRI_EPSILON);
        const bool allPositive = (qrTest >  TRI_EPSILON) && (rpTest >  TRI_EPSILON) && (pqTest >  TRI_EPSILON);
        if (allNegative || allPositive) return PointInTriangle::INTERIOR;
        const bool allNonPositive = (qrTest <  TRI_EPSILON) && (rpTest <  TRI_EPSILON) && (pqTest <  TRI_EPSILON);
        const bool allNonNegative = (qrTest > -TRI_EPSILON) && (rpTest > -TRI_EPSILON) && (pqTest > -TRI_EPSILON);
        if (allNonPositive || allNegative) return PointInTriangle::BOUNDARY;
        return PointInTriangle::EXTERIOR;
    }


    static double triangleSignedArea(const Vec2& p, const Vec2& q, const Vec2& r)
    {
        return 0.5 * ((p[0] * q[1] - p[1] * q[0]) + (q[0] * r[1] - q[1] * r[0]) + (r[0] * p[1] - r[1] * p[0]));
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



