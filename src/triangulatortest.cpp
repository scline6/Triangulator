#include "triangulatortest.h"




using namespace Triangulator;
using namespace TriangulatorTest;




int TriangulatorTest::exportContourToImage(const std::vector<Vec2>& contour,
                                           const QString& imageFilePath)
{
    // Set some image parameters
    const int IMAGE_WIDTH      = 1000;
    const int LINE_WIDTH       = 3;
    const double GROWTH_FACTOR = 1.01;

    // Compute the bounding box for the polygon and triangulation
    Vec2 lower = contour[0];
    Vec2 upper = contour[0];
    for (std::size_t i = 0; i < contour.size(); i++)
    {
        lower[0] = std::min(lower[0], contour[i][0]);
        lower[1] = std::min(lower[1], contour[i][1]);
        upper[0] = std::max(upper[0], contour[i][0]);
        upper[1] = std::max(upper[1], contour[i][1]);
    }
    const Vec2 center = {0.5 * (lower[0] + upper[0]), 0.5 * (lower[1] + upper[1])};
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
    painter.setFont( QFont("Arial", 30, QFont::Bold) );
    std::array<QColor,12> palette = {Qt::red, Qt::green, Qt::blue,
                                     Qt::cyan, Qt::magenta, Qt::yellow,
                                     Qt::darkRed, Qt::darkGreen, Qt::darkBlue,
                                     Qt::darkCyan, Qt::darkMagenta, Qt::darkYellow};

    // Draw contour
    for (std::size_t i = 0; i < contour.size(); i++)
    {
        const QColor color = palette[i % 12];
        painter.setBrush( QBrush(color, Qt::BrushStyle::SolidPattern) );
        QPointF points[2];
        points[0].setX( imageHalfWidth + scaleFactor * (contour[i][0] - center[0]) );
        points[0].setY( imageHalfWidth - scaleFactor * (contour[i][1] - center[1]) );
        const std::size_t j = (i+1)%contour.size();
        points[1].setX( imageHalfWidth + scaleFactor * (contour[j][0] - center[0]) );
        points[1].setY( imageHalfWidth - scaleFactor * (contour[j][1] - center[1]) );
        painter.drawLine(points[0], points[1]);
        painter.drawEllipse(points[0], 9, 9);
        //painter.drawText(points[0], QString::number(i));
    }

    // Finalize painter and save image
    painter.end();
    int result = image.save(imageFilePath);
    return result;
}




int TriangulatorTest::exportTriangulationToImage(const std::vector<Vec2>& contour,
                                                 std::vector<std::array<std::size_t,3> >& triangles,
                                                 const QString& imageFilePath)
{
    // Set some image parameters
    const int IMAGE_WIDTH      = 1000;
    const int LINE_WIDTH       = 3;
    const double GROWTH_FACTOR = 1.01;

    // Compute the bounding box for the polygon and triangulation
    Vec2 lower = contour[0];
    Vec2 upper = contour[0];
    for (std::size_t i = 0; i < contour.size(); i++)
    {
        lower[0] = std::min(lower[0], contour[i][0]);
        lower[1] = std::min(lower[1], contour[i][1]);
        upper[0] = std::max(upper[0], contour[i][0]);
        upper[1] = std::max(upper[1], contour[i][1]);
    }
    const Vec2 center = {0.5 * (lower[0] + upper[0]), 0.5 * (lower[1] + upper[1])};
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




int TriangulatorTest::exportTriangulationToImage(const std::vector<std::vector<Vec2> >& contours,
                                                 std::vector<std::array<std::array<std::size_t,2>,3> >& triangles,
                                                 const QString& imageFilePath)
{
    std::size_t numVerticesTotal = 0;
    for (std::size_t i = 0; i < contours.size(); i++)
    {
        numVerticesTotal += contours[i].size();
    }

    std::vector<Vec2> singleContour(numVerticesTotal);
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

    return exportTriangulationToImage(singleContour, trianglesFlat, imageFilePath);
}




std::vector<std::array<std::size_t,2> >
TriangulatorTest::testDelaunayProperty(const std::vector<Vec2>& contour,
                                       const std::vector<std::array<std::size_t,3> >& triangles)
{
    std::vector<std::array<std::size_t,2> > encroachedList;
    for (std::size_t i = 0; i < triangles.size(); i++)
    {
        const Vec2& p = contour[triangles[i][0]];
        const Vec2& q = contour[triangles[i][1]];
        const Vec2& r = contour[triangles[i][2]];
        Circumcircle cc = triangleCircumcircle(p, q, r);
        for (std::size_t j = 0; j < contour.size(); j++)
        {
            if (j == triangles[i][0] || j == triangles[i][1] || j == triangles[i][2]) continue;
            const double distSq = distanceSquared(cc.center, contour[j]);
            if (distSq < cc.radSq) encroachedList.push_back({i, j});
        }
    }
    return encroachedList;
}




std::vector<Vec2>
TriangulatorTest::createRandomPolygon(const int& seed, bool& success)
{
    const std::size_t CONTOUR_MAX_SIZE = 16;
    const std::size_t CONTOUR_MIN_SIZE =  8;
    const std::size_t MAX_ATTEMPTS     = 10000;
    const double      MIN_DISTANCE     =  1.0;
    const double      MAX_DISTANCE     = 10.0;
    const double      RANGE            = MAX_DISTANCE - MIN_DISTANCE;
    const Vec2        STARTING_VERTEX  = {0.0, 0.0};

    success = false;
    std::vector<Vec2> contour;
    contour.reserve(CONTOUR_MAX_SIZE);
    contour.push_back(STARTING_VERTEX);
    bool tryToGoBackToStart = false;
    std::size_t tryCount = 0;
    srand(seed);
    while (contour.size() < CONTOUR_MAX_SIZE && tryCount < MAX_ATTEMPTS)
    {
        Vec2 nextVertex;
        if (tryToGoBackToStart)
        {
            nextVertex = STARTING_VERTEX;
        }
        else {
            const double angle = 2.0 * M_PI * double(rand() % 100000) / double(100000);
            const double distance = MIN_DISTANCE + RANGE * double(rand() % 100000) / double(100000);
            nextVertex = {contour.back()[0] + distance * cos(angle), contour.back()[1] + distance * sin(angle)};
        }

        bool intersectsSomething = false;
        for (std::size_t i = 1; i < contour.size()-1; i++)
        {
            IntersectionResult result = edgeEdgeIntersectionTest2D(contour.back(), nextVertex, contour[i-1], contour[i], TRI_EPSILON);
            if (result == IntersectionResult::INTERIOR)
            {
                intersectsSomething = true;
                break;
            }
        }

        if (intersectsSomething)
        {
            tryToGoBackToStart = false;
            tryCount++;
            continue;
        }
        else {
            if (tryToGoBackToStart)
            {
                success = true;
                break;    // This means it got back to the start, so quit
            }
            contour.push_back(nextVertex);
            tryToGoBackToStart = (contour.size() < CONTOUR_MIN_SIZE ? false : true);
        }
    }
    return std::move(contour);
}



