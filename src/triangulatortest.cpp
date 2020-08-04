#include "triangulatortest.h"




using namespace Triangulator;
using namespace TriangulatorTest;




int TriangulatorTest::exportContourToImage(const std::vector<Vec2>& contour,
                                           const QString& imageFilePath,
                                           const int& imageWidth)
{
    // Set some image parameters
    const int IMAGE_WIDTH      = imageWidth;
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
        painter.drawEllipse(points[0], 5, 5);
        //painter.drawText(points[0], QString::number(i));
    }

    // Finalize painter and save image
    painter.end();
    int result = image.save(imageFilePath);
    return result;
}




int TriangulatorTest::exportTriangulationToImage(const std::vector<Vec2>& contour,
                                                 std::vector<std::array<std::size_t,3> >& triangles,
                                                 const QString& imageFilePath,
                                                 const int& imageWidth)
{
    // Set some image parameters
    const int IMAGE_WIDTH      = imageWidth;
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
                                                 const QString& imageFilePath,
                                                 const int& imageWidth)
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

    return exportTriangulationToImage(singleContour, trianglesFlat, imageFilePath, imageWidth);
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




std::vector<Vec2>
TriangulatorTest::createRandomPolygon2(const int& startingSeed, int& nextSeed, const std::size_t& fixedSize, const int& timeout)
{
    nextSeed = startingSeed;
    bool success = false;
    while (!success && nextSeed - startingSeed < timeout)
    {
        std::vector<Vec2> contour = createRandomPolygon(nextSeed, success);
        if (8 <= fixedSize && fixedSize <=16 && contour.size() != fixedSize) success = false;
        nextSeed++;
        if (success) return contour;
    }
    return std::vector<Vec2>();
}




int TriangulatorTest::generateDeepLearningDataset1(const QString& pathName,
                                                   const QString& imageFolderName,
                                                   const QString& csvFileName,
                                                   const std::size_t& numTriangulations,
                                                   const int& startingSeed,
                                                   int& nextSeed,
                                                   const std::size_t& fixedSize)
{
    // Name and open output files
    QFile csvFile( QDir::cleanPath(pathName + QDir::separator() + csvFileName) );
    if (csvFile.open(QIODevice::WriteOnly))
    {
        std::cerr << "CSV Output File opened successfully\n";
    }
    else {
        std::cerr << "CSV Output File could not be opened in WriteOnly mode\n";
        std::cerr << csvFile.errorString().toStdString() <<"\n";
        return -1;
    }
    QDir imageDir( QDir::cleanPath(pathName + QDir::separator() + imageFolderName) );
    if (!imageDir.exists()) QDir().mkpath(imageDir.absolutePath());

    // Randomly generate contours and triangulate them
    nextSeed = startingSeed;
    std::size_t t = 0;
    while (t < numTriangulations)
    {
        // Create random polygon and then triangulate it
        std::vector<Vec2> contour = createRandomPolygon2(nextSeed, nextSeed, fixedSize);
        Diagnostics diagnostics;
        std::vector<std::array<std::size_t,3> > triangles = triangulate1(contour, diagnostics, true, true);
        if (std::abs(diagnostics.areaDiff) > TRI_EPSILON || diagnostics.triDiff != 0) continue;

        // Scale the contour to the [0,1]x[0,1] interval and pad with (-1,-1)
        Vec2 lower = contour[0];
        Vec2 upper = contour[0];
        for (std::size_t i = 1; i < contour.size(); i++)
        {
            lower[0] = std::min(lower[0], contour[i][0]);
            lower[1] = std::min(lower[1], contour[i][1]);
            upper[0] = std::max(upper[0], contour[i][0]);
            upper[1] = std::max(upper[1], contour[i][1]);
        }
        const double width = std::max(upper[0] - lower[0], upper[1] - lower[1]);    // Has to use the same x and y scale to retain Delaunay property
        for (std::size_t i = 0; i < contour.size(); i++)
        {
            contour[i][0] = (contour[i][0] - lower[0]) / width;
            contour[i][1] = (contour[i][1] - lower[1]) / width;
        }
        const std::size_t N = (8 <= fixedSize && fixedSize <= 16 ? fixedSize : 16);    //Bigger number, often 16
        const std::size_t M = contour.size();    // Smaller number, between 8 and 16 - Need to save this value before contour.size changes
        for (std::size_t i = contour.size(); i < N; i++)
        {
            contour.push_back({-1.0, -1.0});    // If the contour is shorter than 16, then pad with (-1,-1)
        }

        // Find the third vertex that makes the triangle with contour[i] and contour[(i+1)%N]
        std::vector<Vec2> thirdVertex(M);
        for (std::size_t i = 0; i < M; i++)
        {
            for (std::size_t j = 0; j < triangles.size(); j++)
            {
                std::size_t match = 0;
                std::size_t other = TRI_UNDEFINED_INDEX;
                for (std::size_t k = 0; k < 3; k++)
                {
                    if (triangles[j][k] == i || triangles[j][k] == (i+1)%M)
                    {
                        match++;
                    }
                    else {
                        other = triangles[j][k];
                    }
                }
                if (match == 2)
                {
                    thirdVertex[i] = contour[other];
                    break;
                }
            }
        }
        for (std::size_t i = thirdVertex.size(); i < N; i++)
        {
            thirdVertex.push_back({-1.0, -1.0});    // If the thirdVertex vector is shorter than 16, then pad with (-1,-1)
        }

        // Write the contour and thirdVertex to the csv file
        QTextStream csvStream(&csvFile);
        for (std::size_t i = 0; i < contour.size(); i++)
        {
            csvStream << contour[i][0] <<","<< contour[i][1] <<",";
        }
        for (std::size_t i = 0; i < thirdVertex.size(); i++)
        {
            csvStream << thirdVertex[i][0] <<","<< thirdVertex[i][1] <<",";
        }
        csvStream << "\n";

        if (false)
        {
            // Export triangulation to image to help with debugging the deep learning model
            QString imageFilePath = QDir::cleanPath(pathName + QDir::separator() +
                                                    imageFolderName + QDir::separator() +
                                                    "Triangulation" + QString::number(t) + ".png");
            int result = exportTriangulationToImage(contour, triangles, imageFilePath, 256);
            if (result < 0)
            {
                csvFile.close();
                return result;
            }
        }

        t++;    // Increment triangulation counter
    }

    csvFile.close();
    return 0;
}




int TriangulatorTest::generateDeepLearningDataset(const QString& baseDir,
                                                  const std::size_t& nTrain,
                                                  const std::size_t& nTest,
                                                  const std::size_t& nVal)
{
    int seed = 1;
    int result1 = TriangulatorTest::generateDeepLearningDataset1(QDir::cleanPath(baseDir + QDir::separator() + "Training/"),
                                                                 "Images/", "Triangulations.csv",
                                                                 nTrain, seed, seed, 8);
    if (result1 < 0) return result1;

    int result2 = TriangulatorTest::generateDeepLearningDataset1(QDir::cleanPath(baseDir + QDir::separator() + "Test/"),
                                                                 "Images/", "Triangulations.csv",
                                                                 nTest, seed, seed, 8);
    if (result2 < 0) return result2;

    int result3 = TriangulatorTest::generateDeepLearningDataset1(QDir::cleanPath(baseDir + QDir::separator() + "Validation/"),
                                                                 "Images/", "Triangulations.csv",
                                                                 nVal, seed, seed, 8);
    if (result3 < 0) return result3;

    return 0;
}



