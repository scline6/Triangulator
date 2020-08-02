#ifndef TRIANGULATORTEST_H
#define TRIANGULATORTEST_H




#include "triangulator.h"
#include <iostream>
#include <stdlib.h>
#include <QImage>
#include <QPainter>
#include <QFile>
#include <QFileInfo>




using namespace Triangulator;




namespace TriangulatorTest
{

    int exportContourToImage(const std::vector<Vec2>& contour,
                             const QString& imageFilePath);

    int exportTriangulationToImage(const std::vector<Vec2>& contour,
                                  std::vector<std::array<std::size_t,3> >& triangles,
                                  const QString& imageFilePath);


    int exportTriangulationToImage(const std::vector<std::vector<Vec2> >& contours,
                                  std::vector<std::array<std::array<std::size_t,2>,3> >& triangles,
                                  const QString& imageFilePath);


    std::vector<std::array<std::size_t,2> >
    testDelaunayProperty(const std::vector<Vec2>& contour,
                         const std::vector<std::array<std::size_t,3> >& triangles);


    std::vector<Vec2>
    createRandomPolygon(const int& seed, bool& success);


}




#endif // TRIANGULATORTEST_H



