#ifndef TRIANGULATORTEST_H
#define TRIANGULATORTEST_H




#include "triangulator.h"
#include <iostream>
#include <QImage>
#include <QPainter>
#include <QFile>
#include <QFileInfo>




using namespace Triangulator;




namespace TriangulatorTest
{

    int DEBUG__WRITE_IMAGE(const std::vector<Vec2>& contour,
                           std::vector<std::array<std::size_t,3> >& triangles,
                           const QString& imageFilePath);


    int DEBUG__WRITE_IMAGE(const std::vector<std::vector<Vec2> >& contours,
                           std::vector<std::array<std::array<std::size_t,2>,3> >& triangles,
                           const QString& imageFilePath);


    std::vector<std::array<std::size_t,2> >
    DEBUG__TEST_DELAUNAY_PROPERTY(const std::vector<Vec2>& contour,
                                  const std::vector<std::array<std::size_t,3> >& triangles);

}




#endif // TRIANGULATORTEST_H



