#ifndef TRIANGULATORTEST_H
#define TRIANGULATORTEST_H




#include "triangulator.h"
#include <iostream>
#include <stdlib.h>
#include <QImage>
#include <QPainter>
#include <QFile>
#include <QFileInfo>
#include <QDir>
#include <QTextStream>




using namespace Triangulator;




namespace TriangulatorTest
{

    int exportContourToImage(const std::vector<Vec2>& contour,
                             const QString& imageFilePath,
                             const int& imageWidth = 1000);


    int exportTriangulationToImage(const std::vector<Vec2>& contour,
                                  std::vector<std::array<std::size_t,3> >& triangles,
                                  const QString& imageFilePath,
                                  const int& imageWidth = 1000);


    int exportTriangulationToImage(const std::vector<std::vector<Vec2> >& contours,
                                  std::vector<std::array<std::array<std::size_t,2>,3> >& triangles,
                                  const QString& imageFilePath,
                                  const int& imageWidth = 1000);


    std::vector<std::array<std::size_t,2> >
    testDelaunayProperty(const std::vector<Vec2>& contour,
                         const std::vector<std::array<std::size_t,3> >& triangles);


    std::vector<Vec2>
    createRandomPolygon(const int& seed,
                        bool& success);


    std::vector<Vec2>
    createRandomPolygon2(const int& startingSeed,
                         int& nextSeed,
                         const std::size_t& fixedSize = 0,
                         const int& timeout = 100000);


    int generateDeepLearningDataset1(const QString& pathName,
                                     const QString& imageFolderName,
                                     const QString& csvFileName,
                                     const std::size_t& numTriangulations,
                                     const int& startingSeed,
                                     int& nextSeed,
                                     const std::size_t& fixedSize = 0);


    int generateDeepLearningDataset(const QString& baseDir,
                                    const std::size_t& nTrain = 50000,
                                    const std::size_t& nTest  = 15000,
                                    const std::size_t& nVal   =  5000);

}




#endif // TRIANGULATORTEST_H



