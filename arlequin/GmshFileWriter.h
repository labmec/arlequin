/**
 * @file GmshFileWriter.h
 * @brief Contains the GmshFileWriter class writes a .geo file to be used as an input by gmsh
 */

#ifndef GMSHWRITER_H
#define GMSHWRITER_H

#include <iostream>
#include "pzgeoelside.h"
#include "pzgnode.h"
#include "pzeltype.h"

class GmshFileWriter{

public:
    //! Default constructor
    GmshFileWriter() = default;

    void WriteGeoPoint(std::ofstream &file, TPZGeoNode &point, REAL refin=1);
    
    void WriteGeoLine(std::ofstream &file, TPZGeoElSide &side, bool isPhysicalEntity, int nDivisions = 3);
    
    void WritePhysicalLine(std::ofstream &file, TPZGeoElSide &side);

    void WriteGeoSurface(std::ofstream &file, TPZGeoElSide &side, int nfacets, MElementType &elType, bool isPhysicalEntity);

    void WritePhysicalSurface(std::ofstream &file, TPZGeoElSide &side);

private:

    int fLineCounter = 1;

};

#endif