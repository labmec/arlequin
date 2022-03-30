#include "GmshFileWriter.h"
#include <fstream>


void GmshFileWriter::WriteGeoPoint(std::ofstream &file, TPZGeoNode &point, REAL refin){
    // std::cout << "WRITING POINT " << point.Id() << std::endl;
    file << "Point(" << point.Id() << ") = {" << point.Coord(0) << ", " 
                                                << point.Coord(1) << ", " 
                                                << point.Coord(2) << ", " << refin << "};\n\n"; 
}
void GmshFileWriter::WriteGeoLine(std::ofstream &file, TPZGeoElSide &side, bool isPhysicalEntity, int nDivisions){
    // std::cout << "WRITING LINE " << side.Id() << std::endl;
    if (isPhysicalEntity){
        file << "Line(" << side.Id() << ") = {" << side.SideNodeIndex(0) << ", " << side.SideNodeIndex(1) << "};\n\n"; 
    } else {
        file << "Line(" << fLineCounter << ") = {" << side.SideNodeIndex(0) << ", " << side.SideNodeIndex(1) << "};\n\n"; 
        file << "Transfinite Curve {" << fLineCounter << "} = " << nDivisions << " Using Progression 1; \n";
        fLineCounter++;
    }
    
}
void GmshFileWriter::WritePhysicalLine(std::ofstream &file, TPZGeoElSide &side){
    std::cout << "WRITING Physical LINE " << side.Id() << std::endl;
    file << "Physical Curve(\"" << 100+side.Id() << "\") = {" << side.Id() << "};\n\n"; 
}

void GmshFileWriter::WriteGeoSurface(std::ofstream &file, TPZGeoElSide &side, int nfacets, MElementType &elType, bool isPhysicalEntity){
    // std::cout << "WRITING LINE " << side.Id() << std::endl;
    if (isPhysicalEntity){
        file << "Curve Loop(" << side.Id() << ") = { ";
        for (int i = 0; i < nfacets; i++)
        {
            file << fLineCounter - nfacets + i << ", " ;
        }
        file.seekp(-2, std::ios_base::cur);
        file << " };\n\n"; 
        file << "Plane Surface(" << side.Id() << ") = {" << side.Id()<< "};\n";
        if (elType == EQuadrilateral){
            file << "Recombine Surface{" << side.Id() << "};\n";
        }
    } else {
        DebugStop();
        // file << "Line(" << fLineCounter << ") = {" << side.SideNodeIndex(0) << ", " << side.SideNodeIndex(1) << "};\n\n"; 
        // fLineCounter++;
    }
    
}

void GmshFileWriter::WritePhysicalSurface(std::ofstream &file, TPZGeoElSide &side){
    std::cout << "WRITING Physical SURFACE " << side.Id() << std::endl;
    file << "Physical Surface(\"" << 100+side.Id() << "\") = {" << side.Id() << "};\n\n"; 
}