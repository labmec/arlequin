#ifdef HAVE_CONFIG_H
  #include <pz_config.h>
#endif

#include "TPZGenGrid2D.h"

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdio>

#include <TPZGeoMeshTools.h>
#include <TPZVTKGeoMesh.h>
#include "TPZMultiphysicsCompMesh.h"
#include "pzlog.h"
#include "pzbuildmultiphysicsmesh.h"
#include <TPZLinearAnalysis.h>
#include <pzstepsolver.h> //for TPZStepSolver
#include <TPZSSpStructMatrix.h> //symmetric sparse matrix storage
#include <TPZGmshReader.h>
#include <pzskylstrmatrix.h> //symmetric skyline matrix storage
#include "arlequin_config.h"
#include "TPZRefPattern.h"
#include "TPZRefLinear.h"
#include "pzgeoel.h"
#include "ArlequinApproxSpaceCreator.h"
#include "TPZRefPatternDataBase.h"



/**
   @brief Reads the test mesh from gmsh
   @param[in] file_name the .msh mesh file.
*/
TPZGeoMesh* ReadMeshFromGmsh(std::string file_name);

TPZCompMesh* CreateModelCMesh(TPZGeoMesh* gmesh,std::set<int> &allMat);

TPZCompMesh* CreateGluingCMesh(TPZGeoMesh* gmesh);

TPZMultiphysicsCompMesh* CreateArlequinModel(TPZGeoMesh* gmesh, TPZVec<TPZCompMesh *> &meshvector);
void PrintGeoMesh(TPZGeoMesh *gmesh);

void ChangeInternalOrder(TPZCompMesh *cmesh, int pOrder);

//Analytical solution
constexpr int solOrder{2};
auto exactSol = [](const TPZVec<REAL> &loc,
    TPZVec<STATE>&u,
    TPZFMatrix<STATE>&gradU){
    const auto &x=loc[0];

    u[0] = x/2*(3-x);
    gradU(0,0) = 1/2*(3-x) - x/2;
    
};

//-------------------------------------------------------------------------------------------------
//   __  __      _      _   _   _     
//  |  \/  |    / \    | | | \ | |
//  | |\/| |   / _ \   | | |  \| |
//  | |  | |  / ___ \  | | | |\  |
//  |_|  |_| /_/   \_\ |_| |_| \_|
//-------------------------------------------------------------------------------------------------
using namespace std;

int main(int argc, char* argv[]){
    const int dim{2};
    const int pOrder{1};

#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    // gRefDBase.InitializeRefPatterns();

    gRefDBase.InitializeUniformRefPattern(EOned);


    //Create Geometric meshes
    TPZGeoMesh* gmeshCoarse = ReadMeshFromGmsh(string(MESHDIR)+"bar2d.msh");
   
    ArlequinApproxSpaceCreator arlequinCreator(gmeshCoarse,pOrder);
    arlequinCreator.CreateGeoElements();
    gmeshCoarse = arlequinCreator.GeoMesh();
    // arlequinCreator.RefineLocalModel();

    // Cria a malha multifísica com as 3 malhas usando interface CompEl.
    TPZMultiphysicsCompMesh *cmeshmulti = arlequinCreator.ArlequinCompMesh();
    std::string cmesh3 = "CMeshArlequin.txt";
    std::ofstream myfile3(cmesh3);
    cmeshmulti->Print(myfile3);

    std::string vtk_name = "compMesh.vtk";
    std::ofstream vtkfile(vtk_name.c_str());
    TPZVTKGeoMesh::PrintCMeshVTK(cmeshmulti, vtkfile, true);

    // //Create analysis environment
    TPZLinearAnalysis an;
    an.SetCompMesh(cmeshmulti,false);
    TPZSkylineStructMatrix<REAL> matskl(cmeshmulti);
    an.SetStructuralMatrix(matskl);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);
    an.Run();
    return 0;
}


TPZGeoMesh* ReadMeshFromGmsh(std::string file_name)
{
    //read mesh from gmsh
    TPZGeoMesh *gmesh;
    gmesh = new TPZGeoMesh();
    {
        TPZGmshReader reader;
        // essa interface permite voce mapear os nomes dos physical groups para
        // o matid que voce mesmo escolher
        TPZManVector<std::map<std::string,int>,4> stringtoint(4);
        stringtoint[1]["Global"] = EGlobal;
        stringtoint[1]["Gluing"] = EGluing;
        stringtoint[1]["Local"] = ELocal;

        reader.SetDimNamePhysical(stringtoint);
        reader.GeometricGmshMesh(file_name,gmesh);

        //Prints gmesh mesh properties
        std::string vtk_name = "geoMesh.vtk";
        std::ofstream vtkfile(vtk_name.c_str());

        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, true);
    }

    return gmesh;
}

