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

#include <pzgmesh.h> //for TPZGeoMesh
#include <pzcmesh.h> //for TPZCompMesh
#include <TPZVTKGeoMesh.h>
#include <TPZGeoMeshTools.h>
#include <TPZVTKGeoMesh.h>
#include <Poisson/TPZMatPoisson.h>
#include "TPZNullMaterial.h"
#include "TPZMultiphysicsCompMesh.h"
#include "pzlog.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZMatPoissonCS.h"
#include <TPZLinearAnalysis.h>
#include <pzstepsolver.h> //for TPZStepSolver
#include <TPZSSpStructMatrix.h> //symmetric sparse matrix storage
#include <TPZGmshReader.h>
#include "TPZCompElDisc.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzgeoelside.h"
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "TPZMixedDarcyFlowArlequin.h"
#include "TPZLagrangeMultiplierCS.h"
#include <pzskylstrmatrix.h> //symmetric skyline matrix storage
#include "arlequin_config.h"
#include "TPZRefPattern.h"
#include "TPZGeoElement.h"
#include "TPZRefLinear.h"
#include "pzgeoel.h"
#include "ArlequinGeoMeshCreator.h"

enum EMatId {ENone, EGlobal, EGluing, ELocal};

/**
   @brief Reads the test mesh from gmsh
   @param[in] file_name the .msh mesh file.
*/
TPZGeoMesh* ReadMeshFromGmsh(std::string file_name);

TPZCompMesh* CreateModelCMesh(TPZGeoMesh* gmesh,std::set<int> &allMat);

TPZCompMesh* CreateGluingCMesh(TPZGeoMesh* gmesh);

TPZMultiphysicsCompMesh* CreateArlequinModel(TPZGeoMesh* gmesh, TPZVec<TPZCompMesh *> &meshvector);
void PrintGeoMesh(TPZGeoMesh *gmesh);



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

    //Create Geometric meshes
    TPZGeoMesh* gmeshCoarse = ReadMeshFromGmsh(string(MESHDIR)+"bar2d.msh");
    // PrintGeoMesh(gmeshCoarse);
    
    ArlequinGeoMeshCreator ArlequinCreator;
    ArlequinCreator.SetMaterialIds(EGluing,EGlobal,ELocal);

    std::set<int> overlapMatId = {EGluing};
    TPZGeoMesh* gmeshFine = ArlequinCreator.CreateFineGeoMesh(gmeshCoarse, overlapMatId);
    // PrintGeoMesh(gmeshFine);
    TPZGeoMesh* gmeshArlequin = ArlequinCreator.AssociateModels(gmeshCoarse,gmeshFine,overlapMatId);
    delete gmeshCoarse, gmeshFine;
    PrintGeoMesh(gmeshArlequin);

    TPZVec<TPZCompMesh *> meshvector(3); 

    //Cria uma malha computacional apenas com os material ids dos elementos globais.
    std::set<int> allMat={EGlobal};
    meshvector[0] = CreateModelCMesh(gmeshArlequin,allMat);
    std::string cmesh0 = "CMesh0.txt";
    std::ofstream myfile0(cmesh0);
    meshvector[0]->Print(myfile0);

    //Cria uma malha computacional apenas com os material ids dos elementos locais.
    std::set<int> allMat2={ELocal};
    meshvector[1] = CreateModelCMesh(gmeshArlequin,allMat2);
    std::string cmesh1 = "CMesh1.txt";
    std::ofstream myfile1(cmesh1);
    meshvector[1]->Print(myfile1);

    //Cria uma malha de pressão com o material id da zona de colagem
    meshvector[2] = CreateGluingCMesh(gmeshArlequin);
    std::string cmesh2 = "CMesh2.txt";
    std::ofstream myfile2(cmesh2);
    meshvector[2]->Print(myfile2);

    //Cria a malha multifísica com as 3 malhas usando interface CompEl.
    TPZMultiphysicsCompMesh *cmesh = CreateArlequinModel(gmeshArlequin,meshvector);
    std::string cmesh3 = "CMeshArlequin.txt";
    std::ofstream myfile3(cmesh3);
    cmesh->Print(myfile3);

    //Create analysis environment
    TPZLinearAnalysis an(meshvector[0],false);
    TPZSkylineStructMatrix<REAL> matskl(cmesh);
    an.SetStructuralMatrix(matskl);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);
    an.Assemble();
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




TPZCompMesh* CreateModelCMesh(TPZGeoMesh* gmesh, std::set<int> &allMat){

    gmesh->ResetReference();
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(1);
    cmesh->SetDimModel(gmesh->Dimension());

    for (std::set<int>::iterator it=allMat.begin(); it!=allMat.end(); ++it){
        TPZNullMaterial<> *mat = new TPZNullMaterial<>(*it);
        cmesh->InsertMaterialObject(mat);
        mat->SetDimension(1);
    } 

    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild();

    return cmesh;
}

TPZCompMesh* CreateGluingCMesh(TPZGeoMesh* gmesh){

    gmesh->ResetReference();
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);

    TPZNullMaterial<> *mat = new TPZNullMaterial<>(EGluing);
    mat->SetDimension(1);
    cmesh->InsertMaterialObject(mat);

    cmesh->SetDefaultOrder(1);
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    cmesh->SetDimModel(1);
    cmesh->AutoBuild();

    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i]; 
        newnod.SetLagrangeMultiplier(1);
    }

    int nel = cmesh->NElements();
    for(int i=0; i<nel; i++){
        TPZCompEl *cel = cmesh->ElementVec()[i];
        TPZCompElDisc *celdisc = dynamic_cast<TPZCompElDisc *>(cel);
        if(!celdisc) continue;
        celdisc->SetConstC(1.);
        celdisc->SetTrueUseQsiEta();
        // espera-se elemento de pressao apenas para o contorno
        auto aaa = celdisc->Reference()->Dimension();
        // if(celdisc && celdisc->Reference()->Dimension() == cmesh->Dimension())
        // {
        //     DebugStop();
        // }
    }

    return cmesh;
}

TPZMultiphysicsCompMesh* CreateArlequinModel(TPZGeoMesh* gmesh, TPZVec<TPZCompMesh *> &meshvector){

    auto cmesh = new TPZMultiphysicsCompMesh(gmesh);
    cmesh->SetDefaultOrder(1);
    cmesh->SetDimModel(1);

    auto mat = new TPZMixedDarcyFlowArlequin(EGlobal,1);
    cmesh->InsertMaterialObject(mat);
    auto mat2 = new TPZMixedDarcyFlowArlequin(ELocal,1);
    cmesh->InsertMaterialObject(mat2);
    // auto mat3 = new TPZMatPoissonCS(EGluing,1);
    // cmesh->InsertMaterialObject(mat3);
    // auto mat4 = new TPZMatPoissonCS(10,1);
    // cmesh->InsertMaterialObject(mat4);

    // auto mat3 = new TPZLagrangeMultiplierCS<STATE>(EGluing, 1);
    // cmesh->InsertMaterialObject(mat3);

    TPZManVector<int> active(meshvector.size(),1);
    // active[2]=0;
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    cmesh->BuildMultiphysicsSpace(active, meshvector);
    cmesh->LoadReferences();
    cmesh->CleanUpUnconnectedNodes(); 

    return cmesh;
}


void PrintGeoMesh(TPZGeoMesh *gmesh){
    
    for (int i = 0; i < gmesh->NElements(); i++)
    {
        auto *gel = gmesh->Element(i);
        if (!gel) continue;

        int matid = gel->MaterialId();
        auto nsides = gel->NSides();
        // auto nconnects = gel->Reference()->NConnects();
        std::cout << "ELGeometric = " << i << ", dim= " << gel->Dimension() << ",mat = " << gel->MaterialId() << std::endl;
  
        nsides = gel->NSides();
        int ncorner = gel->NCornerNodes();
        for (int side = 0; side < nsides; side++) {
            // if(gel->SideDimension(side) != 1) continue;
            TPZGeoElSide gelside(gel,side);
            TPZGeoElSide neighbour = gelside.Neighbour();
            
            std::cout << "Element = " << i << ", side = " << side  
                    << ", NEL = " << neighbour.Element()->Index() 
                    << ", Nmatid = " << neighbour.Element()->MaterialId()
                    << ", NNEL = " << neighbour.Neighbour().Element()->Index() 
                    << ", NNmatid = " << neighbour.Neighbour().Element() -> MaterialId() << std::endl;
        }
    }

    //Prints gmesh mesh properties
    std::string vtk_name = "geoMesh.vtk";
    std::ofstream vtkfile(vtk_name.c_str());

    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, true);

}


