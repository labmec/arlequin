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

void CreateGeoMeshes(TPZGeoMesh *gmesh0, TPZGeoMesh *gmesh1);
void CreateCompMeshes(const int &dim, const int &pOrder, TPZGeoMesh *gmesh0, TPZGeoMesh *gmesh1, TPZCompMesh *cmesh0, TPZCompMesh *cmesh1);

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
enum EMatid {ENone, EDomain, EPointLeft0, EPointRight0, EPointLeft1, EPointRight1};


int main(int argc, char* argv[]){
    const int dim{1};
    const int pOrder{1};

#ifdef PZ_LOG
TPZLogger::InitializePZLOG();
#endif

    //Create Geometric meshes
    //Create GeoMesh 1
    const REAL xMin0=0.;
    const REAL xMax0=2.;
    const int nels0 = 4;
    TPZVec<int> matids0 = {EDomain, EPointLeft0, EPointRight0}; 
    TPZGeoMesh *gmesh0 = TPZGeoMeshTools::CreateGeoMesh1D(xMin0,xMax0,nels0,matids0,true);

    std::ofstream out0("gmesh0.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh0, out0);

    //Create GeoMesh 1
    const REAL xMin1=1.;
    const REAL xMax1=3.;
    const int nels1 = 4;
    TPZVec<int> matids1 = {EDomain, EPointLeft1, EPointRight1}; 
    TPZGeoMesh *gmesh1 = TPZGeoMeshTools::CreateGeoMesh1D(xMin1,xMax1,nels1,matids1,true);

    std::ofstream out1("gmesh1.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh1, out1);


    //Create Computational meshes
    TPZCompMesh *cmesh0 = new TPZCompMesh(gmesh0);
    TPZNullMaterial<> *mat0 = new TPZNullMaterial<>(EDomain);
    cmesh0->InsertMaterialObject(mat0);
    TPZNullMaterial<> *mat01 = new TPZNullMaterial<>(EPointLeft0);
    cmesh0->InsertMaterialObject(mat01);
    cmesh0->SetDefaultOrder(pOrder);
    mat0->SetDimension(dim);
    cmesh0->SetAllCreateFunctionsContinuous();
    cmesh0->AutoBuild();

    TPZCompMesh *cmesh1 = new TPZCompMesh(gmesh1);
    TPZNullMaterial<> *mat1 = new TPZNullMaterial<>(EDomain);
    cmesh1->InsertMaterialObject(mat1);
    TPZNullMaterial<> *mat11 = new TPZNullMaterial<>(EPointRight1);
    cmesh1->InsertMaterialObject(mat11);
    cmesh1->SetDefaultOrder(pOrder);
    mat1->SetDimension(dim);
    cmesh1->SetAllCreateFunctionsContinuous();
    cmesh1->AutoBuild();


    //Creates multiphysics mesh with
    //Multiphysics mesh
    TPZManVector< TPZCompMesh *, 1> meshvector(1);
    meshvector[0] = cmesh0;
    // meshvector[1] = cmesh1;

    gmesh0->ResetReference();
    auto cmesh = new TPZMultiphysicsCompMesh(gmesh0);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    auto mat = new TPZMatPoissonCS(EDomain, dim);
    cmesh->InsertMaterialObject(mat);

    //Boundary Conditions
    TPZFMatrix<STATE> val1(1,1,1.);
    TPZManVector<STATE> val2(1,0);
    auto * BCond0 = mat->CreateBC(mat, EPointLeft0, 0, val1, val2);
    // auto * BCond1 = mat->CreateBC(mat, EPointRight1, 0, val1, val2);
    BCond0->SetForcingFunctionBC(exactSol);
    // BCond1->SetForcingFunctionBC(exactSol);
    cmesh->InsertMaterialObject(BCond0);
    // cmesh->InsertMaterialObject(BCond1);
    
    TPZManVector<int> active(1,1);
    active[0]=1;
    // active[1]=0;
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    cmesh->BuildMultiphysicsSpace(active, meshvector);
    cmesh->LoadReferences();
    cmesh->CleanUpUnconnectedNodes(); 

    // Prints Multiphysics mesh
    std::ofstream myfile("MultiPhysicsMesh.txt");
    cmesh->Print(myfile);

    //Analysis and solve
    TPZLinearAnalysis an(cmesh,false);
    //sets number of threads to be used by the solver
    constexpr int nThreads{0};
    TPZSSpStructMatrix<STATE,TPZStructMatrixOR<STATE>> matskl(cmesh);   
    matskl.SetNumThreads(nThreads);  
    an.SetStructuralMatrix(matskl);

    ///Setting a direct solver
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);//ELU //ECholesky // ELDLt
    an.SetSolver(step);
    //assembles the system
    an.Assemble();

    ///solves the system
    an.Solve();


    
    return 0;
}



void CreateGeoMeshes(TPZGeoMesh *gmesh0, TPZGeoMesh *gmesh1){
    
    //Create GeoMesh 1
    const REAL xMin0=0.;
    const REAL xMax0=2.;
    const int nels0 = 4;
    TPZVec<int> matids0 = {EDomain, EPointLeft0, EPointRight0}; 
    gmesh0 = TPZGeoMeshTools::CreateGeoMesh1D(xMin0,xMax0,nels0,matids0,true);

    std::ofstream out0("gmesh0.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh0, out0);

    //Create GeoMesh 1
    const REAL xMin1=1.;
    const REAL xMax1=3.;
    const int nels1 = 4;
    TPZVec<int> matids1 = {EDomain, EPointLeft1, EPointRight1}; 
    gmesh1 = TPZGeoMeshTools::CreateGeoMesh1D(xMin1,xMax1,nels1,matids1,true);

    std::ofstream out1("gmesh1.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh1, out1);

}



