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
#include "TPZLagrangeMultiplierCS.h"
#include <pzskylstrmatrix.h> //symmetric skyline matrix storage
#include "arlequin_config.h"
#include "TPZRefPattern.h"
#include "TPZGeoElement.h"
#include "TPZRefLinear.h"
#include "pzgeoel.h"

enum EMatId {ENone, EPoint0, EPoint1, EPoint2, EPoint3, EGlobal, EGluing, ELocal};

/**
   @brief Reads the test mesh from gmsh
   @param[in] file_name the .msh mesh file.
*/
TPZGeoMesh* ReadMeshFromGmsh(std::string file_name);

TPZGeoMesh* CreateFineGeoMesh(TPZGeoMesh* gmesh, std::set<int> &overlapMatId);

TPZCompMesh* CreateModelCMesh(TPZGeoMesh* gmesh,std::set<int> &allMat);

TPZCompMesh* CreateGluingCMesh(TPZGeoMesh* gmesh);

TPZMultiphysicsCompMesh* CreateArlequinModel(TPZGeoMesh* gmesh, TPZVec<TPZCompMesh *> &meshvector);
void PrintGeoMesh(TPZGeoMesh *gmesh);

void WriteGeoPoint(std::ofstream &file, TPZGeoNode &point, REAL refin=0.25);
void WriteGeoLine(std::ofstream &file, TPZGeoElSide &side);
void WritePhysicalLine(std::ofstream &file, TPZGeoElSide &side);
TPZGeoMesh* AssociateModels(TPZGeoMesh *gmeshCoarse,TPZGeoMesh *gMeshFine, std::set<int> &overlapMatId);//Pode ser alterado no futuro para um man vector. 

void CopyGeoElCoarseToArlequinGMesh(TPZGeoEl *gel, TPZGeoMesh *gmesh, std::map<int64_t,int64_t> &gl2auxNdIdx, std::map<int64_t,int64_t> &gl2auxElIdx, std::map<int64_t,int64_t> &node2model);
void CopyGeoElFineToArlequinGMesh(TPZGeoEl *gel, TPZGeoMesh *gmesh, std::map<int64_t,int64_t> &lc2auxNdIdx, std::map<int64_t,int64_t> &lc2auxElIdx, std::map<int64_t,int64_t> &node2model);
void CopyGeoElGluingToArlequinGMesh(TPZGeoMesh *gmeshAux, TPZGeoMesh *gmeshArlequin, std::map<int64_t,int64_t> &gl2auxNdIdx, std::map<int64_t,int64_t> &gl2auxElIdx, std::map<int64_t,int64_t> &lc2auxNdIdx, std::map<int64_t,int64_t> &lc2auxElIdx, std::map<int64_t,int64_t> &node2model);

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
    const int dim{1};
    const int pOrder{1};

#ifdef PZ_LOG
TPZLogger::InitializePZLOG();
#endif

    //Create Geometric meshes
    TPZGeoMesh* gmeshCoarse = ReadMeshFromGmsh(string(MESHDIR)+"bar.msh");
    // PrintGeoMesh(gmeshCoarse);
    std::set<int> overlapMatId = {EGluing};
    TPZGeoMesh* gmeshFine = CreateFineGeoMesh(gmeshCoarse, overlapMatId);
    // PrintGeoMesh(gmeshFine);
    TPZGeoMesh* gmeshArlequin = AssociateModels(gmeshCoarse,gmeshFine,overlapMatId);
    delete gmeshCoarse, gmeshFine;
    PrintGeoMesh(gmeshArlequin);

    TPZVec<TPZCompMesh *> meshvector(3); 

    // //Cria uma malha computacional apenas com os material ids dos elementos globais.
    // std::set<int> allMat={1,5,6};
    // meshvector[0] = CreateModelCMesh(gmeshArlequin,allMat);
    // std::string cmesh0 = "CMesh0.txt";
    // std::ofstream myfile0(cmesh0);
    // meshvector[0]->Print(myfile0);

    // //Cria uma malha computacional apenas com os material ids dos elementos locais.
    // std::set<int> allMat2={105};
    // meshvector[1] = CreateModelCMesh(gmeshArlequin,allMat2);
    // std::string cmesh1 = "CMesh1.txt";
    // std::ofstream myfile1(cmesh1);
    // meshvector[1]->Print(myfile1);

    // //Cria uma malha de pressão com o material id da zona de colagem
    // meshvector[2] = CreateGluingCMesh(gmeshArlequin);
    // std::string cmesh2 = "CMesh2.txt";
    // std::ofstream myfile2(cmesh2);
    // meshvector[2]->Print(myfile2);

    // //Cria a malha multifísica com as 3 malhas usando interface CompEl.
    // TPZMultiphysicsCompMesh *cmesh = CreateArlequinModel(gmeshArlequin,meshvector);
    // std::string cmesh3 = "CMeshArlequin.txt";
    // std::ofstream myfile3(cmesh3);
    // cmesh->Print(myfile3);

    //Create analysis environment
    // TPZLinearAnalysis an(meshvector[0],false);
    // TPZSkylineStructMatrix<REAL> matskl(cmesh);
    // an.SetStructuralMatrix(matskl);
    // TPZStepSolver<STATE> step;
    // step.SetDirect(ELDLt);
    // an.SetSolver(step);
    // an.Run();
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
        stringtoint[0]["Point0"] = EPoint0;
        stringtoint[0]["Point1"] = EPoint1;
        stringtoint[0]["Point2"] = EPoint2;
        stringtoint[0]["Point3"] = EPoint3;
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


TPZGeoMesh* CreateFineGeoMesh(TPZGeoMesh* gmesh, std::set<int> &overlapMatId)
{
    //GEO file
    std::string geo_name = string(MESHDIR)+"fine.geo";
    std::ofstream geofile(geo_name.c_str());

    std::set<int> nodeIds;
    std::set<int> lineIds;

    for (auto gel:gmesh->ElementVec()){
        if (gel->Dimension() != gmesh->Dimension()) continue;
        
        if (overlapMatId.find(gel->MaterialId()) != overlapMatId.end()){
            
            int nSides = gel->NSides();
            int nSideNodes = gel->NSideNodes(nSides-1);
            for (int inode = 0; inode < nSideNodes; inode++)
            {
                int nodeID = gel->Node(inode).Id();
                if (nodeIds.find(nodeID) == nodeIds.end()){
                    nodeIds.insert(nodeID);
                    WriteGeoPoint(geofile,gel->Node(inode));
                } 
            }
            TPZGeoElSide side(gel,nSides-1);
            WriteGeoLine(geofile,side);
            WritePhysicalLine(geofile,side);
            lineIds.insert(side.Id()+100);
        }
    }

    geofile.close();

    std::string mshfile = string(MESHDIR)+"fine.msh";
	std::string cmd = "gmsh -1 " + geo_name + " -o " + mshfile;
    system(cmd.c_str());

    //read mesh from gmsh
    TPZGeoMesh *gmeshFine;
    gmeshFine = new TPZGeoMesh();
    {
        TPZGmshReader reader;
        // essa interface permite voce mapear os nomes dos physical groups para
        // o matid que voce mesmo escolher
        TPZManVector<std::map<std::string,int>,4> stringtoint(4);
        for(auto line : lineIds) {
            std::string s = std::to_string(line);
            stringtoint[1][s] = line;
        }    

        reader.SetDimNamePhysical(stringtoint);
        reader.GeometricGmshMesh(mshfile,gmeshFine);
    }
    return gmeshFine;
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

    TPZNullMaterial<> *mat = new TPZNullMaterial<>(105);
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

    auto mat = new TPZMatPoissonCS(5,1);
    cmesh->InsertMaterialObject(mat);
    // auto mat2 = new TPZMatPoissonCS(105,1);
    // cmesh->InsertMaterialObject(mat2);
    // auto mat3 = new TPZMatPoissonCS(7,1);
    // cmesh->InsertMaterialObject(mat3);
    // auto mat4 = new TPZMatPoissonCS(10,1);
    // cmesh->InsertMaterialObject(mat4);


    TPZManVector<int> active(meshvector.size(),1);
    // active[2]=0;
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    cmesh->BuildMultiphysicsSpace(active, meshvector);
    cmesh->LoadReferences();
    cmesh->CleanUpUnconnectedNodes(); 


    // auto mat5 = new TPZLagrangeMultiplierCS<STATE>(105, 1);
    // cmesh->InsertMaterialObject(mat5);

    // cmesh->LoadReferences();
    
    // for (auto gel : gmesh->ElementVec())
    // {
    //     if (!gel || gel->MaterialId() != 10) continue;
    //     auto nsides = gel->NSides();
    //     TPZGeoElSide gelside(gel,nsides-1);
    //     // // here I generalized - an interface is created whenever a wrap element exists
    //     auto gelsidepr = gelside.HasNeighbour(6);
    //     if (!gelsidepr)
    //     {
    //         DebugStop();
    //     }

    //     TPZCompElSide celside = gelside.Reference();
    //     TPZCompElSide celneigh = gelsidepr.Reference();
    //     if (!celside || !celneigh) {
    //         DebugStop();
    //     }

    //     TPZGeoEl *gelIntface = gel->Neighbour(2).Element();
    //     if (gelIntface->MaterialId() != 11) DebugStop();
        
    //     // std::cout << "WRAP " << gel->Index() << ", normal left = " << celneigh.Element()->Reference()->NormalOrientation(6)
    //     //                                      << ", normal right = " << celside.Element()->Reference()->NormalOrientation(6)
    //     //                                      << std::endl;

    //     // Creates Multiphysics Interface element
    //     int64_t index;
    //     TPZMultiphysicsInterfaceElement *intf = new TPZMultiphysicsInterfaceElement(*cmesh,gelIntface,celneigh,celside);
    // }

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



void WriteGeoPoint(std::ofstream &file, TPZGeoNode &point, REAL refin){
    // std::cout << "WRITING POINT " << point.Id() << std::endl;
    file << "Point(" << point.Id() << ") = {" << point.Coord(0) << ", " 
                                                << point.Coord(1) << ", " 
                                                << point.Coord(2) << ", " << refin << "};\n\n"; 
}
void WriteGeoLine(std::ofstream &file, TPZGeoElSide &side){
    // std::cout << "WRITING LINE " << side.Id() << std::endl;
    file << "Line(" << side.Id() << ") = {" << side.SideNodeIndex(0) << ", " << side.SideNodeIndex(1) << "};\n\n"; 
}
void WritePhysicalLine(std::ofstream &file, TPZGeoElSide &side){
    std::cout << "WRITING Physical LINE " << side.Id() << std::endl;
    file << "Physical Curve(\"" << 100+side.Id() << "\") = {" << side.Id() << "};\n\n"; 
}


TPZGeoMesh* AssociateModels(TPZGeoMesh *gmeshCoarse,TPZGeoMesh *gmeshFine, std::set<int> &overlapMatId){

    std::map<int64_t,int64_t> gl2auxNdIdx;//global to aux node index
    std::map<int64_t,int64_t> gl2auxElIdx;//global to aux element index
    std::map<int64_t,int64_t> lc2auxNdIdx;//local to aux node index
    std::map<int64_t,int64_t> lc2auxElIdx;//local to aux element index
    std::map<int64_t,int64_t> node2model;//indicates to which model the node belongs, 0-coarse, 1-fine, 2-gluing

    TPZGeoMesh *gmeshArlequin;
    gmeshArlequin = new TPZGeoMesh();
    //Loop over all the fine model geometric elements
    for (auto gelFine:gmeshFine->ElementVec())
    {
        CopyGeoElFineToArlequinGMesh(gelFine,gmeshArlequin,lc2auxNdIdx,lc2auxElIdx,node2model);
    }

    //Loop over all the coarse model geometric elements
    for (auto gelCoarse:gmeshCoarse->ElementVec())
    {
        if (gelCoarse->Dimension() != gmeshCoarse->Dimension()) 
        {
            // Just copy the coarse gelEl to the Arlequin geomesh
            CopyGeoElCoarseToArlequinGMesh(gelCoarse,gmeshArlequin,gl2auxNdIdx,gl2auxElIdx,node2model);

        } else {
            if (overlapMatId.find(gelCoarse->MaterialId()) == overlapMatId.end())
            {
                // Just copy the coarse gelEl to the Arlequin geomesh
                CopyGeoElCoarseToArlequinGMesh(gelCoarse,gmeshArlequin,gl2auxNdIdx,gl2auxElIdx,node2model);

            } else {
                // continue;
                // The element is in the overlapping zone, so we need to create the refinement pattern

                TPZGeoMesh *gmeshAux;
                gmeshAux = new TPZGeoMesh();
                
                int count = 0;
                int countNode = 0;
                int ncorner = gelCoarse->NCornerNodes();
                
                //Global mesh information
                gmeshAux->ElementVec().Resize(1);
                gmeshAux->NodeVec().Resize(ncorner);

                gl2auxElIdx[gelCoarse->Index()] = count;
                count++;

                for (int icorn = 0; icorn < ncorner; icorn++)
                {
                    auto id = gmeshAux->CreateUniqueNodeId();
                    gmeshAux->NodeVec()[countNode] = gelCoarse->Node(icorn);
                    gmeshAux->NodeVec()[countNode].SetNodeId(id);
                    gl2auxNdIdx[gelCoarse->NodeIndex(icorn)] = countNode;
                    // node2model[countNode] = 0;
                    countNode++;            
                }
                gelCoarse->ClonePatchEl(*gmeshAux,gl2auxNdIdx,gl2auxElIdx);
                //

                for (auto gelFine:gmeshFine->ElementVec())
                {
                    if(gelCoarse->Id() != gelFine->MaterialId()-100) continue;
                    
                    gmeshAux->ElementVec().Resize(count+1);
                    
                    lc2auxElIdx[gelFine->Index()] = count;
                    count++;

                    for (int icorner = 0; icorner < ncorner; icorner++){
                        auto x = gelFine->Node(icorner).Coord(0);
                        auto y = gelFine->Node(icorner).Coord(1);
                        bool createNode = true;
                        auto NNodes = gmeshAux->NNodes();
                        
                        //Avoid duplicate nodes - se quiser que os nós das extremidades sejam duplicados é so comecar esse loop de ncorner ao inves de zero
                        for (int k = ncorner; k < NNodes; k++){
                            auto xMesh = gmeshAux->NodeVec()[k].Coord(0);
                            auto yMesh = gmeshAux->NodeVec()[k].Coord(1);
                            REAL dist = sqrt((x-xMesh)*(x-xMesh) + (y-yMesh)*(y-yMesh));
                            if (dist <= 1.e-5){
                                createNode = false;
                                break;
                            }
                        }
                        if (createNode){
                            gmeshAux->NodeVec().Resize(NNodes+1);
                            auto id = gmeshAux->CreateUniqueNodeId();
                            gmeshAux->NodeVec()[countNode] = gelFine->Node(icorner);
                            gmeshAux->NodeVec()[countNode].SetNodeId(id);
                            lc2auxNdIdx[gelFine->NodeIndex(icorner)] = countNode;
                            countNode++;
                        }
                    }
                    gelFine->ClonePatchEl(*gmeshAux,lc2auxNdIdx,lc2auxElIdx);
                    std::cout << "gelFine " << gelFine->Index() << " associated to gelCoarse " << gelCoarse->Index() << std::endl; 
                }

                gmeshAux->ResetConnectivities();
                gmeshAux->BuildConnectivity();

                TPZRefPattern ref(*gmeshAux);
                gmeshAux->ElementVec()[0]->SetRefPattern(&ref);

                for (int i=0; i<ref.NSubElements(); i++){
                    auto id = gmeshAux->ElementVec()[i+1]->Index()-1;
                    auto &gelAux = gmeshAux->ElementVec()[i+1];
                    // gmeshAux->ElementVec()[0]->SetSubElement(id,gelAux);
                    // gelAux->SetFather(gmeshAux->ElementVec()[0]);
                }

                std::string file = "gmeshAux.txt";
                std::ofstream myfile2(file);
                gmeshAux->Print(myfile2);
                
                std::string fileref = "refPattern.txt";
                std::ofstream myfileref(fileref);
                ref.Print(myfileref);
                // PrintGeoMesh(gmeshAux);

                CopyGeoElGluingToArlequinGMesh(gmeshAux,gmeshArlequin,gl2auxNdIdx,gl2auxElIdx,lc2auxNdIdx,lc2auxElIdx,node2model);

                // gmeshArlequin = gmeshAux;
            }//Else GelCoarse material id in overlapping zone
        }//Else GelCoarse Dimension
    }

    gmeshArlequin->ResetConnectivities();
    gmeshArlequin->BuildConnectivity();

    for (int i = 0; i < gmeshArlequin->NodeVec().NElements(); i++)
    {
        std::cout << "Node " << i<< " , id = " << gmeshArlequin->NodeVec()[i].Id() << ", model = " << node2model[i] << " , Coord = " <<  gmeshArlequin->NodeVec()[i].Coord(0) << " " <<  gmeshArlequin->NodeVec()[i].Coord(1) << std::endl;
    }
    
    return gmeshArlequin;
}


void CopyGeoElCoarseToArlequinGMesh(TPZGeoEl *gelCoarse, TPZGeoMesh *gmesh, std::map<int64_t,int64_t> &gl2auxNdIdx, std::map<int64_t,int64_t> &gl2auxElIdx, std::map<int64_t,int64_t> &node2model){

    // if (gelCoarse->MaterialId() == ELocal) return;
    auto nElements = gmesh->NElements();
    int ncorner = gelCoarse->NCornerNodes();

    gmesh->ElementVec().Resize(nElements+1);
    
    auto elId = gmesh->CreateUniqueElementId();
    gl2auxElIdx[gelCoarse->Index()] = elId;    

    for (int icorner = 0; icorner < ncorner; icorner++){
        auto x = gelCoarse->Node(icorner).Coord(0);
        auto y = gelCoarse->Node(icorner).Coord(1);
        bool createNode = true;
        auto NNodes = gmesh->NNodes();
        
        //Avoid duplicate nodes - se quiser que os nós das extremidades sejam duplicados é so comecar esse loop de ncorner ao inves de zero
        for (int k = 0; k < NNodes; k++){
            if (node2model[k] != 0) continue;
            auto xMesh = gmesh->NodeVec()[k].Coord(0);
            auto yMesh = gmesh->NodeVec()[k].Coord(1);
            REAL dist = sqrt((x-xMesh)*(x-xMesh) + (y-yMesh)*(y-yMesh));
            if (dist <= 1.e-5){
                createNode = false;
                // gelCoarse->SetNodeIndex(icorner,gmesh->NodeVec()[k].Id());
                gl2auxNdIdx[gelCoarse->NodeIndex(icorner)] = gmesh->NodeVec()[k].Id();
                break;
            }
        }
        if (createNode){
            gmesh->NodeVec().Resize(NNodes+1);
            auto nodeId = gmesh->CreateUniqueNodeId();
            gmesh->NodeVec()[nodeId] = gelCoarse->Node(icorner);
            gmesh->NodeVec()[nodeId].SetNodeId(nodeId);
            gl2auxNdIdx[gelCoarse->NodeIndex(icorner)] = nodeId;
            node2model[nodeId] = 0;
        }
    }
    gelCoarse->ClonePatchEl(*gmesh,gl2auxNdIdx,gl2auxElIdx);
    gmesh->ElementVec()[nElements]->SetIndex(elId);
    gmesh->ElementVec()[nElements]->SetMaterialId(EGlobal);

}

void CopyGeoElFineToArlequinGMesh(TPZGeoEl *gelFine, TPZGeoMesh *gmesh, std::map<int64_t,int64_t> &lc2auxNdIdx, std::map<int64_t,int64_t> &lc2auxElIdx, std::map<int64_t,int64_t> &node2model){

    auto nElements = gmesh->NElements();
    int ncorner = gelFine->NCornerNodes();

    gmesh->ElementVec().Resize(nElements+1);
    
    auto elId = gmesh->CreateUniqueElementId();
    lc2auxElIdx[gelFine->Index()] = elId;    

    for (int icorner = 0; icorner < ncorner; icorner++){
        auto x = gelFine->Node(icorner).Coord(0);
        auto y = gelFine->Node(icorner).Coord(1);
        bool createNode = true;
        auto NNodes = gmesh->NNodes();
        
        //Avoid duplicate nodes - se quiser que os nós das extremidades sejam duplicados é so comecar esse loop de ncorner ao inves de zero
        for (int k = 0; k < NNodes; k++){
            auto xMesh = gmesh->NodeVec()[k].Coord(0);
            auto yMesh = gmesh->NodeVec()[k].Coord(1);
            REAL dist = sqrt((x-xMesh)*(x-xMesh) + (y-yMesh)*(y-yMesh));
            if (dist <= 1.e-5 && node2model[k] == 1){
                createNode = false;
                break;
            }
        }
        if (createNode){
            gmesh->NodeVec().Resize(NNodes+1);
            auto nodeId = gmesh->CreateUniqueNodeId();
            gmesh->NodeVec()[nodeId] = gelFine->Node(icorner);
            gmesh->NodeVec()[nodeId].SetNodeId(nodeId);
            lc2auxNdIdx[gelFine->NodeIndex(icorner)] = nodeId;
            node2model[nodeId] = 1;
        }
    }
    gelFine->ClonePatchEl(*gmesh,lc2auxNdIdx,lc2auxElIdx);
    gmesh->ElementVec()[nElements]->SetIndex(elId);
    gmesh->ElementVec()[nElements]->SetMaterialId(ELocal);

}


void CopyGeoElGluingToArlequinGMesh(TPZGeoMesh *gmeshAux, TPZGeoMesh *gmeshArlequin, 
                                    std::map<int64_t,int64_t> &gl2auxNdIdx, std::map<int64_t,int64_t> &gl2auxElIdx, 
                                    std::map<int64_t,int64_t> &lc2auxNdIdx, std::map<int64_t,int64_t> &lc2auxElIdx, 
                                    std::map<int64_t,int64_t> &node2model){

    int nElemAux = gmeshAux->NElements();

    auto nElements = gmeshArlequin->NElements();
    gmeshArlequin->ElementVec().Resize(nElements+nElemAux);

    for (int i = 0; i < gmeshArlequin->NodeVec().NElements(); i++)
    {
        std::cout << "Node " << i<< " , id = " << gmeshArlequin->NodeVec()[i].Id() << ", model = " << node2model[i] << " , Coord = " <<  gmeshArlequin->NodeVec()[i].Coord(0) << " " <<  gmeshArlequin->NodeVec()[i].Coord(1) << std::endl;
    }

    for (int iEl = 0; iEl < nElemAux; iEl++)
    {
        auto gel = gmeshAux->ElementVec()[iEl];
        auto elId = gmeshArlequin->CreateUniqueElementId();
        int ncorner = gel->NCornerNodes();
        
        for (int icorner = 0; icorner < ncorner; icorner++){
            auto x = gel->Node(icorner).Coord(0);
            auto y = gel->Node(icorner).Coord(1);
            bool createNode = true;
            auto NNodes = gmeshArlequin->NNodes();
            // Check for the created nodes
            int domain = iEl == 0 ? 0 : 2;

            //Avoid duplicate nodes - se quiser que os nós das extremidades sejam duplicados é so comecar esse loop de ncorner ao inves de zero
            for (int k = 0; k < NNodes; k++){
                if (domain != node2model[k]) continue;
                auto xMesh = gmeshArlequin->NodeVec()[k].Coord(0);
                auto yMesh = gmeshArlequin->NodeVec()[k].Coord(1);
                REAL dist = sqrt((x-xMesh)*(x-xMesh) + (y-yMesh)*(y-yMesh));
                if (dist <= 1.e-5){
                    createNode = false;
                    gel->SetNodeIndex(icorner,gmeshArlequin->NodeVec()[k].Id());
                    break;
                }
            }
            if (createNode){
                gmeshArlequin->NodeVec().Resize(NNodes+1);
                auto nodeId = gmeshArlequin->CreateUniqueNodeId();
                gmeshArlequin->NodeVec()[nodeId] = gel->Node(icorner);
                gmeshArlequin->NodeVec()[nodeId].SetNodeId(nodeId);
                gel->SetNodeIndex(icorner,nodeId);
                node2model[nodeId] = domain;
            }
        }
        gel->SetIndex(elId);
        gmeshArlequin->ElementVec()[nElements+iEl] = gel->Clone(*gmeshArlequin);
        if (iEl == 0){
            gmeshArlequin->ElementVec()[nElements+iEl]->SetMaterialId(EGlobal);
        } else {
            gmeshArlequin->ElementVec()[nElements+iEl]->SetMaterialId(EGluing);
        }
    }

}