#include "ArlequinGeoMeshCreator.h"
#include <string>
#include "GmshFileWriter.h"
#include <TPZGmshReader.h>
#include "TPZRefPattern.h"

TPZGeoMesh* ArlequinGeoMeshCreator::CreateFineGeoMesh(TPZGeoMesh* gmesh, std::set<int> &overlapMatId)
{
    //GEO file
    std::string geo_name = std::string(MESHDIR)+"fine.geo";
    std::ofstream geofile(geo_name.c_str());

    std::set<int> nodeIds;
    std::set<int> lineIds;
    GmshFileWriter writer;

    for (auto gel:gmesh->ElementVec()){
        if (gel->Dimension() != gmesh->Dimension()) continue;
        
        if (overlapMatId.find(gel->MaterialId()) != overlapMatId.end()){
            
            int nSides = gel->NSides();            
            int nSideNodes = gel->NSideNodes(nSides-1);
            int nCorner = gel->NCornerNodes();
            for (int inode = 0; inode < nSideNodes; inode++)
            {
                int nodeID = gel->Node(inode).Id();
                if (nodeIds.find(nodeID) == nodeIds.end()){
                    nodeIds.insert(nodeID);
                    writer.WriteGeoPoint(geofile,gel->Node(inode),0.25);
                } 
            }
            for (int iSide = nCorner; iSide < nSides; iSide++)
            {
                TPZGeoElSide side(gel,iSide);
                
                if (side.Dimension() == 1){
                    
                    if (gmesh->Dimension() == 1) {
                        writer.WriteGeoLine(geofile,side,true);
                        writer.WritePhysicalLine(geofile,side);
                        lineIds.insert(side.Id()+100);
                    } else {
                        writer.WriteGeoLine(geofile,side,false);
                    }
                } else if (side.Dimension() == 2) {
                    int nfacets = nSides - nCorner - 1;
                    if (gmesh->Dimension() == 2) {
                        auto sideType = side.Element()->Type(iSide);
                        writer.WriteGeoSurface(geofile,side, nfacets, sideType, true);
                        writer.WritePhysicalSurface(geofile,side);
                    }
                }
            }
            
        }
    }

    geofile.close();

    std::string mshfile = std::string(MESHDIR)+"fine.msh";
	std::string cmd = "gmsh -" + std::to_string(gmesh->Dimension()) + " " + geo_name + " -o " + mshfile;
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


void ArlequinGeoMeshCreator::CopyGeoElCoarseToArlequinGMesh(TPZGeoEl *gelCoarse, TPZGeoMesh *gmesh, std::map<int64_t,int64_t> &gl2auxNdIdx, std::map<int64_t,int64_t> &gl2auxElIdx, std::map<int64_t,int64_t> &node2model){

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
    gmesh->ElementVec()[nElements]->SetId(elId);
    gmesh->ElementVec()[nElements]->SetMaterialId(EGlobal);

}

void ArlequinGeoMeshCreator::CopyGeoElFineToArlequinGMesh(TPZGeoEl *gelFine, TPZGeoMesh *gmesh, std::map<int64_t,int64_t> &lc2auxNdIdx, std::map<int64_t,int64_t> &lc2auxElIdx, std::map<int64_t,int64_t> &node2model){

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
    gmesh->ElementVec()[nElements]->SetId(elId);
    gmesh->ElementVec()[nElements]->SetMaterialId(ELocal);

}


void ArlequinGeoMeshCreator::CopyGeoElGluingToArlequinGMesh(TPZGeoMesh *gmeshAux, TPZGeoMesh *gmeshArlequin, 
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
        gel->SetId(elId);
        gmeshArlequin->ElementVec()[nElements+iEl] = gel->Clone(*gmeshArlequin);
        if (iEl == 0){
            gmeshArlequin->ElementVec()[nElements+iEl]->SetMaterialId(EGlobal);
        } else {
            gmeshArlequin->ElementVec()[nElements+iEl]->SetMaterialId(EGluing);
        }
    }

// "Mudar onde tem EGLuing para EGlobal e atribuir um material novo, parecido com o MixedDarcyFlow, só que ao inves de usar 
// Deformed Directions, tem que usar phi."


}



TPZGeoMesh* ArlequinGeoMeshCreator::AssociateModels(TPZGeoMesh *gmeshCoarse,TPZGeoMesh *gmeshFine, std::set<int> &overlapMatId){

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


