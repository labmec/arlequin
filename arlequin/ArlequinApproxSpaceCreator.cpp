#include "ArlequinApproxSpaceCreator.h"

#include "TPZVTKGeoMesh.h"
#include "TPZNullMaterial.h"
#include "TPZNullMaterialCS.h"
#include "TPZMatPoissonCS.h"
#include "TPZLagrangeMultiplierCS.h"
#include "pzintel.h"
#include "pzgeoquad.h"
#include "pzshapequad.h"
#include "TPZCompElH1.h"
#include "TPZMultiphysicsInterfaceElArlequin.h"


ArlequinApproxSpaceCreator::ArlequinApproxSpaceCreator(TPZGeoMesh *gmesh, int order){
    fGeoMesh = gmesh;
    fDimension = fGeoMesh->Dimension();
    fOrder = order;
}

void ArlequinApproxSpaceCreator::CreateGeoElements(){
    
    int64_t nel = fGeoMesh->NElements();
    for (int64_t el = 0; el<nel; el++){
        TPZGeoEl* gel = fGeoMesh->ElementVec()[el];
        if (gel->MaterialId() == EGluing){
            TPZVec<int64_t> corner(gel->NCornerNodes());
            for (int i = 0; i < gel->NCornerNodes(); i++)
            {
                corner[i] = gel->SideNodeIndex(i,0);
                // std::cout << "Corner i " << i << " " << corner[i]<< std::endl;
            }
            int nsides = gel->NSides();
            auto newId = fGeoMesh->CreateUniqueElementId();
            fGeoMesh->CreateGeoElement(gel->Type(),corner,EGlobal,newId);
            // TPZGeoEl* gelGlobal = fGeoMesh->ElementVec()[newId];
            TPZGeoElSide geosideglo(gel,nsides-1);
            TPZGeoElBC gelbcWrapglo(geosideglo, EWrapGlob);
            TPZGeoElSide gelWrapSideglo(gelbcWrapglo.CreatedElement(),gelbcWrapglo.CreatedElement()->NSides()-1);
            TPZGeoElBC gelbcglo(gelWrapSideglo, EInter);
            TPZGeoElSide geosideloc(gelbcglo.CreatedElement(),nsides-1);
            TPZGeoElBC gelbcWraploc(geosideloc, EWrapLoc);
            TPZGeoElSide gelWrapSideloc(gelbcWraploc.CreatedElement(),gelbcWraploc.CreatedElement()->NSides()-1);
            TPZGeoElBC gelbcloc(gelWrapSideloc, ELocal);
        }
        // if (gel->MaterialId() == ELocal){
        //     auto newId = fGeoMesh->CreateUniqueElementId();
        //     TPZVec<int64_t> corner(gel->NCornerNodes());
        //     for (int i = 0; i < gel->NCornerNodes(); i++)
        //     {
        //         corner[i] = gel->SideNodeIndex(i,0);
        //         std::cout << "Corner i " << i << " " << corner[i]<< std::endl;
        //     }
            
        //     fGeoMesh->CreateGeoElement(gel->Type(),corner,EGlobal,newId);
        //     // gel->ResetReference();
        // }
    }
    fGeoMesh->BuildConnectivity();
    // PrintGeoMesh(fGeoMesh);
}

void ArlequinApproxSpaceCreator::RefineLocalModel(){

    int64_t nel = fGeoMesh->NElements();
    for (int i = 0; i < nel; i++){
        TPZGeoEl* gel = fGeoMesh->ElementVec()[i];
        //This makes the Lagrange Multiplier to be defined in the Fine mesh
        // if (gel->MaterialId() == EGluing || gel->MaterialId() == ELocal || gel->MaterialId() == EWrapLoc || gel->MaterialId() == EWrapGlob || gel->MaterialId() == EInter){
        //This makes the Lagrange Multiplier to be defined in the Coarse mesh
        if (gel->MaterialId() == ELocal){

            TPZManVector<TPZGeoEl*,10> children;
            fGeoMesh->ElementVec()[i]->Divide(children);
            //Divide BC
            for(int64_t el = 0; el<nel; el++) {
                TPZGeoEl *gel = fGeoMesh->Element(el);
                if(gel->Dimension() != fDimension-1) continue;
                if(gel->HasSubElement()) continue;
                TPZGeoElSide gelside(gel);
                TPZGeoElSide neighbour = gelside.Neighbour();
                if(neighbour.HasSubElement()) {
                    TPZManVector<TPZGeoEl*,10> children;
                    gel->Divide(children);
                }
            }
        }
    }
    // PrintGeoMesh(gmeshCoarse);

}

void ArlequinApproxSpaceCreator::PrintGeoMesh(TPZGeoMesh *gmesh){
    
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

    std::string gmesh2 = "GMesh.txt";
    std::ofstream myfile2(gmesh2);
    gmesh->Print(myfile2);

    //Prints gmesh mesh properties
    std::string vtk_name = "geoMesh.vtk";
    std::ofstream vtkfile(vtk_name.c_str());

    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, true);

}

TPZCompMesh * ArlequinApproxSpaceCreator::GlobalCompMesh(){

    TPZCompMesh *cmesh = new TPZCompMesh(fGeoMesh); 
    
    cmesh->SetDefaultOrder(fOrder);
    cmesh->SetDimModel(fDimension);

    auto mat = new TPZNullMaterial(EGlobal,fDimension);
    cmesh->InsertMaterialObject(mat);
    auto mat2 = new TPZNullMaterial(EWrapGlob,fDimension);
    cmesh->InsertMaterialObject(mat2);

    TPZFMatrix<STATE> val1(1,1,0.);
    TPZManVector<STATE> val2(1,0.);
    TPZBndCondT<STATE> *BCond1 = mat->CreateBC(mat, EDirichlet1, 0, val1, val2);
    cmesh->InsertMaterialObject(BCond1);


    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild();

    std::string cmesh1 = "CMeshGlobal.txt";
    std::ofstream myfile1(cmesh1);
    cmesh->Print(myfile1);

    return cmesh;
};

TPZCompMesh * ArlequinApproxSpaceCreator::LocalCompMesh(){

    TPZCompMesh *cmesh = new TPZCompMesh(fGeoMesh); 
    
    cmesh->SetDefaultOrder(fOrder);
    cmesh->SetDimModel(fDimension);

    auto mat = new TPZNullMaterial(ELocal,fDimension);
    cmesh->InsertMaterialObject(mat);
    auto mat2 = new TPZNullMaterial(EWrapLoc,fDimension);
    cmesh->InsertMaterialObject(mat2);

    TPZFMatrix<STATE> val1(1,1,0.);
    TPZManVector<STATE> val2(1,1.);
    TPZBndCondT<STATE> *BCond1 = mat->CreateBC(mat, EDirichlet2, 0, val1, val2);
    cmesh->InsertMaterialObject(BCond1);

    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild();

    std::string cmesh1 = "CMeshLocal.txt";
    std::ofstream myfile1(cmesh1);
    cmesh->Print(myfile1);

    return cmesh;
};

TPZCompMesh * ArlequinApproxSpaceCreator::GluingCompMesh(){

    fGeoMesh->ResetReference();
    TPZCompMesh *cmesh = new TPZCompMesh(fGeoMesh);

    TPZNullMaterial<> *mat = new TPZNullMaterial<>(EGluing,fDimension);
    // mat->SetDimension();
    cmesh->InsertMaterialObject(mat);

    cmesh->SetDefaultOrder(fOrder);
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    cmesh->SetDimModel(fDimension);
    cmesh->AutoBuild();

    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i]; 
        newnod.SetLagrangeMultiplier(1);
    }

    // int nel = cmesh->NElements();
    // for(int i=0; i<nel; i++){
    //     TPZCompEl *cel = cmesh->ElementVec()[i];
    //     TPZCompElDisc *celdisc = dynamic_cast<TPZCompElDisc *>(cel);
    //     if(!celdisc) continue;
    //     celdisc->SetConstC(1.);
    //     celdisc->SetTrueUseQsiEta();
    //     // espera-se elemento de pressao apenas para o contorno
    //     auto aaa = celdisc->Reference()->Dimension();
    //     // if(celdisc && celdisc->Reference()->Dimension() == cmesh->Dimension())
    //     // {
    //     //     DebugStop();
    //     // }
    // }

    // ChangeInternalOrder(cmesh,2);

    return cmesh;

};


void ArlequinApproxSpaceCreator::ChangeInternalOrder(TPZCompMesh *cmesh, int pOrder) {

    int64_t nel = cmesh->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        if (!cel) continue;

        TPZGeoEl *gel = cel->Reference();
        if (gel->Dimension() != cmesh->Dimension()) {//Only elements with the same dimension of the mesh
            continue;
        }
        int nc = cel->NConnects();
        //Gets the volumetric connect
        int64_t conIndex = cel->ConnectIndex(nc-1);
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
        if (!intel) DebugStop();

        TPZConnect &c = cmesh->ConnectVec()[conIndex];
        int64_t index;
        //Sets the new connect order
        c.SetOrder(pOrder,index);

        //Gets connect information to update block size (stiffness matrix data structure)
        int64_t seqnum = c.SequenceNumber();
        int nvar = 1;
        TPZMaterial * mat = cel->Material();
        if (mat) nvar = mat->NStateVariables();
        int nshape = intel->NConnectShapeF(nc-1,pOrder);
        c.SetNShape(nshape);
        // c.SetNState(nvar);
        cmesh->Block().Set(seqnum, nvar * nshape);

        cel->SetIntegrationRule(2*pOrder);
    }
    cmesh->InitializeBlock();
}

TPZMultiphysicsCompMesh* ArlequinApproxSpaceCreator::ArlequinCompMesh(){

    TPZVec<TPZCompMesh *> meshvector(3); 
    meshvector[0] = GlobalCompMesh();
    meshvector[1] = LocalCompMesh();
    meshvector[2] = GluingCompMesh();

    auto cmesh = new TPZMultiphysicsCompMesh(fGeoMesh);
    cmesh->SetDefaultOrder(fOrder);
    cmesh->SetDimModel(fDimension);

    auto mat = new TPZMatPoissonCS(EGlobal,fDimension);
    cmesh->InsertMaterialObject(mat);
    TPZFMatrix<STATE> val1(1,1,0.);
    TPZManVector<STATE> val2(1,0.);
    TPZBndCondT<STATE> *BCond2 = mat->CreateBC(mat, EDirichlet1, 0, val1, val2);
    cmesh->InsertMaterialObject(BCond2);

    auto mat2 = new TPZMatPoissonCS(ELocal,fDimension);
    cmesh->InsertMaterialObject(mat2);
    val2[0] = 1.;
    TPZBndCondT<STATE> *BCond1 = mat2->CreateBC(mat2, EDirichlet2, 0, val1, val2);
    cmesh->InsertMaterialObject(BCond1);

    auto mat3 = new TPZNullMaterialCS<>(EWrapLoc);
    cmesh->InsertMaterialObject(mat3);
    auto mat5 = new TPZNullMaterialCS<>(EWrapGlob);
    cmesh->InsertMaterialObject(mat5);
    auto mat4 = new TPZNullMaterialCS<STATE>(EGluing);
    cmesh->InsertMaterialObject(mat4);

    TPZManVector<int> active(meshvector.size(),1);
    // active[2]=0;
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    cmesh->BuildMultiphysicsSpace(active, meshvector);
    cmesh->LoadReferences();

    auto mat6 = new TPZLagrangeMultiplierCS<STATE>(EInter, fDimension);
    cmesh->InsertMaterialObject(mat6);

    for (auto gel : fGeoMesh->ElementVec())
    {
        if (!gel || (gel->MaterialId() != EWrapLoc && gel->MaterialId() != EWrapGlob)) continue;
        auto nsides = gel->NSides();
        TPZGeoElSide gelside(gel,nsides-1);
        // here I generalized - an interface is created whenever a wrap element exists
        auto gelsidepr = gelside.HasNeighbour(EGluing);
        if (!gelsidepr) DebugStop();
        if (!gelsidepr.Reference()) continue;
        
        auto glu = gelsidepr.Element();
        TPZCompElSide celside = gelside.Reference();
        TPZCompElSide celneigh = gelsidepr.Reference();
        if (!celside || !celneigh) DebugStop();

        auto gelint = gelside.HasNeighbour(EInter);
        if (!gelint.Element()) DebugStop();
        
        // Creates Multiphysics Interface element
        int64_t index;
        TPZMultiphysicsInterfaceElementArlequin *intf = new TPZMultiphysicsInterfaceElementArlequin(*cmesh,gelint.Element(),celside,celneigh);
       
    }

    return cmesh;
}
