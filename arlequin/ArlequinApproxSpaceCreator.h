/**
 * @file ArlequinGeoMeshCreator.h
 * @brief Contains the ArlequinGeoMeshCreator which creates a geometric mesh based on an original coarse scale geoMesh
 * with different material ids for each zone (Coarse, Fine and Gluing)
 */

#ifndef ArlequinGeoMeshCreator_H
#define ArlequinGeoMeshCreator_H


#include <pzgmesh.h> //for TPZGeoMesh
#include "TPZMultiphysicsCompMesh.h"
#include "arlequin_config.h"

enum EArlequinMatId {ENone, EGlobal, EGluing, ELocal, EWrapGlob, EWrapLoc, EInter, EDirichlet1, EDirichlet2, ENeumann, EDirichlet3};

class ArlequinApproxSpaceCreator{

private:
//Material ids
    TPZGeoMesh *fGeoMesh;

    int fOrder;

    int fDimension;

public:
    ArlequinApproxSpaceCreator(TPZGeoMesh *gmesh, int order);


    void PrintGeoMesh(TPZGeoMesh *gmesh);

    void CreateGeoElements();

    void RefineLocalModel();

    TPZGeoMesh *GeoMesh(){return fGeoMesh;}

    TPZCompMesh *GlobalCompMesh();

    TPZCompMesh *LocalCompMesh();

    TPZCompMesh *GluingCompMesh();

    TPZMultiphysicsCompMesh* ArlequinCompMesh();

    void ChangeInternalOrder(TPZCompMesh *cmesh, int pOrder);
};



#endif