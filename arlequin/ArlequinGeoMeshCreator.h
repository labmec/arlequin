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

class ArlequinGeoMeshCreator{

private:
//Material ids
    int EGluing, EGlobal, ELocal;

public:
    TPZGeoMesh* CreateFineGeoMesh(TPZGeoMesh* gmesh, std::set<int> &overlapMatId);

    void CopyGeoElCoarseToArlequinGMesh(TPZGeoEl *gel, TPZGeoMesh *gmesh, std::map<int64_t,int64_t> &gl2auxNdIdx, std::map<int64_t,int64_t> &gl2auxElIdx, std::map<int64_t,int64_t> &node2model);
    void CopyGeoElFineToArlequinGMesh(TPZGeoEl *gel, TPZGeoMesh *gmesh, std::map<int64_t,int64_t> &lc2auxNdIdx, std::map<int64_t,int64_t> &lc2auxElIdx, std::map<int64_t,int64_t> &node2model);
    void CopyGeoElGluingToArlequinGMesh(TPZGeoMesh *gmeshAux, TPZGeoMesh *gmeshArlequin, std::map<int64_t,int64_t> &gl2auxNdIdx, std::map<int64_t,int64_t> &gl2auxElIdx, std::map<int64_t,int64_t> &lc2auxNdIdx, std::map<int64_t,int64_t> &lc2auxElIdx, std::map<int64_t,int64_t> &node2model);

    TPZGeoMesh* AssociateModels(TPZGeoMesh *gmeshCoarse,TPZGeoMesh *gMeshFine, std::set<int> &overlapMatId);//Pode ser alterado no futuro para um man vector. 
    void SetMaterialIds(int gluing, int global, int local){
        EGluing = gluing;
        EGlobal = global;
        ELocal = local;
    };

};



#endif