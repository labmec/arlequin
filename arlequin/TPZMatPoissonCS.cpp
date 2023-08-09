#include "TPZMatPoissonCS.h"
#include "TPZMaterialDataT.h"
#include "TPZBndCondT.h"
#include "pzaxestools.h"

template<class TVar>
TPZMatPoissonCS<TVar>::TPZMatPoissonCS(int id, int dim) :
    TPZRegisterClassId(&TPZMatPoissonCS::ClassId),
    TBase(id), fDim(dim), fSol(0)
{
}

template<class TVar>
TPZMaterial * TPZMatPoissonCS<TVar>::NewMaterial() const{
	return new TPZMatPoissonCS(*this);
}

template<class TVar>
void TPZMatPoissonCS<TVar>::Contribute(const TPZVec<TPZMaterialDataT<TVar>> &datavec,
                                       REAL weight,
                                       TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef){
	
	const auto nLoads = this->fNumLoadCases;
    TPZManVector<TVar,10> force(nLoads,0.); 
    if(this->HasForcingFunction()){
        this->fForcingFunction(datavec[0].x,force);
    }
    const auto &phiCoarse = datavec[0].phi;
    const auto &dphiCoarse = datavec[0].dphix;
    const auto nshapeCoarse = datavec[0].phi.Rows();
	
    // if (datavec[0].x[0] >= 1. && datavec[0].x[0] <= 2.){
    //     fScale = 0.5;
    // } else {
    //     fScale = 1.;
    // }
    
    // std::cout << "Xcoarse = " << datavec[0].x[0] << " , fScale = " << fScale << std::endl;
    // std::cout << "XFine   = " << datavec[1].x[0] << " , fScale = " << fScale << std::endl;
    for(int i = 0; i < nshapeCoarse; i++){
		for(int j = 0; j < nshapeCoarse; j++){
            STATE dphiIdphiJ = 0;
            for(int x = 0; x < fDim; x++){
                dphiIdphiJ += dphiCoarse.GetVal(x,i) * dphiCoarse.GetVal(x,j);
            }
            ek(i, j) += weight*fScale*dphiIdphiJ;
        }//forj
        for(auto l = 0; l < nLoads; l++)
            ef(i,l) += weight*fScale*phiCoarse.GetVal(i,0)*force[l];
    }//for i

    const auto &phiFine = datavec[1].phi;
    const auto &dphiFine = datavec[1].dphix;
    const auto nshapeFine = datavec[1].phi.Rows();
	// if (datavec[1].x[0] >= 1. && datavec[1].x[0] <= 2.){
    //     fScale = 0.5;
    // } else {
    //     fScale = 1.;
    // }
    for(int i = 0; i < nshapeFine; i++){
		for(int j = 0; j < nshapeFine; j++){
            STATE dphiIdphiJ = 0;
            for(int x = 0; x < fDim; x++){
                dphiIdphiJ += dphiFine.GetVal(x,i) * dphiFine.GetVal(x,j);
            }
            ek(i, j) += weight*fScale*dphiIdphiJ;
        }//forj
        for(auto l = 0; l < nLoads; l++)
            ef(i,l) += weight*fScale*phiFine.GetVal(i,0)*force[l];
    }//for i

    
    // const auto &phiFine = datavec[1].phi;
    // const auto &dphiFine = datavec[1].dphix;
    // const auto nshapeFine = datavec[1].phi.Rows();
	
    // for(int i = 0; i < nshapeFine; i++){
	// 	for(int j = 0; j < nshapeFine; j++){
    //         STATE dphiIdphiJ = 0;
    //         for(int x = 0; x < fDim; x++){
    //             dphiIdphiJ += dphiFine.GetVal(x,i) * dphiFine.GetVal(x,j);
    //         }
    //         ek(i, j) += weight*fScale*dphiIdphiJ;
    //     }//forj
    //     for(auto l = 0; l < nLoads; l++)
    //         ef(i,l) += weight*fScale*phiFine.GetVal(i,0)*force[l];
    // }//for i


}

template<class TVar>
void TPZMatPoissonCS<TVar>::ContributeBC(const TPZVec<TPZMaterialDataT<TVar>> &datavec,
                                         REAL weight,
                                         TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef,
                                         TPZBndCondT<TVar> &bc)
{
	
	const auto &phiGlo = datavec[0].phi;
    const auto &dphiGlo = datavec[0].dphix;
	const auto &phiLoc = datavec[1].phi;
    const auto &dphiLoc = datavec[1].dphix;
	const int phrG = phiGlo.Rows();
	const int phrL = phiLoc.Rows();
    constexpr int nvars = 1;
    const auto nloads = this->fNumLoadCases;
    const auto &bcNumLoads =
        dynamic_cast<TPZMatLoadCasesBC<TVar>&>(bc);

    TPZManVector<TVar,10> v2(nvars*nloads);
    TPZFNMatrix<30,TVar> v1(nvars,1);
	[&bc = std::as_const(bc),
     &bcNumLoads = std::as_const(bcNumLoads),
     &datavec = std::as_const(datavec),
     nvars,nloads]( TPZFMatrix<TVar> &v1, TPZVec<TVar> &v2) {
        if(bc.HasForcingFunctionBC()){
            bc.ForcingFunctionBC()(datavec[0].x,v2,v1);
        }else {
            for(auto l = 0; l < nloads; l++){
                const auto &val2 = bcNumLoads.GetBCRhsVal(l);
                for(auto i = 0; i < nvars; i++)
                    v2[nvars*l+i] = val2[i];
            }
            v1 = bc.Val1();
        }
    }(v1,v2);
    
	switch (bc.Type()){		
        // Dirichlet condition
    case 0 : {      
        for(auto iv = 0; iv < nvars; iv++){
            for(auto in = 0 ; in < phrG; in++) {
                for(auto l=0; l < nloads; l++)
                    ef(nvars*in+iv,l) +=
                        (TVar)TPZMaterial::fBigNumber * v2[nvars*l+iv] * (TVar)phiGlo.GetVal(in,0) * (TVar)weight;
                for (auto jn = 0 ; jn < phrG; jn++) {
                    ek(nvars*in+iv,nvars*jn+iv) += TPZMaterial::fBigNumber * phiGlo.GetVal(in,0) * phiGlo.GetVal(jn,0) * weight;
                }//jn
            }//in
        }//iv
        for(auto iv = 0; iv < nvars; iv++){
            for(auto in = 0 ; in < phrL; in++) {
                for(auto l=0; l < nloads; l++)
                    ef(nvars*in+iv,l) +=
                        (TVar)TPZMaterial::fBigNumber * v2[nvars*l+iv] * (TVar)phiLoc.GetVal(in,0) * (TVar)weight;
                for (auto jn = 0 ; jn < phrL; jn++) {
                    ek(nvars*in+iv,nvars*jn+iv) += TPZMaterial::fBigNumber * phiLoc.GetVal(in,0) * phiLoc.GetVal(jn,0) * weight;
                }//jn
            }//in
        }//iv
        break;
    }
		// Neumann condition
    case 1 : {
        for(auto iv = 0; iv < nvars; iv++){
            for(auto in = 0 ; in < phrG; in++) {
                for(auto l = 0; l < nloads; l++)
                    ef(nvars*in+iv,l) += v2[nvars*l+iv] * (TVar)fScale * (TVar)phiGlo.GetVal(in,0) * (TVar)weight;
            }//in
        }//iv
        for(auto iv = 0; iv < nvars; iv++){
            for(auto in = 0 ; in < phrL; in++) {
                for(auto l = 0; l < nloads; l++)
                    ef(nvars*in+iv,l) += v2[nvars*l+iv] * (TVar)fScale * (TVar)phiLoc.GetVal(in,0) * (TVar)weight;
            }//in
        }//iv
        break;
    }
        //Robin condition
    case 2 : {
        for(auto iv = 0; iv < nvars; iv++){
            for(auto in = 0 ; in < phrG; in++) {
                for(auto l=0; l < nloads; l++)
                    ef(nvars*in+iv,l) +=
                        (TVar)TPZMaterial::fBigNumber * v2[nvars*l+iv] * (TVar)phiGlo.GetVal(in,0) * (TVar)weight;
                for (auto jn = 0 ; jn < phrG; jn++) {
                    ek(nvars*in+iv,nvars*jn+iv) +=
                        TPZMaterial::fBigNumber * v1.GetVal(iv,0) * dphiGlo.GetVal(0,in) * dphiGlo.GetVal(0,jn) * weight;
                }//jn
            }//in
        }//iv
        for(auto iv = 0; iv < nvars; iv++){
            for(auto in = 0 ; in < phrL; in++) {
                for(auto l=0; l < nloads; l++)
                    ef(nvars*in+iv,l) +=
                        (TVar)TPZMaterial::fBigNumber * v2[nvars*l+iv] * (TVar)phiLoc.GetVal(in,0) * (TVar)weight;
                for (auto jn = 0 ; jn < phrL; jn++) {
                    ek(nvars*in+iv,nvars*jn+iv) +=
                        TPZMaterial::fBigNumber * v1.GetVal(iv,0) * dphiLoc.GetVal(0,in) * dphiLoc.GetVal(0,jn) * weight;
                }//jn
            }//in
        }//iv
        break;
    }		
    default:{
        std::cout << __PRETTY_FUNCTION__ << " at line " << __LINE__ << " not implemented\n";
    }
	}//switch
	
}
template<class TVar>
void TPZMatPoissonCS<TVar>::GetSolDimensions(uint64_t &u_len,
                                             uint64_t &du_row,
                                             uint64_t &du_col) const
{
    u_len=1;
    du_row=fDim;
    du_col=1;
}


template<class TVar>
int TPZMatPoissonCS<TVar>::VariableIndex(const std::string &name) const{
	if(!strcmp("Solution",name.c_str())) return ESolution;
    if(!strcmp("Derivative",name.c_str())) return EDerivative;
    // if(!strcmp("LagrangeMultiplier",name.c_str())) return ELagrangeMult;
	return TPZMaterial::VariableIndex(name);
}

template<class TVar>
int TPZMatPoissonCS<TVar>::NSolutionVariables(int var) const{
	if(var == ESolution) return 1;
	if(var == ELagrangeMult) return 1;
    else if (var == EDerivative) {
        return fDim;
    }
	
    return TPZMaterial::NSolutionVariables(var);
}

template<class TVar>
void TPZMatPoissonCS<TVar>::Solution(const TPZVec<TPZMaterialDataT<TVar>> &datavec,
                                     int var, TPZVec<TVar> &solOut)
{
    const auto &solGlo = datavec[0].sol[this->fPostProcIndex];
    const auto &dsolGlo = datavec[0].dsol[this->fPostProcIndex];
    const auto &solLoc = datavec[1].sol[this->fPostProcIndex];
    const auto &dsolLoc = datavec[1].dsol[this->fPostProcIndex];
    const auto &solLag = datavec[2].sol[this->fPostProcIndex];
    const auto &dsolLag = datavec[2].dsol[this->fPostProcIndex];
    fScale = 1.;
    
	if (var == ESolution){
        if (solGlo.size() > 0){
            // if (datavec[0].x[0] >= 1. && datavec[0].x[0] <= 2. && solGlo.size() > 0){
            //     fScale = 0.5;
            // } else {
            //     fScale = 1.;
            // }
            solOut.Resize(solGlo.size());
            // std::cout << "solGlo " << solGlo[0] << std::endl;
            for (int i=0; i<solGlo.size(); i++) {
                solOut[i] = solGlo[i]*fScale;
            }
            return; 
        } else if (solLoc.size() > 0){
            // if (datavec[1].x[0] >= 1. && datavec[1].x[0] <= 2. && solLoc.size() > 0){
            //     fScale = 0.5;
            // } else {
            //     fScale = 1.;
            // }
            solOut.Resize(solLoc.size());
            // std::cout << "solLoc " << solLoc[0] << std::endl;
            for (int i=0; i<solLoc.size(); i++) {
                solOut[i] = solLoc[i]*fScale;
            }
            return; 
        } else {
            DebugStop();
        }
	}
    if (var == EDerivative) {
        if (solGlo.size() > 0){
            solOut.Resize(3);
            for (int i=0; i<fDim; i++) {
                // std::cout << "solout 1 " << solOut[i] << std::endl;
                solOut[i] = dsolGlo.GetVal(i,0)/fScale;
                
            }
            // std::cout << "solout 2 " << solOut[0] << std::endl;
            return;
        } else if (solLoc.size() > 0){
            solOut.Resize(fDim);
            for (int i=0; i<fDim; i++) {
                solOut[i] = dsolLoc.GetVal(i,0)/fScale;
            }
            // std::cout << "solloc 2 " << solOut[0] << std::endl;
            return;
        } else {
            DebugStop();
        }
    }
}

template<class TVar>
void TPZMatPoissonCS<TVar>::Errors(const TPZVec<TPZMaterialDataT<TVar>> &datavec,
                                 TPZVec<REAL> &values) {
    const auto &x = datavec[0].x;
    const auto &u = datavec[0].sol[this->fPostProcIndex];
    const auto &dudx = datavec[0].dsol[this->fPostProcIndex];
    const auto &axes = datavec[0].axes;

#ifdef PZDEBUG
    if(!this->HasExactSol()){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"\nThe material has no associated exact solution. Aborting...\n";
        DebugStop();
    }
#endif

    TPZManVector<TVar,1> u_exact={0.};
    TPZFNMatrix<3,TVar> du_exact(3,1,0.);
    this->ExactSol()(x,u_exact,du_exact);
    values.Resize(this->NEvalErrors());
    values.Fill(0.0);
    TPZManVector<TVar> sol(1),dsol(3,0.);
    TPZFNMatrix<3,TVar> gradu(3,1);
    TPZAxesTools<TVar>::Axes2XYZ(dudx,gradu,axes);
    
    //values[0] : error in H1 norm
    //values[1] : eror in L2 norm
    //values[2] : erro in H1 semi-norm
    TVar diff = (u[0] - u_exact[0]);
    if constexpr (is_complex<TVar>::value){
        values[1]  = std::real((diff*std::conj(diff)));
    }else{
        values[1]  = diff*diff;
    }
  
    values[2] = 0.;

    for(auto id=0; id<fDim; id++) {
      diff = (gradu(id) - du_exact(id,0));
      if constexpr(is_complex<TVar>::value){
          values[2]  += std::real(diff*std::conj(diff));
      }else{
          values[2]  += diff*diff;
      }
    }
    values[0]  = values[1]+values[2];
}

template<class TVar>
int TPZMatPoissonCS<TVar>::ClassId() const{
    return Hash("TPZMatPoissonCS") ^ TBase::ClassId() << 1;
}


template class TPZMatPoissonCS<STATE>;
template class TPZMatPoissonCS<CSTATE>;