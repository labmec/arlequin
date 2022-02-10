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
    const auto &phi = datavec[0].phi;
    const auto &dphi = datavec[0].dphix;
    const auto nshape = datavec[0].phi.Rows();
	for(int i = 0; i < nshape; i++){
		for(int j = 0; j < nshape; j++){
            STATE dphiIdphiJ = 0;
            for(int x = 0; x < fDim; x++){
                dphiIdphiJ += dphi.GetVal(x,i) * dphi.GetVal(x,j);
            }
            ek(i, j) += weight*fScale*dphiIdphiJ;
        }//forj
        for(auto l = 0; l < nLoads; l++)
            ef(i,l) += weight*fScale*phi.GetVal(i,0)*force[l];
    }//for i
}

template<class TVar>
void TPZMatPoissonCS<TVar>::ContributeBC(const TPZVec<TPZMaterialDataT<TVar>> &datavec,
                                         REAL weight,
                                         TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef,
                                         TPZBndCondT<TVar> &bc)
{
	
	const auto &phi = datavec[0].phi;
    const auto &dphi = datavec[0].dphix;
	const int phr = phi.Rows();
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
            for(auto in = 0 ; in < phr; in++) {
                for(auto l=0; l < nloads; l++)
                    ef(nvars*in+iv,l) +=
                        (TVar)TPZMaterial::fBigNumber * v2[nvars*l+iv] * (TVar)phi.GetVal(in,0) * (TVar)weight;
                for (auto jn = 0 ; jn < phr; jn++) {
                    ek(nvars*in+iv,nvars*jn+iv) += TPZMaterial::fBigNumber * phi.GetVal(in,0) * phi.GetVal(jn,0) * weight;
                }//jn
            }//in
        }//iv
        break;
    }
		// Neumann condition
    case 1 : {
        for(auto iv = 0; iv < nvars; iv++){
            for(auto in = 0 ; in < phr; in++) {
                for(auto l = 0; l < nloads; l++)
                    ef(nvars*in+iv,l) += v2[nvars*l+iv] * (TVar)fScale * (TVar)phi.GetVal(in,0) * (TVar)weight;
            }//in
        }//iv
        break;
    }
        //Robin condition
    case 2 : {
        for(auto iv = 0; iv < nvars; iv++){
            for(auto in = 0 ; in < phr; in++) {
                for(auto l=0; l < nloads; l++)
                    ef(nvars*in+iv,l) +=
                        (TVar)TPZMaterial::fBigNumber * v2[nvars*l+iv] * (TVar)phi.GetVal(in,0) * (TVar)weight;
                for (auto jn = 0 ; jn < phr; jn++) {
                    ek(nvars*in+iv,nvars*jn+iv) +=
                        TPZMaterial::fBigNumber * v1.GetVal(iv,0) * dphi.GetVal(0,in) * dphi.GetVal(0,jn) * weight;
                }//jn
            }//in
        }//iv
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
	return TPZMaterial::VariableIndex(name);
}

template<class TVar>
int TPZMatPoissonCS<TVar>::NSolutionVariables(int var) const{
	if(var == ESolution) return 1;
    else if (var == EDerivative) {
        return fDim;
    }
	
    return TPZMaterial::NSolutionVariables(var);
}

template<class TVar>
void TPZMatPoissonCS<TVar>::Solution(const TPZVec<TPZMaterialDataT<TVar>> &datavec,
                                     int var, TPZVec<TVar> &solOut)
{
    const auto &sol = datavec[0].sol[this->fPostProcIndex];
    const auto &dsol = datavec[0].dsol[this->fPostProcIndex];
	if (var == ESolution){
        solOut.Resize(sol.size());
        for (int i=0; i<sol.size(); i++) {
            solOut[i] = sol[i]/fScale;
        }
		return;
	}
    if (var == EDerivative) {
        solOut.Resize(fDim);
        for (int i=0; i<fDim; i++) {
            solOut[i] = dsol.GetVal(i,0)/fScale;
        }
        return;
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