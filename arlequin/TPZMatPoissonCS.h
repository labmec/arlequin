/**
 * @file TPZMatPoissonCS.h
 * @brief Contains the TPZMatPoissonCS class which a H1 formulation of the Poisson equation
 */

#ifndef TPZMatPoissonCS_H
#define TPZMatPoissonCS_H

#include "TPZMatBase.h"
#include "TPZMatCombinedSpaces.h"
#include "TPZMatLoadCases.h"
#include "TPZMatErrorCombinedSpaces.h"


/**
 * @brief Implements a H1 formulation of the poisson equation.
 * This material uses a scalar approximation space in a geometric domain of
 * dimension defined by the user (1, 2 or 3D).
 * It can solve for multiple values of rhs at once, but for calculating the error,
 * it will only consider one solution at a time.
 */
template<class TVar=STATE>
class TPZMatPoissonCS :
    public TPZMatBase<TVar,
                      TPZMatCombinedSpacesT<TVar>,
                      TPZMatErrorCombinedSpaces<TVar>,
                      TPZMatLoadCases<TVar>>{
	using TBase = TPZMatBase<TVar,
                             TPZMatCombinedSpacesT<TVar>,
                             TPZMatErrorCombinedSpaces<TVar>,
                             TPZMatLoadCases<TVar>>;
public:
    //! Default constructor
    TPZMatPoissonCS() = default;
    /**
	 * @brief Class constructor 
	 * @param id material id
	 * @param dim problem dimension
	 * @param nstate number of state variables
	 */
	TPZMatPoissonCS(int id, int dim);

    std::string Name() const override { return "TPZMatPoissonCS"; }
	
	/** @brief Solution indices of post-processing */
	enum ESolutionVars { ENone = 0, ESolution = 1 , EDerivative = 2, ELagrangeMult = 3};
	
    /**
     * @brief Set a scale factor for the stiffness matrix and right hand side
     * the default value of the scale factor is 1
     */
    void SetScaleFactor(REAL scale)
    {fScale = scale;}
    
    [[nodiscard]] REAL ScaleFactor() const
    {return fScale;}

	int Dimension() const  override { return this->fDim; }

    int NStateVariables() const override{ return 1;}
    //H1 norm L2 norm H1 seminorm
    int NEvalErrors() const override { return 3;}
    /** @brief Sets problem dimension */
    virtual void SetDimension(int dim) { this->fDim = dim; }

    void Contribute(const TPZVec<TPZMaterialDataT<TVar>> &datavec,
                    REAL weight,
                    TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef) override;

    void ContributeBC(const TPZVec<TPZMaterialDataT<TVar>> &datavec, REAL weight,
                      TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef,
                      TPZBndCondT<TVar> &bc) override;
    /** @brief To create another material of the same type */
	TPZMaterial * NewMaterial() const override;
	
	/** @brief It returns the variable index associated with the name */
	int VariableIndex(const std::string &name)const override;
	
	int NSolutionVariables(int var) const override;

    void Solution(const TPZVec<TPZMaterialDataT<TVar>> &datavec,
                  int var, TPZVec<TVar> &solOut) override;
    
    void GetSolDimensions(uint64_t &u_len,
                          uint64_t &du_row,
                          uint64_t &du_col) const;

    void Errors(const TPZVec<TPZMaterialDataT<TVar>> &datavec,
                TPZVec<REAL> &errors) override;
    
    virtual int ClassId() const override;
protected:
	
	/** @brief Problem dimension */
	int fDim;
	
	/** @brief Constant solution vector. Ignored if forcing function is set. */
	TPZVec<TVar> fSol;
    
    /** @brief Scale factor applied to the stiffness matrix and right hand side */
    REAL fScale{1};
};


extern template class TPZMatPoissonCS<STATE>;
extern template class TPZMatPoissonCS<CSTATE>;
#endif
