/**
 * @file TPZLagrangeMultiplier.h
 * @brief Contains the TPZLagrangeMultiplierCSArlequin class 
 * which implements a lagrange multiplier to be used in FEM formulations
 * with combined spaces.
 */

#ifndef TPZLagrangeMultiplierCSArlequin_H
#define TPZLagrangeMultiplierCSArlequin_H
#include "TPZMatBase.h"
#include "TPZMatCombinedSpaces.h"
#include "TPZMatInterfaceCombinedSpaces.h"

#include "TPZLagrangeMultiplier.h"

/// Material which implements a Lagrange Multiplier for combined spaces
template<class TVar = STATE>
class TPZLagrangeMultiplierCSArlequin :
    public TPZMatBase<TVar,
                      TPZMatCombinedSpacesT<TVar>,
                      TPZMatInterfaceCombinedSpaces<TVar>>,
    public TPZLagrangeMultiplierBase
{
    using TBase = TPZMatBase<TVar,
                             TPZMatCombinedSpacesT<TVar>,
                             TPZMatInterfaceCombinedSpaces<TVar>>;
    bool IsLagrangeMult() override{return true;};
    public :
	/** @brief Simple constructor */
	TPZLagrangeMultiplierCSArlequin() : 
        TBase(),TPZRegisterClassId(&TPZLagrangeMultiplierCSArlequin::ClassId)
    {
        
    }
	/** @brief Constructor with the index of the material object within the vector */
	TPZLagrangeMultiplierCSArlequin(int nummat, int dimension, int nstate=1) :
        TPZRegisterClassId(&TPZLagrangeMultiplierCSArlequin::ClassId),
        TBase(nummat), fNStateVariables(nstate), fDimension(dimension)
    {
        
    }
    
    TPZMaterial *NewMaterial() const override
    {
        return new TPZLagrangeMultiplierCSArlequin(*this);
    }
    
    /** @brief Returns the integrable dimension of the material */
    int Dimension() const override
    {return fDimension;}

    //! @name Lagrange
    /** @{*/
    //! Sets the multiplier.
    virtual void SetMultiplier(const TVar mult)
    {fMultiplier = mult;}
    //! Gets the multiplier.
    TVar Multiplier()
    {return fMultiplier;}
	/** @}*/
	std::string Name() const override
    {return "TPZLagrangeMultiplierCSArlequin";}
    
    void ContributeInterface(const TPZMaterialDataT<TVar> &data,
                             const std::map<int, TPZMaterialDataT<TVar>> &dataleft,
                             const std::map<int, TPZMaterialDataT<TVar>> &dataright,
                             REAL weight,
                             TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef) override;

    
    void FillDataRequirementsInterface(TPZMaterialDataT<TVar> &data,
                                       std::map<int, TPZMaterialDataT<TVar>> &datavec_left,
                                       std::map<int, TPZMaterialDataT<TVar>> &datavec_right) override;


        void Contribute(const TPZVec<TPZMaterialDataT<TVar>> &datavec,
                    REAL weight,
                    TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef) override
    {}
    
    void ContributeBC(const TPZVec<TPZMaterialDataT<TVar>> &datavec,
                      REAL weight, TPZFMatrix<TVar> &ek,
                      TPZFMatrix<TVar> &ef,
                      TPZBndCondT<TVar> &bc) override
    {}
    void Solution(const TPZVec<TPZMaterialDataT<TVar>> &datavec,
                  int var, TPZVec<TVar> &sol) override;
    // {
    // }
    void ContributeBCInterface(const TPZMaterialDataT<TVar> &data,
                               const std::map<int, TPZMaterialDataT<TVar>> &dataleft,
                               REAL weight,
                               TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef,
                               TPZBndCondT<TVar> &bc) override
    {}

    void SolutionInterface(const TPZMaterialDataT<TVar> &data,
                           const std::map<int, TPZMaterialDataT<TVar>> &dataleftvec,
                           const std::map<int, TPZMaterialDataT<TVar>> &datarightvec,
                           int var, TPZVec<TVar> &Solout) override
    {}


    void SolutionInterface(const TPZMaterialDataT<TVar> &data,
                           const std::map<int, TPZMaterialDataT<TVar>> &dataleftvec,
                           const std::map<int, TPZMaterialDataT<TVar>> &datarightvec,
                           int var, TPZVec<TVar> &Solout,
                           TPZCompEl *left,TPZCompEl *right) override
    {}
    
    int NStateVariables() const override
    {return fNStateVariables;}

    /** @brief It returns the variable index associated with the name */
	int VariableIndex(const std::string &name)const override;

    int NSolutionVariables(int var) const override;


    void Print(std::ostream &out) const override;


    /** @{*/
    int ClassId() const override;
	
	void Write(TPZStream &buf, int withclassid) const override;

	void Read(TPZStream &buf, void *context) override;
	
    /** @}*/
protected:
    //! Number of state variables
    int fNStateVariables{0};
    
    //! Dimension associated with the material
    int fDimension;
    
    TVar fMultiplier{1.};
};

#endif /* defined(__PZ__TPZLagrangeMultiplier__) */
