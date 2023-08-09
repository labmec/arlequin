//
// Created by Jeferson Fernandes on 12/05/22.
//

#ifndef TPZ_ARLEQUIN_SOLVER_H
#define TPZ_ARLEQUIN_SOLVER_H

#include "pzmatred.h"
#include "TPZArlequinMatRed.h"
#include "TPZLinearAnalysis.h"

template <class TVar>
class TPZArlequinSolver {
public:
    enum SolverType{EDefault, ESparse, ENoCondense}; 

    TPZArlequinSolver() = default;

    TPZArlequinSolver(TPZLinearAnalysis &an, SolverType sType = EDefault, std::function<TPZManVector<STATE,3>(const TPZVec<REAL> &coord)> permFunction = nullptr){
        fAnalysis = &an;
        fSolverType = sType;
        fPermFunction = permFunction;
    };

    void Solve();

    void SolveProblemDefault();

    void SolveProblemSparse(); 
    void SolveProblemNoCondense(); 
    
    // void ComputeConditionNumber(TPZSparseMatRed<STATE> &matRed, TPZAutoPointer<TPZMatrix<REAL>> precond);
    void ComputeConditionNumber(TPZMatRed<STATE,TPZFMatrix<STATE>> &matRed, TPZAutoPointer<TPZMatrix<REAL>> precond);
    void ThresholdPermeability(REAL threshold);

protected:
    SolverType fSolverType;

    TPZLinearAnalysis *fAnalysis;

    TPZVec<int64_t> fActiveEquations;

    std::function<TPZManVector<STATE,3>(const TPZVec<REAL> &coord)> fPermFunction;
};

#endif