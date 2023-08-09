#include "TPZArlequinSolver.h"
#include "pzcmesh.h"
#include <TPZSSpStructMatrix.h> //symmetric sparse matrix storage
#include <pzstepsolver.h> //for TPZStepSolver
#include "TPZTimer.h"
#include "pzblockdiag.h"
#include "pzbdstrmatrix.h"
#include "TPZEigenSolver.h"
#include "TPZLapackEigenSolver.h"
#include "pzspblockdiagpivot.h"
#include "pzintel.h"
#include "TPZMultiphysicsCompMesh.h"
#ifdef PZ_USING_MKL
#include "TPZSYSMPPardiso.h"
#endif

template<class TVar>
void TPZArlequinSolver<TVar>::Solve(){
  
  switch (fSolverType)
  {
    case EDefault:
      SolveProblemDefault();
      break;
    case ESparse:
      SolveProblemSparse();
      break;
    case ENoCondense:
      SolveProblemNoCondense();
      break;

    default:
      DebugStop();
      break;
  }
  
}


template<class TVar>
void TPZArlequinSolver<TVar>::SolveProblemNoCondense(){
  auto cmesh = fAnalysis->Mesh();
  TPZMultiphysicsCompMesh *cmeshmulti = dynamic_cast<TPZMultiphysicsCompMesh *>(cmesh);

  TPZSSpStructMatrix<REAL> matskl(cmeshmulti);
  fAnalysis->SetStructuralMatrix(matskl);
  TPZStepSolver<STATE> step;
  step.SetDirect(ELDLt);
  fAnalysis->SetSolver(step);
  fAnalysis->Run();
}


template<class TVar>
void TPZArlequinSolver<TVar>::SolveProblemDefault(){
  
    auto cmesh = fAnalysis->Mesh();
    TPZMultiphysicsCompMesh *cmeshmulti = dynamic_cast<TPZMultiphysicsCompMesh *>(cmesh);

    int64_t nMeshsEq = cmeshmulti->MeshVector()[0]->NEquations() 
                     + cmeshmulti->MeshVector()[1]->NEquations();
    int64_t nLagMuEq = cmeshmulti->MeshVector()[2]->NEquations();
    int64_t nTotalEq = cmeshmulti->NEquations();
    TPZMatRed<STATE, TPZFMatrix<REAL>> *matRed = new TPZMatRed<STATE, TPZFMatrix<REAL>>(nTotalEq,nMeshsEq);
    TPZFMatrix<REAL> K00(nMeshsEq,nMeshsEq,0.);
    TPZFMatrix<STATE> rhsFull(nTotalEq,1,0.);
    TPZAutoPointer<TPZGuiInterface> guiInterface;
    TPZSSpStructMatrix<REAL> matskl(cmeshmulti);
    
    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky);
    step.SetMatrix(&K00);
    matRed->SetSolver(&step);

    matskl.Assemble(*matRed,rhsFull,guiInterface);
    matRed->SetF(rhsFull);
    matRed->SetReduced();
    TPZFMatrix<STATE> rhsHigh(nMeshsEq,1,0.),rhsLow(nLagMuEq,1,0.);
    matRed->F1Red(rhsLow);
     
    TPZFMatrix<STATE> solution(nLagMuEq,1,0.);
    TPZFMatrix<REAL> *K11Red = new TPZFMatrix<REAL>(nLagMuEq,nLagMuEq);
    matRed->K11Reduced(*K11Red,rhsLow);
    K11Red->MultiplyByScalar(-1.,*K11Red);
    rhsLow.MultiplyByScalar(-1.,rhsLow);
    // std::cout <<"RHS = " << rhsLow << std::endl;
    // K11Red->Print(std::cout);
    // std::cout <<"Solution = " << solution << std::endl;
    K11Red->SolveDirect(rhsLow,ELDLt);
    TPZFMatrix<STATE> result(nTotalEq,1,0.);
    // std::cout <<"Solution = " << rhsLow << std::endl;    
    matRed->UGlobal(rhsLow,result);
    fAnalysis->Solution()=result;
    fAnalysis->LoadSolution();
    fAnalysis->Solution().Print("Sol",std::cout);
}

template<class TVar>
void TPZArlequinSolver<TVar>::SolveProblemSparse(){
  
  //HERE STARTS THE ITERATIVE SOLVER SET
  auto cmesh = fAnalysis->Mesh();
  int dimension = cmesh->Dimension();
  //Primeiro cria a matriz auxiliar K00 - que será decomposta
  //    TPZSYsmpMatrix<REAL> K00;
//  
  
}



// template<class TVar>
// void TPZArlequinSolver<TVar>::ComputeConditionNumber(TPZSparseMatRed<STATE> &matRed, TPZAutoPointer<TPZMatrix<REAL>> precond){
  
//   TPZFMatrix<REAL> KBDInv;
//   TPZAutoPointer<TPZFMatrix<REAL>> Res = new TPZFMatrix<REAL>;
//   auto dim = precond->Rows();
//   Res->Redim(dim,dim);
//   precond->Inverse(KBDInv,ELU);
// //   KBDInv.Print("KBD=",std::cout,EMathematicaInput);
//   // KBDInv.Identity();
//   // KBDInv.Print("KBDInv=",std::cout,EMathematicaInput);
//   // matRed.K11().Print("K11=",std::cout,EMathematicaInput);
//   // matRed.MultAdd(KBDInv,*Res,*Res,1.,0.);
//   TPZFMatrix<REAL> F1,Aux(dim,dim,0);
//   matRed.K11Reduced(Aux,F1);
//   Aux.Multiply(KBDInv,Res);
  
// //   Res->Print("Res=",std::cout,EMathematicaInput);
  
//   TPZLapackEigenSolver<REAL> eigSolver;
  
//   TPZVec<std::complex<REAL>> eigenvalues;
//   eigSolver.SetMatrixA(Res);
//   auto a1 = eigSolver.SolveEigenProblem(eigenvalues);
  
//   std::ofstream rprint3,rprint4;
//   rprint3.open("REAL_EIGEN_ALL.txt",std::ios_base::app);
//   rprint4.open("REAL_EIGEN.txt",std::ios_base::app);
  
//   REAL maxEig = 0.;
//   REAL minEig = 1e3;
//   REAL minAbs = 1e3;
//   REAL maxAbs = 0.;
//   int nonzeroEigenvalues = 0;
//   REAL tol = 1e-8;
//   for (int i = 0; i < eigenvalues.size(); i++)
//   {
//     rprint3 << eigenvalues[i].real() << std::endl;
//     if (eigenvalues[i].real() > maxEig) maxEig = eigenvalues[i].real();
//     if (eigenvalues[i].real() < minEig) minEig = eigenvalues[i].real();
//     if (fabs(eigenvalues[i].real()) < minAbs) minAbs = fabs(eigenvalues[i].real());
//     if (fabs(eigenvalues[i].real()) > maxAbs) maxAbs = fabs(eigenvalues[i].real());
//     if (fabs(eigenvalues[i].real()) > tol) nonzeroEigenvalues++;
//   }
//   rprint3 << std::endl;
  
//   rprint4 << maxEig << " " << minEig << " " << maxAbs << " " << minAbs << " " << nonzeroEigenvalues << "/" << dim<< std::endl;
//   std::cout << maxEig << " " << minEig << " " << maxAbs << " " << minAbs << " " << nonzeroEigenvalues << "/" << dim<< std::endl;
  
  
// }

template<class TVar>
void TPZArlequinSolver<TVar>::ComputeConditionNumber(TPZMatRed<STATE,TPZFMatrix<STATE>> &matRed, TPZAutoPointer<TPZMatrix<REAL>> precond){
  
  TPZFMatrix<REAL> KBDInv;
  TPZAutoPointer<TPZFMatrix<REAL>> Res = new TPZFMatrix<REAL>;
  auto dim = precond->Rows();
  // Res->(dim,dim,true);
  precond->Inverse(KBDInv,ELU);
  // KBDInv.Identity();
//   KBDInv.Print("KBDInv=",std::cout,EMathematicaInput);
//   matRed.Print("MatRed=",std::cout,EMathematicaInput);
  
  matRed.Multiply(KBDInv,Res);
  // KBDInv.Multiply(*matRed,Res);
  
  
  
//   Res->Print("Res=",std::cout,EMathematicaInput);
  
  TPZLapackEigenSolver<REAL> eigSolver;
  
  TPZVec<std::complex<REAL>> eigenvalues;
  eigSolver.SetMatrixA(Res);
  
  // auto a1 = eigSolver.SolveEigenProblem(eigenvalues);
  
  std::ofstream rprint3,rprint4;
  rprint3.open("REAL_EIGEN_ALL.txt",std::ios_base::app);
  rprint4.open("REAL_EIGEN.txt",std::ios_base::app);
  
  REAL maxEig = 0.;
  REAL minEig = 1e3;
  REAL minAbs = 1e3;
  REAL maxAbs = 0.;
  int nonzeroEigenvalues = 0;
  REAL tol = 1e-10;
  for (int i = 0; i < eigenvalues.size(); i++)
  {
    rprint3 << eigenvalues[i].real() << std::endl;
    if (eigenvalues[i].real() > maxEig) maxEig = eigenvalues[i].real();
    if (eigenvalues[i].real() < minEig) minEig = eigenvalues[i].real();
    if (fabs(eigenvalues[i].real()) < minAbs) minAbs = fabs(eigenvalues[i].real());
    if (fabs(eigenvalues[i].real()) > maxAbs) maxAbs = fabs(eigenvalues[i].real());
    if (fabs(eigenvalues[i].real()) > tol) nonzeroEigenvalues++;
  }
  rprint3 << std::endl;
  
  rprint4 << maxEig << " " << minEig << " " << maxAbs << " " << minAbs << " " << nonzeroEigenvalues << "/" << dim<< std::endl;
  std::cout << maxEig << " " << minEig << " " << maxAbs << " " << minAbs << " " << nonzeroEigenvalues << "/" << dim<< std::endl;
  
  
}



template<class TVar>
void TPZArlequinSolver<TVar>::ThresholdPermeability(REAL threshold){
  
  auto cmeshanalysis = fAnalysis->Mesh();
  TPZMultiphysicsCompMesh *mCmesh = dynamic_cast<TPZMultiphysicsCompMesh *> (cmeshanalysis);
  auto cmesh = mCmesh->MeshVector()[0];
  const int64_t nEquations = fAnalysis->Mesh()->NEquations();
  
  
  std::set<int64_t> remove_eq;
  
  for (TPZCompEl* cel : cmesh->ElementVec()){
    if (!cel) continue;
    int dim = cel->Dimension();
    if (dim != cmesh->Dimension()) continue;
    
    TPZGeoEl* gel = cel->Reference();
    if (!gel) DebugStop();
    
    int nFaces = gel->NSides(dim-1);
    int nSides = gel->NSides();
    
    for (int iFace = 0; iFace < nFaces; iFace++){
      // Pega o geoElSide de cada face, depois pega CenterX.
      
      TPZGeoElSide gelside(gel,nSides-1-nFaces+iFace);
      TPZManVector<REAL,3> xCenter(3,0.), perm(3,0.);
      gelside.CenterX(xCenter);
      
      perm = fPermFunction(xCenter);
      
      if (perm[0] < threshold){
        
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
        if (!intel) DebugStop();
        int64_t conindex1 = cel->ConnectIndex(2*iFace);
        int64_t conindex2 = cel->ConnectIndex(2*iFace+1);
        
        TPZConnect &c1 = cmesh->ConnectVec()[conindex1];
        TPZConnect &c2 = cmesh->ConnectVec()[conindex2];
        
        // std::cout << "Dim = " << dim << std::endl;
        
        int64_t seqnum = c1.SequenceNumber();
        const auto pos = cmesh->Block().Position(seqnum);
        const auto blocksize = cmesh->Block().Size(seqnum);
        // if (blocksize == 0){
        //     continue;
        // }
        // const auto vs = fActiveEquations.size();
        if (pos < nEquations){
          // fActiveEquations.Resize(vs + blocksize);
          for (auto ieq = 0; ieq < blocksize; ieq++) {
            remove_eq.insert(pos + ieq);
          }
        }
        
        // std::cout << "Connect index = " << conindex1 << " " << blocksize << std::endl;
        
        // int64_t seqnum2 = c2.SequenceNumber();
        // const auto pos2 = cmesh->Block().Position(seqnum2);
        // const auto blocksize2 = cmesh->Block().Size(seqnum2);
        // if (blocksize == 0){
        //     continue;
        // }
        // if (c2.IsCondensed()){
        //     std::cout << "Is Condensed " << seqnum2 << std::endl;
        //     continue;
        // }
        // // std::cout << "Connect index = " << conindex2 << " " << blocksize2 << std::endl;
        // const auto vs2 = fActiveEquations.size();
        // if (pos2 < nEquations){
        //     fActiveEquations.Resize(vs2 + blocksize2);
        //     for (auto ieq = 0; ieq < blocksize2; ieq++) {
        //         fActiveEquations[vs2 + ieq] = pos2 + ieq;
        //     }
        // }
        
        
        // int64_t index;
        // //Sets the new connect order
        // c2.SetOrder(0,index);//Div constant function
        // c1.SetOrder(0,index);//High order div-free functions
        
        // //Gets connect information to update block size (stiffness matrix data structure)
        // int64_t seqnum = c2.SequenceNumber();
        // int nvar = 1;
        // TPZMaterial * mat = cel->Material();
        // if (mat) nvar = mat->NStateVariables();
        // int nshape = intel->NConnectShapeF(2*iFace,c1.Order());
        // int nshape2 = intel->NConnectShapeF(2*iFace+1,c2.Order());
        // c2.SetNShape(nshape2);
        // // c.SetNState(nvar);
        // cmesh->Block().Set(seqnum, nvar * nshape2);
        
        
        // seqnum = c1.SequenceNumber();
        // // nshape = 1;
        // cmesh->Block().Set(seqnum, nvar * nshape);
        
      }
    }
  }
  cmesh->InitializeBlock();
  // cmesh->ExpandSolution();
  
  std::cout << "remove_eq = " ;
  for (auto it:remove_eq)
  {
    std::cout << it << " " ;
  }
  
  fActiveEquations.Resize(nEquations - remove_eq.size());
  std::cout << "NEquations = " << nEquations << std::endl;
  std::cout << "remove_eq = " << remove_eq.size() << std::endl;
  std::cout << "fActiveEquations = " << fActiveEquations.size() << std::endl;
  
  int64_t count = 0;
  for (int i = 0; i < nEquations; i++)
  {
    if (remove_eq.find(i) == remove_eq.end()) {
      fActiveEquations[count] = i;
      std::cout << "Active equation [" << count << "] = " << fActiveEquations[count] << std::endl;
      count++;
    }
  }
  
}



template class TPZArlequinSolver<STATE>;
