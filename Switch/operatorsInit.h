#ifndef OPERATORSINIT
#define OPERATORSINIT

#include <vector>
class QI_Model;
class QI_Matrix;

void GetSelfHamiltonian(QI_Model& model, QI_Matrix& pH_self);
void GetRelaxationMatrices(QI_Model& m, std::vector<QI_Matrix*>& L_matrices);
void GetInteractionHamiltonian(QI_Model& m, QI_Matrix& H_int);

#endif

