#ifndef SLAU_H
#define SLAU_H
#include "MyMatrix.h"
#include "MyColumn.h"

enum SLAUMethod {Gauss, Krammer, Marching};
	
class SLAU
{
private:
	Matrix A;
	Column b;
	Column x;
public:
	SLAU();
	SLAU(Matrix m);
	SLAU(Matrix m, Column c);

	Column &Solve(SLAUMethod method);
	void SolveGauss();
	void SolveKrammer();
	void SolveMarching();
	
	void Init(Matrix m);
	void Init(Matrix m, Column c);
};

#endif