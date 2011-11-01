#include "stdafx.h"
#include "SLAU.h"

SLAU::SLAU()
{
}

SLAU::SLAU(Matrix m)
{
	A = m;
}
SLAU::SLAU(Matrix m, Column c)
{
	A = m;
	b = c;
}

Column &SLAU::Solve(SLAUMethod method)
{
	switch (method)
	{
	case SLAUMethod::Gauss: SolveGauss();
		break;
	case SLAUMethod::Krammer: SolveKrammer();
		break;
	case SLAUMethod::Marching: SolveMarching();
		break;
	default:{}
	}

	return x;
}

void SLAU::SolveGauss()
{
	int n = A.N();
	double *res = new double[n];
	for (int h=0; h<n; h++) 
	{
		for (int i=h; i<n; i++) 
		{
			if (i==h)
				b[h] = b[h]/A.Cell(h, h);
			else
				b[i] = b[i] - A.Cell(i, h)*b[h];
			for (int j=n-1; j>=h; j--) 
			{
				if (i==h) 
					A.Cell(h, j) = A.Cell(h, j) / A.Cell(h, h); 
				else
					A.Cell(i, j) = A.Cell(i, j) - A.Cell(i, h)*A.Cell(h, j);
			}
			
				
		}
	}

	for (int i=n-1; i>=0; i--) 
	{
		res[i] = b[i];
		for (int j=n-1; j>=i+1; j--)
			res[i] = res[i] - A.Cell(i, j) * res[j];	
	}
	
	x = Column(res, n);
}

void SLAU::SolveKrammer()
{
	int n = A.N();
	double *res = new double[n];
	for (int i=0; i<n; i++)
	{
		Matrix C(A);
		C.SetColumn(b, i);
		res[i] = C.Determinant() / A.Determinant();
	}

	x = Column(res, n);
}

void SLAU::SolveMarching()
{
	if (A.IsTridiagonal())
	{
		int n = A.N();
		Matrix Kf(2, A.N());
		Kf.Cell(0, 0) = A.Cell(0, 1) / A.Cell(0, 0);
		Kf.Cell(1, 0) = b[0] / A.Cell(0, 0);

		for (int i=1; i<n; i++)
		{
			double z = A.Cell(i, i) - A.Cell(i, i-1)*Kf.Cell(0, i-1); 
			Kf.Cell(0, i) = A.Cell(i, i+1) / z;
			Kf.Cell(1, i) = (b[i] + A.Cell(i, i-1)*Kf.Cell(1, i-1)) / z;
		}

		double *res = new double[n];
		res[n-1] = (b[n-1] + A.Cell(n-1, n-2)*Kf.Cell(1, n-1)) / (A.Cell(n-1, n-1) - A.Cell(n-1, n-2)*Kf.Cell(0, n-1));

		for (int i = n-2; i>=0; i--)
		{
			res[i] = Kf.Cell(0, i)*res[i+1] + Kf.Cell(1, i);
		}

		x = Column(res, n);
	}
}	
void SLAU::Init(Matrix m)
{
	A = m;
}

void SLAU::Init(Matrix m, Column c)
{
	A = m;
	b = c;
}