/*********************************************
Object-Oriented Numerical Analysis
Copyright(c) 2000-10, S. D. Rajan
All rights reserved

Contains matrix-related template helper functions

*********************************************/
#include <iomanip>
#include <algorithm>
#include <functional>
#include <utility>
#include <cmath>
#include "..\library\vectortemplate.h"
#include "..\library\matrixtemplate.h"

const int NUM_ELEMENTS_PER_LINE = 4;

/* ==================================================================
   ======================= helper (non-class) functions =============
   ================================================================== */

template <class T>
void PrintVector (const CVector<T>& A, const std::string& szMessage,
                  std::ostream& Out)
// ==================================================================
// Function: displays a message and the elements of a vector
//    Input: vector, a heading, output stream object 
//   Output: None
// ==================================================================
{
    Out << '\n' << szMessage << '\n';
    Out.setf(std::ios::left);
    for (int i=1; i <= A.GetSize(); i+=NUM_ELEMENTS_PER_LINE)
    {
        for (int k=i; k <= std::min(i+NUM_ELEMENTS_PER_LINE-1, A.GetSize());
             k++)
        {
            Out << "[" << std::setw (4) << k << "]";
            Out << std::setw (15) << A(k) << " ";
        }
        Out << '\n';
    }
}

template <class T>
void Scale (CMatrix<T>& A, T s)
// ---------------------------------------------------------------------------
// Function: Scales matrix A with value s
// ---------------------------------------------------------------------------
{
    int nRows, nColumns;
    nRows = A.GetRows(); nColumns = A.GetColumns();
    for (int i=1; i <= nRows; i++)
    {
        for (int j=1; j <= nColumns; j++)
        {
            A(i,j) *= s;
        }
    }
}

template <class T>
void Shift (CMatrix<T>& A, T s, CMatrix<T> B)
// ---------------------------------------------------------------------------
// Function: Shifts matrix A. A = A + sB 
// ---------------------------------------------------------------------------
{
    int nRows, nColumns;
    nRows = A.GetRows(); nColumns = A.GetColumns();
    for (int i=1; i <= nRows; i++)
    {
        for (int j=1; j <= nColumns; j++)
        {
            A(i,j) += s*B(i,j);
        }
    }
}

template <class T>
void PrintMatrixRowWise (CMatrix<T>& A, const std::string& heading,
                         std::ostream& Out)
// ---------------------------------------------------------------------------
// Function: outputs a matrix into stream Out
// Input:    matrix, a heading, output stream object
// Output:   none
// ---------------------------------------------------------------------------
{
    int i, j, k, nRows, nColumns;
    
    Out << '\n' << heading << '\n';
    Out << setiosflags (std::ios::left);

    nRows = A.GetRows(); nColumns = A.GetColumns();
    for (i=1; i <= nRows; i++)
    {
        Out << "Row No: " << i << '\n';
        for (j=1; j <= nColumns; j=j+NUM_ELEMENTS_PER_LINE)
        {
            for (k=j; k <= std::min(j+NUM_ELEMENTS_PER_LINE-1, nColumns); k++)
            {
                Out << "[" << std::setw (4) << k << "]";
                Out << std::setw (15) << A(i,k) << " ";
            }
            Out << '\n';
        }
    }
}

template <class T>
void PrintMatrixColumnWise (CMatrix<T>& A, const std::string& heading,
                            std::ostream& Out)
// ---------------------------------------------------------------------------
// Function: outputs a matrix into stream Out
// Input:    matrix, a heading, output stream object
// Output:   none
// ---------------------------------------------------------------------------
{
    int i, j, k, nRows, nColumns;
    
    Out << '\n' << heading << '\n';
    Out << setiosflags (std::ios::left);

    nRows = A.GetRows(); nColumns = A.GetColumns();
    for (i=1; i <= nColumns; i++)
    {
        Out << "Column No: " << i << '\n';
        for (j=1; j <= nRows; j=j+NUM_ELEMENTS_PER_LINE)
        {
            for (k=j; k <= std::min(j+NUM_ELEMENTS_PER_LINE-1, nRows); k++)
            {
                Out << "[" << std::setw (4) << k << "]";
                Out << std::setw (15) << A(k,i) << " ";
            }
            Out << '\n';
        }
    }
}

template <class T>
void PrintMatrixColumn (CMatrix<T>& A, const std::string& heading, int i,
                        std::ostream& Out)
// ---------------------------------------------------------------------------
// Function: outputs a column of the matrix into stream Out
// Input:    matrix, a heading, column index, output stream object
// Output:   none
// ---------------------------------------------------------------------------
{
    int j, k, nRows;
    
    Out << '\n' << heading << '\n';
    Out << setiosflags (std::ios::left);

    nRows = A.GetRows();
    for (j=1; j <= nRows; j=j+NUM_ELEMENTS_PER_LINE)
    {
        for (k=j; k <= std::min(j+NUM_ELEMENTS_PER_LINE-1, nRows); k++)
        {
            Out << "[" << std::setw (4) << k << "]";
            Out << std::setw (15) << A(k,i) << " ";
        }
        Out << '\n';
    }
}

template <class T>
bool Transpose (CMatrix<T>& A, const CMatrix<T>& B)
// ---------------------------------------------------------------------------
// Function: computes transpose of a matrix, i.e. A = B(T)
// Input:    matrix B
// Output:   matrix A
// ---------------------------------------------------------------------------
{
    int nRows = A.GetRows();
    int nCols = A.GetColumns();
    if (nRows != B.GetColumns() || nCols != B.GetRows())
        return false;

    for (int i=1; i <= nRows; i++)
    {
        for (int j=1; j <= nCols; j++)
        {
            A(i,j) = B(j,i);
        }
    }

    return true;
}
                
template <class T>
bool Mult (CMatrix<T>& A, const CMatrix<T>& B,
           const CMatrix<T>& C)
// ---------------------------------------------------------------------------
// Function: computes product of two matrices, A = B*C
// Input:    matrices B and C
// Output:   matrix A
// ---------------------------------------------------------------------------
{
    T sum;
    int nRowsA = A.GetRows();
    int nColsA = A.GetColumns();
    int nColsB = B.GetColumns();

    if (nRowsA != B.GetRows() || nColsA != C.GetColumns()
        || nColsB != C.GetRows())
        return false;

    for (int i=1; i <= nRowsA; i++)
    {
        for (int j=1; j <= nColsA; j++)
        {
            sum = 0.0;
            for (int k=1; k <= nColsB; k++)
            {
                sum += B(i,k) * C(k,j);
            }
            A(i,j) = sum;
        }
    }

    return true;
}
                
template <class T>
bool MultMatVec (CVector<T>& A, const CMatrix<T>& B,
                 const CVector<T>& C)
// ---------------------------------------------------------------------------
// Function: computes product of a matrix and a vector, a = B*c
// Input:    Matrix B and vector c
// Output:   Vector a
// ---------------------------------------------------------------------------
{
    T sum;
    int nRowsB = B.GetRows();
    int nColsB = B.GetColumns();

    if (nRowsB != A.GetSize() || nColsB != C.GetSize())
        return false;

    for (int j=1; j <= nRowsB; j++)
    {
        sum = 0.0;
        for (int k=1; k <= nColsB; k++)
        {
            sum += B(j,k) * C(k);
        }
        A(j) = sum;
    }

    return true;
}
                
template <class T>
bool MultMatVec (CVector<T>& A, const CVector<T>& B, const CMatrix<T>& C)
// ---------------------------------------------------------------------------
// Function: computes product of a vector and a matrix, a = b*C
// Input:    vector b and Matrix C
// Output:   vector a
// ---------------------------------------------------------------------------
{
    T sum;
    int nRowsC = C.GetRows();
    int nColsC = C.GetColumns();

    if (nRowsC != A.GetSize() || nColsC != B.GetSize())
        return false;

    for (int j=1; j <= nRowsC; j++)
    {
        sum = 0.0;
        for (int k=1; k <= nColsC; k++)
        {
            sum += B(k)*C(k,j);
        }
        A(j) = sum;
    }

    return true;
}
                
// ============================================================================
// ============================================================================
// ============================= GAUSS ELIMINATION ============================
// ============================================================================
// ============================================================================
template <class T>
int GaussElimination (CMatrix<T>& A, CMatrix<T>& x,
                      const CMatrix<T>& b, T TOL)
// ---------------------------------------------------------------------------
// Function: Solves A x = b
// Input:    A, x, and b
// Output:   A and x are modified. return value is zero if a solution exists
//           Otherwise the eqn. number is returned.
// ---------------------------------------------------------------------------
{
    // solves A x = b
    int i, j, k, ii;
    double c;

    // number of equations to solve
    int n = A.GetRows();
    if (n != A.GetColumns() || n != x.GetRows() || n != b.GetRows() ||
        x.GetColumns() != b.GetColumns())
        return 1;

    // x initially contains b
    x = b;

    // forward elimination
    for (k=1; k <= n-1; k++)
    {
        for (i=k+1; i <= n; i++)
        {
            // singular matrix?
            if (fabs(A(k,k)) <= TOL)
            {
                return k;
            }
            c = A(i,k)/A(k,k);
            for (j=k+1; j <= n; j++)
            {
                A(i,j) -= c * A(k,j);
            }
            x(i,1) -= c * x(k,1);
        }
    }

    // back substitution
    x(n,1) /= A(n,n);

    for (ii=1; ii <= n-1; ii++)
    {
        i = n - ii;
        double sum = 0.0;
        for (j=i+1; j <= n; j++)
        {
            sum += A(i,j) * x(j,1);
        }
        x(i,1) = (x(i,1) - sum)/A(i,i);
    }

    return 0;
}

// ============================================================================
// ============================================================================
// ============================= INVERSE ITERATION ============================
// ============================================================================
// ============================================================================
template <class T>
int InverseIteration (CMatrix<T>& K, CMatrix<T>& M,
                      T& dLambda, CMatrix<T>& u, 
                      int IterationMax, T DTOL,
                      bool bVerbose, std::ostream& Out)
// ---------------------------------------------------------------------------
// Function: Solves K u = (lambda) M u using the Inverse Iteration Method
//           to obtain the lowest eigenvalue and eigenvector
// Input:    K, M, lambda, u (eigenvector), max. # of iterations and
//           convergence tolerance (also used to detect singular K)
//           bVerbose is true if intermediate values are to be printed
//           on the output stream Out
// Output:   lambda and u are modified. 
//           return value is zero if a solution is found
// ---------------------------------------------------------------------------
{
    const int SIZE = K.GetRows();
    CMatrix<double> KSAVE(SIZE,SIZE);
    CMatrix<double> v(SIZE,1), vk(SIZE,1);
    double dLambdaOld;
    int nIterationNum = 0;

    do
    {
        if (bVerbose)
            Out << "\nIteration: " << nIterationNum+1
                << "\n----------\n";
        // v = M u
        Mult (v, M, u);

        // solve K u = v. since K is destroyed, save K.
        KSAVE = K;
		GaussElimination (KSAVE, u, v, DTOL);
        if (bVerbose)
            PrintMatrixColumnWise (u, "Step 3:xhat(k)", Out);

        // v = M u
        Mult (vk, M, u);
        if (bVerbose)
            PrintMatrixColumnWise (vk, "Step 4: yhat(k)", Out);

        // u.v / u.vk
        double dNum = 0.0;
        double dDen = 0.0;
        for (int i=1; i <= SIZE; i++)
        {
            dNum += u(i,1) * v(i,1);
            dDen += u(i,1) * vk(i,1);
        }
        // compute Rayleigh Quotient
        dLambda = dNum/dDen;
        if (bVerbose)
            Out << "Lambda: " << dLambda << '\n';

        // normalize eigenvector
        for (int i=1; i <= SIZE; i++)
            u(i,1) /= dDen;
        if (bVerbose)
            PrintMatrixColumnWise (u, "Step 5: x(k)", Out);

        // check convergence
        if (nIterationNum > 1)
        {
            if (fabs((dLambda - dLambdaOld)/dLambdaOld) <= DTOL)
                return 0;
        }

        // update iteration counter and save current eigenvalue
        ++nIterationNum;
        dLambdaOld = dLambda;
    } while (nIterationNum < IterationMax);

    return 1;
}

// ============================================================================
// ============================================================================
// ============================= JACOBI =======================================
// ============================================================================
// ============================================================================
template <class T>
void Rotate (CMatrix<T>& a, T alpha, T gamma,
             int j, int k, int nq, CVector<T>& temp)
// ---------------------------------------------------------------------------
// function: sets up the rotation matrix to zero out an off-diagonal element
//           used with generalized jacobi method
// ---------------------------------------------------------------------------
{
    int i;

    for (i=1; i <= nq; i++)
        temp(i)=a(i,j)+gamma*a(i,k);
    for (i=1; i <= nq; i++)
        a(i,k)=a(i,k)+alpha*a(i,j);
    for (i=1; i <= nq; i++)
        a(i,j)=temp(i);
    for (i=1; i <= nq; i++)
      temp(i)=a(j,i)+gamma*a(k,i);
    for (i=1; i <= nq; i++)
        a(k,i)=a(k,i)+alpha*a(j,i);
    for (i=1; i <= nq; i++)
        a(j,i)=temp(i);
}

template <class T>
int GeneralizedJacobi (CMatrix<T>& a, CMatrix<T>& b,
                       CMatrix<T>& x, CVector<T>& eigv, int NSMAX, T RTOL,
                       bool bVerbose, std::ostream& Out) 
// ---------------------------------------------------------------------------
// Function: Solves A X = (lambda) B X using the Generalized Jacobi Method
//           to obtain all the eigenpairs
//           a and b are symmetric, positive definite matrices
// Input:    a, b, x, eigv, convergence tolerance and max. # of iterations
//           bVerbose is true if intermediate values are to be printed
//           on the output stream Out
// Output:   x and eigv are modified. 
//           return value is zero if a solution is found
// Source:   Modified version of the FORTRAN subroutine in
//           Bathe and Wilson, Numerical Methods in Finite Element Analysis,
//           Prentice-Hall
//           Equation numbers are from OO Numerical Analysis
// ---------------------------------------------------------------------------
{
    const int n = a.GetRows();     // # of eigenvalues/eigenvectors
    if (n != a.GetColumns() || n != b.GetRows() || n != b.GetColumns() ||
        n != x.GetRows()    || n != x.GetColumns() || n != eigv.GetSize())
        return 1;

    int i, j, k;                   // loop counters
    double eps, eptola, eptolb;
    double akk, ajj, ab, check, sqch, bb;
    double d1, d2, den, ca, cg, xj, xk, tol, dif, epsa, epsb;
    CVector<double> d(n), temp(n); // temporary vectors

    // d stores ratio of diagonal elements of k and m.
    // later it stores the current eigenvalues
    d.Set (0.0);

    // stop if diagonal entry is zero signaling a non-positive matrix
    for (i=1; i <= n; i++)
    {
        if (a(i,i) == 0.0 || b(i,i) == 0.0)
            return 1;
        d(i)=a(i,i)/b(i,i);
        eigv(i)=d(i);
    }

    // set initial value of eigenvector matrix to identity
    x.Set (0.0);
    for (i=1; i <= n; i++)
        x(i,i)=1.0;
    // trivial problem?
    if (n == 1)
        return 0;
    
    // initialize sweep counter and begin iterations
    int nsweep=0;
    int nr=n-1;
    do
    {
        // increment sweep counter
        ++nsweep;
        if (bVerbose)
            Out << "\nIteration: " << nsweep
                << "\n----------\n";

        // check if current off-diagonal element is large enough
        // to require zeroing
        eps = pow(0.01,2*nsweep);
        if (bVerbose)
            Out << "Current eps: " << eps << '\n';
        for (j=1; j <= nr; j++)
        {
            for (k=j+1; k <= n; k++)
            {
                // compute coupling factors (16.4.12)
                eptola=(a(j,k)*a(j,k))/(a(j,j)*a(k,k));
                eptolb=(b(j,k)*b(j,k))/(b(j,j)*b(k,k));
                if (bVerbose)
                    Out << "Coupling factors: (" << j << "," << k
                        << ") eptola " << eptola
                        << " and eptolb " << eptolb << '\n';

                // if elements are small then no zeroing is required
                if (eptola < eps && eptolb < eps) 
                    continue;

                // if zeroing is required, calculate the elements
                // of the rotation matrix ca and cg
                akk=a(k,k)*b(j,k)-b(k,k)*a(j,k); // (16.4.17a)
                ajj=a(j,j)*b(j,k)-b(j,j)*a(j,k); // (16.4.17b)
                ab=a(j,j)*b(k,k)-a(k,k)*b(j,j);  // (16.4.17c)
                check=(ab*ab+4.0*akk*ajj)/4.0;   // (16.4.18)
                if (check < 0.0)
                    return 1;
                sqch=sqrt(check);                // (16.4.18)
                d1=ab/2.0+sqch;                  // (16.4.18)
                d2=ab/2.0-sqch;                  // (16.4.18)
                den=d1;
                den = (fabs(d2) > fabs(d1)? d2 : den);
                if (den == 0.0)
                {
                    ca=0.0;
                    cg=-a(j,k)/a(k,k);           // (16.4.20)
                }
                else
                {
                    ca=akk/den;                  // (16.4.18)
                    cg=-ajj/den;                 // (16.4.18)
                }
                // perform the generalized rotation to zero the current
                // off-diagonal element
                if (bVerbose)
                    Out << "Zeroing out: (" << j << "," << k << ") "
                        << "alpha and beta are " << ca << " " << cg << '\n';
                Rotate (a,ca,cg,j,k,n,temp);
                Rotate (b,ca,cg,j,k,n,temp);

                // update the eigenvector matrix after each sweep
                for (i=1; i <= n; i++)
                {
                    xj=x(i,j);
                    xk=x(i,k);
                    x(i,j)=xj+cg*xk;
                    x(i,k)=xk+ca*xj;
                }
            }
        }

        // update the eigenvalues after each sweep
        for (i=1; i <= n; i++)
        {
            if (a(i,i) <= 0.0 || b(i,i) <= 0.0)
                return 1;
            eigv(i)=a(i,i)/b(i,i);
        }
        if (bVerbose)
            PrintVector (eigv, "Eigenvalues ", Out);

        // convergence check
        bool bhit = false;
        for (i=1; i <= n; i++)
        {
            tol=RTOL*d(i);
            dif=fabs(eigv(i)-d(i));
            if (dif > tol)
            {
                bhit = true;
                break;
            }
        }

        // check all off-diagonal elements to see if another sweep is required
        if (!bhit)
        {
            eps=RTOL*RTOL;
            for (j=1; j <= nr; j++)
            {
                for (k=j+1; k <= n; k++)
                {
                    epsa=(a(j,k)*a(j,k))/(a(j,j)*a(k,k));
                    epsb=(b(j,k)*b(j,k))/(b(j,j)*b(k,k));
                    if (epsa >= eps || epsb >= eps)
                    {
                        bhit = true;
                        break;
                    }
                }
            }
        }

        // zeroing required?
        if (!bhit)
            break;

        // yes. update d matrix and start a new sweep
        for (i=1; i <= n; i++)
            d(i)=eigv(i);
        if (bVerbose)
        {
            PrintMatrixRowWise (a, "Matrix K", Out); 
            PrintMatrixRowWise (b, "Matrix M", Out); 
            PrintMatrixColumnWise (x, "Matrix X", Out); 
        }
    } while (nsweep < NSMAX);  // sweep limit exceeded?

    // scale the eigenvectors
    for (j=1; j <= n; j++)
    {
        bb=sqrt(b(j,j));
        for (k=1; k <= n; k++)
             x(k,j)=x(k,j)/bb;  // (16.4.8)
    }

    // sort the eigenvalues and eigenvectors in ascending order
    for (i=1; i <= nr; i++)
    {
        if (eigv(i+1) < eigv(i))
        {
            // swap eigenvalues
            double dtemp = eigv(i+1);
            eigv(i+1) = eigv(i);
            eigv(i) = dtemp;
            // swap eigenvectors
            for (j=1; j <= n; j++)
            {
                dtemp = x(j,i+1);
                x(j,i+1) = x(j,i);
                x(j,i) = dtemp;
            }
        }
    }

    return 0;
}

template <class T>
void Rotate (CMatrix<T>& a, int i, int j, int k, int l, T tau, T s)
// ---------------------------------------------------------------------------
// Function: Sets up the rotation matrix to zero out an off-diagonal element
//           Used with Jacobi Method
// ---------------------------------------------------------------------------
{
    T g = a(i,j);
    T h = a(k,l);
    a(i,j) = g-s*(h+g*tau);
	a(k,l) = h+s*(g-h*tau);
}

template <class T>
int Jacobi (CMatrix<T>& a, CVector<T>& d, CMatrix<T>&v, T RTOL, int NSMAX,
            bool bVerbose, std::ostream& Out)
// ---------------------------------------------------------------------------
// Function: Solves A v = (d) v using the Jacobi Method
//           to obtain all the eigenpairs
//           A is a symmetric matrix
// Input:    A, v, d, convergence tolerance and max. # of iterations
//           bVerbose is true if intermediate values are to be printed
//           on the output stream Out
// Output:   v and d are modified. 
//           return value is zero if a solution is found
// Source:   Numerical Recipes
// ---------------------------------------------------------------------------
{
    int n = a.GetRows();
    if (n != a.GetColumns() || n != v.GetRows() || n != v.GetColumns() ||
        n != d.GetSize())
        return 1;

	int j,iq,ip,i;
	T tresh,theta,tau,t,sm,s,h,g,c;
    CVector<T> b(n), z(n);

	for (ip=1; ip <= n; ip++)
    {
		for (iq=1; iq <= n; iq++)
            v(ip,iq) = T(0);
		v(ip,ip) = T(1.0);
	}
	for (ip=1; ip <= n; ip++)
    {
		b(ip)=d(ip)=a(ip,ip);
		z(ip) = T(0.0);
	}

	for (i=1; i <= NSMAX; i++)
    {
		sm=0.0;
		for (ip=1; ip <= n-1; ip++)
        {
			for (iq=ip+1; iq <= n; iq++)
				sm += fabs(a(ip,iq));
		}
		if (sm == 0.0)
			break;
		if (i < 4)
			tresh=T(0.2*sm/(n*n));
		else
			tresh=T(0.0);
		for (ip=1; ip <= n-1; ip++)
        {
			for (iq=ip+1; iq <= n; iq++)
            {
				g=100.0*fabs(a(ip,iq));
				if (i > 4 && fabs(d(ip))+g == fabs(d(ip))
					&& fabs(d(iq))+g == fabs(d(iq)))
					a(ip,iq)=0.0;
				else if (fabs(a(ip,iq)) > tresh)
                {
					h=d(iq)-d(ip);
					if (fabs(h)+g == fabs(h))
						t=(a(ip,iq))/h;
					else
                    {
						theta=0.5*h/(a(ip,iq));
						t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0) t = -t;
					}
					c=1.0/sqrt(1.0+t*t);
					s=t*c;
					tau=s/(1.0+c);
					h=t*a(ip,iq);
					z(ip) -= h;
					z(iq) += h;
					d(ip) -= h;
					d(iq) += h;
					a(ip,iq)=0.0;
					for (j=1; j <= ip-1; j++)
						Rotate(a,j,ip,j,iq,tau,s);
					for (j=ip+1; j <= iq-1; j++)
						Rotate(a,ip,j,j,iq,tau,s);
					for (j=iq+1; j <= n; j++)
						Rotate(a,ip,j,iq,j,tau,s);
					for (j=1; j <= n; j++)
						Rotate(v,j,ip,j,iq,tau,s);
				}
			}
		}
		for (ip=1; ip <= n; ip++)
        {
			b(ip) += z(ip);
			d(ip) = b(ip);
			z(ip) = 0.0;
		}
	}

    // sort the eigenvalues and eigenvectors in ascending order
    for (i=1; i <= n-1; i++)
    {
        if (d(i+1) < d(i))
        {
            // swap eigenvalues
            double dTemp = d(i+1);
            d(i+1) = d(i);
            d(i) = dTemp;
            // swap eigenvectors
            for (j=1; j <= n; j++)
            {
                dTemp = v(j,i+1);
                v(j,i+1) = v(j,i);
                v(j,i) = dTemp;
            }
        }
    }

	return 0;
}

// ============================================================================
// ============================================================================
// ============================= BANDED =======================================
// ============================================================================
// ============================================================================
template <class T>
bool CholeskyFactorizationBanded (CMatrix<T>& a, T tol)
// ---------------------------------------------------------------------------
// function: cholesky decomposition 
// input:    a and tolerance tol
// output:   decomposed a. return value is zero if no error
//           otherwise the eqn. number is returned.
// ---------------------------------------------------------------------------
{
    int mband = a.GetColumns();
    int size = a.GetRows();
    T sum, dummy;
    int m, n, jk;

    for (int i=1; i <= size; i++)
    {
        int l = size-i+1;
        if (mband < l) l=mband;
        for (int j=1; j <= l; j++)
        {
            m=mband-j;
            if ((i-1) < m) m=i-1;
            sum=a(i,j);
            for (int k=1; k <= m; k++)
            {
                n=i-k;
                jk=j+k;
                sum -= a(n,k+1)*a(n,jk);
            }
            if (j <= 1)
            {
                if (sum <= tol)
                    return false;
                dummy=1.0/sqrt(sum);
                a(i,j)=dummy;
            }
            else
                a(i,j)=sum*dummy;
        }
    }

    return true;
}

template <class T>
bool CholeskySolveBanded (const CMatrix<T>& a, CVector<T>& x)
// ---------------------------------------------------------------------------
// function: cholesky decomposition solution
// input:    a and rhs vector x
// output:   solution returned in x. 
// ---------------------------------------------------------------------------
{
    int mband = a.GetColumns();
    int size = a.GetRows();
    if (x.GetSize() != size)
        return false;
    T sum;
    int ii, j, l, lk, m;

    // forward substitution solves system (u)t*(y)=(b) for (y)
    x(1)=x(1)*a(1,1);
    for (int i=2; i <= size; i++)
    {
        l=i+1;
        j=l-mband;
        if (l <= mband) j=1;
        sum=x(i);
        m=i-1;
        for (int k=j; k <= m; k++)
        {
            lk=l-k;
            sum=sum-a(k,lk)*x(k);
        }
        x(i)=sum*a(i,1);
    }

    // backward substitution solves the system (u)*(x)=(y) for (x)
    x(size)=x(size)*a(size,1);
    ii=size+1;
    for (int n=2; n <= size; n++)
    {
        int i=ii-n;
        l=i-1;
        j=l+mband;
        if (j > size) j=size;
        sum=x(i);
        m=i+1;
        for (int k=m; k <= j; k++)
        {
            int kl=k-l;
            sum=sum-a(i,kl)*x(k);
        }
        x(i)=sum*a(i,1);
    }

    return true;
}

// ============================================================================
// ============================================================================
// ============================= SKYLINE ======================================
// ============================================================================
// ============================================================================
inline int SkylineKIJ (int rowno, int colno, const CVector<int>& dloc)
// Finds the index where K(rowno, colno) is stored in the
// the skyline vector
{
	int isky, idg;
	
	// valid values?
	assert (colno > 0 && rowno > 0);

	// yes. upper triangle component.
    if (colno >= rowno)
	{
		// column number gives the diagonal location.
		idg=colno;
		isky=dloc(colno)+colno-rowno;
	}
	else
	{
		// row number gives the diagonal location.
		idg=rowno;
		isky=dloc(rowno)+rowno-colno;
	}

	// valid value?
    assert (isky < dloc(idg+1));

	return (isky);
}

template <class T>
bool LDLTFactorizationSkyline (CVector<T>& a, CVector<int>& dloc, 
                               int n, double tol)
{
/*
C
C     DECOMPOSITION ONLY: LDLT TO SOLVE
C         A     X    =   C
C          NxN   Nx1      Nx1
C
C     SKYLINE STORAGE SCHEME (NOTE: UPPER TRIANGULAR OF "A" IS STORED
C                                   NOT THE LOWER TRIANGULAR)
C
*/
    int i, ic, idiag, j, ki, kl, ku, kh, kk, k, klt, l, nd;
    double b, c;

    //  more than one equation?
    if (a(1) <= tol) 
        return false;

    // no.
    if (n == 1)
        return true;

    // yes. loop thro' all equations.
    for (i=1; i <= n; i++)
    {
        // locate diagonal entry.
	    idiag = dloc(i);

        // kl - ku track off-diagonal terms in column.
	    kl = idiag+1;
	    ku = dloc(i+1)-1;

        // find effective column height.
	    kh=ku-kl;

        // column contains several off-diagonal elements.
        if (kh > 0)
        {
            // handle off-diagonal terms in column.
            k=i-kh;
            // 'k' locates the diagonal element.
	        ic=0;

            // start with the last entry in column.
            klt=ku;

            // number of off-diagonal terms is 'kh'.
            for (j=1; j <= kh; j++)
            {
	            ic++;
	            klt=klt-1;
	            ki=dloc(k);
	            nd=dloc(k+1)-ki-1;
	            if (nd > 0)
                {
                    kk = std::min(ic,nd);
	                c=0.0;
	                for (l=1; l <= kk; l++)
	                    c += a(ki+l)*a(klt+l);
                    a(klt) -= c;
                }
                k++;
            }
        }
        // column contains only one off-diagonal element.
	    if (kh >= 0)
        {
            // reduction.
            k=i;
            b=0.0;
            for (kk=kl; kk <= ku; kk++)
            {
                k--;
	            ki = dloc(k);
	            c = a(kk)/a(ki);
	            b = b + c*a(kk);

                // off-diagonal element.
	            a(kk)=c;
            }

            // diagonal element.
	        a(idiag)=a(idiag)-b;
        }
        // column contains only diagonal element.
        // singular/non-positive definite matrix ?
        if (a(idiag) <= tol) 
               return false;
    }

    return true;
}

template <class T>
bool LDLTSolveSkyline (const CVector<T>& a, const CVector<int>& dloc,
                       CVector<T>& x, int n)
{
/*
C     SOLUTION ONLY: LDLT TO SOLVE
C         A     X    =   C
C          NxN   Nx1      Nx1
C
C     SKYLINE STORAGE SCHEME (NOTE: UPPER TRIANGULAR OF "A" IS STORED
C                                   NOT THE LOWER TRIANGULAR)
*/
    int i, kl, ku, kk, k, l, nn;
    double c;

    if (x.GetSize() != n)
        return false;

    // forward substitution.
    for (i=1; i <= n; i++)
    {
        kl=dloc(i)+1;
        ku=dloc(i+1)-1;
        kk=ku-kl;
        if (kk >= 0)
        {
            k=i;
            c=0.0;
            for (kk=kl; kk <= ku; kk++)
            {
                k--;
                c=c+a(kk)*x(k);
            }
            x(i)=x(i)-c;
        }
    }

    // backward substitution.
    for (i=1; i <= n; i++)
    {
      k=dloc(i);
      x(i)=x(i)/a(k);
    }
    if (n == 1)
        return true;
    else
    {
        nn=n;
        for (l=2; l <= n; l++)
        {
            kl=dloc(nn)+1;
            ku=dloc(nn+1)-1;
            kk=ku-kl;
            if (kk >= 0)
            {
                k=nn;
                for (kk=kl; kk <= ku; kk++)
                {
                    k--;
                    x(k)=x(k)-a(kk)*x(nn);
                }
            }
            nn--;
        }
    }

    return true;
}

/*
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  Description:                                                              //
//     This routine uses Crout's method to decompose a row interchanged       //
//     version of the n x n matrix A into a lower triangular matrix L and a   //
//     unit upper triangular matrix U such that A = LU.                       //
//     The matrices L and U replace the matrix A so that the original matrix  //
//     A is destroyed.                                                        //
//     Note!  In Crout's method the diagonal elements of U are 1 and are      //
//            not stored.                                                     //
//     Note!  The determinant of A is the product of the diagonal elements    //
//            of L.  (det A = det L * det U = det L).                         //
//     The LU decomposition is convenient when one needs to solve the linear  //
//     equation Ax = B for the vector x while the matrix A is fixed and the   //
//     vector B is varied.  The routine for solving the linear system Ax = B  //
//     after performing the LU decomposition for A is                         //
//                      Crout_LU_with_Pivoting_Solve.                         //
//     (see below).                                                           //
//                                                                            //
//     The Crout method with partial pivoting is: Determine the pivot row and //
//     interchange the current row with the pivot row, then assuming that     //
//     row k is the current row, k = 0, ..., n - 1 evaluate in order the      //
//     the following pair of expressions                                      //
//       L[i][k] = (A[i][k] - (L[i][0]*U[0][k] + . + L[i][k-1]*U[k-1][k]))    //
//                                 for i = k, ... , n-1,                      //
//       U[k][j] = A[k][j] - (L[k][0]*U[0][j] + ... + L[k][k-1]*U[k-1][j])    //
//                                                                  / L[k][k] //
//                                      for j = k+1, ... , n-1.               //
//       The matrix U forms the upper triangular matrix, and the matrix L     //
//       forms the lower triangular matrix.                                   //
//                                                                            //
//  Arguments:                                                                //
//     double *A       Pointer to the first element of the matrix A[n][n].    //
//     int    pivot[]  The i-th element is the pivot row interchanged with    //
//                     row i.                                                 //
//     int     n       The number of rows or columns of the matrix A.         //
//                                                                            //
//  Return Values:                                                            //
//     true  Success                                                          //
//     false Failure - The matrix A is singular.                              //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
bool Crout_LU_Decomposition_with_Pivoting (CMatrix<double>& A, 
                                           CVector<int>& pivot, int n)
{
   int row, i, j, k, p;
   double *p_k, *p_row, *p_col;
   double max;

//         For each row and column, k = 0, ..., n-1,

   for (k = 0, p_k = A; k < n; p_k += n, k++)
   {
 
//            find the pivot row
      pivot[k] = k;
      max = fabs( *(p_k + k) );
      for (j = k + 1, p_row = p_k + n; j < n; j++, p_row += n)
      {
         if ( max < fabs(*(p_row + k)) )
         {
            max = fabs(*(p_row + k));
            pivot[k] = j;
            p_col = p_row;
         }
      }

//     and if the pivot row differs from the current row, then
//     interchange the two rows.
      if (pivot[k] != k)
         for (j = 0; j < n; j++)
         {
            max = *(p_k + j);
            *(p_k + j) = *(p_col + j);
            *(p_col + j) = max;
         }

//                and if the matrix is singular, return error
      if (*(p_k + k) == 0.0 ) return false;

//      otherwise find the upper triangular matrix elements for row k. 
      for (j = k+1; j < n; j++)
      {
         *(p_k + j) /= *(p_k + k);
      }

//            update remaining matrix
      for (i = k+1, p_row = p_k + n; i < n; p_row += n, i++)
         for (j = k+1; j < n; j++)
            *(p_row + j) -= *(p_row + k) * *(p_k + j);

   }

   return true;
}

////////////////////////////////////////////////////////////////////////////////
//  Description:                                                              //
//     This routine uses Crout's method to solve the linear equation Ax = B.  //
//     This routine is called after the matrix A has been decomposed into a   //
//     product of a lower triangular matrix L and a unit upper triangular     //
//     matrix U without pivoting.  The argument LU is a pointer to the matrix //
//     the superdiagonal part of which is U and the subdiagonal together with //
//     the diagonal part is L. (The diagonal part of U is 1 and is not        //
//     stored.)   The matrix A = LU.                                          //
//     The solution proceeds by solving the linear equation Ly = B for y and  //
//     subsequently solving the linear equation Ux = y for x.                 //
//                                                                            //
//  Arguments:                                                                //
//     double *LU      Pointer to the first element of the matrix whose       //
//                     elements form the lower and upper triangular matrix    //
//                     factors of A.                                          //
//     double *B       Pointer to the column vector, (n x 1) matrix, B.       //
//     int    pivot[]  The i-th element is the pivot row interchanged with    //
//                     row i.                                                 //
//     double *x       Solution to the equation Ax = B.                       //
//     int     n       The number of rows or columns of the matrix LU.        //
//                                                                            //
//  Return Values:                                                            //
//     true  Success                                                          //
//     false Failure - The matrix A is singular.                              //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
bool Crout_LU_with_Pivoting_Solve(CMatrix<double>& LU,
                                  CVector<double>& B,
                                  CVector<double>& pivot,
                                  CVector<double>& x, int n)
{
   int i, k;
   double *p_k;
   double dum;

//         Solve the linear equation Lx = B for x, where L is a lower
//         triangular matrix.                                      
   
   for (k = 0, p_k = LU; k < n; p_k += n, k++)
   {
      if (pivot[k] != k)
      {
          dum = B[k]; B[k] = B[pivot[k]]; B[pivot[k]] = dum;
      }
      x[k] = B[k];
      for (i = 0; i < k; i++)
          x[k] -= x[i] * *(p_k + i);
      x[k] /= *(p_k + k);
   }

//         Solve the linear equation Ux = y, where y is the solution
//         obtained above of Lx = B and U is an upper triangular matrix.
//         The diagonal part of the upper triangular part of the matrix is
//         assumed to be 1.0.

   for (k = n-1, p_k = LU + n*(n-1); k >= 0; k--, p_k -= n)
   {
      if (pivot[k] != k)
      {
          dum = B[k]; B[k] = B[pivot[k]]; B[pivot[k]] = dum;
      }
      for (i = k + 1; i < n; i++)
          x[k] -= x[i] * *(p_k + i);
      if (*(p_k + k) == 0.0)
          return false;
   }
  
   return true;
}
*/