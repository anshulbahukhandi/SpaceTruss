/* EXAMPLE 10.3.1
   Copyright(c) 2005-10, S. D. Rajan
   Object-Oriented Numerical Analysis

  NOTES
   (1) Illustrates the development of a matrix toolbox
       based on the CVector and CMatrix classes
       developed earlier.
   (2) return value of most functions is either true or false.
       a true value signifies no error
       a false value indicates an input data error

*/
/* MODIFIED AS A PART OF ASSIGNMENT BY:
Name : Anshul Bahukhandi
ID:1207178527
PROJECT 1 MATRIX TOOLBOX
*/

#ifndef __ANSHUL_MATTOOLBOX_H__
#define __ANSHUL_MATTOOLBOX_H__

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include <string>
#include "F:\ASU Courses\cee532\Library\Library\vectortemplate.h"
#include "F:\ASU Courses\cee532\Library\Library\matrixtemplate.h"

const int ELEMENTSPERLINE = 5;// # of vector/matrix elements per line
const int NUM_ELEMENTS_PER_LINE=5;
const int FW = 10;              // field width

template <class T>
class CMatToolBox
{
    public:
        CMatToolBox ();
        ~CMatToolBox ();

        // vector-related functions
        void Display (const std::string& szMessage,const CVector<T>& A) const;
        bool Add (const CVector<T>& A, const CVector<T>& B, CVector<T>& C);
        bool Subtract (const CVector<T>& A, const CVector<T>& B, CVector<T>& C);
        bool DotProduct (const CVector<T>& A,const CVector<T>& B, T& product);
        bool Normalize (CVector<T>& A);
        void Scale (CVector<T>& A, T factor);
        T    MaxValue (const CVector<T>& A);
        T    MinValue (const CVector<T>& A);
        T    TwoNorm (const CVector<T>& A);
        T    MaxNorm (const CVector<T>& A);
        bool CrossProduct (const CVector<T>& A,const CVector<T>& B, CVector<T>& C);

        // matrix-related functions
		void PrintMatrixRowWise (CMatrix<T>& A, const std::string& heading,
                         std::ostream& Out);
		void PrintMatrixColumnWise (CMatrix<T>& A, const std::string& heading,
                            std::ostream& Out);

		void Display (const std::string& szMessage,const CMatrix<T>& A) const;
        bool Add (const CMatrix<T>& A, const CMatrix<T>& B, CMatrix<T>& C);
        bool Subtract (const CMatrix<T>& A,const CMatrix<T>& B, CMatrix<T>& C);
        bool Multiply (const CMatrix<T>& A,const CMatrix<T>& B, CMatrix<T>& C);
        bool Determinant (const CMatrix<T>& A, T& det);
        void Scale (CMatrix<T>& A, T scale);
        T    MaxNorm (const CMatrix<T>& A);
        bool Transpose (const CMatrix<T>& A,CMatrix<T>& B);
        bool Inverse (const CMatrix<T>& A,CMatrix<T>& B, T TOL);
        bool MatMultVec (const CMatrix<T>& A,const CVector<T>& x,CVector<T>& b);   
        bool LUFactorization (CMatrix<T>& A, T TOL);							   
        bool LUSolve (const CMatrix<T>& A,CVector<T>& x, const CVector<T>& b);     
        bool AxEqb (CMatrix<T>& A,CVector<T>& x,CVector<T>& b,T TOL);			   
        bool LDLTFactorization (CMatrix<T>& A, T TOL);                             
        bool LDLTSolve (const CMatrix<T>& A, CMatrix<T>& x,const CMatrix<T>& b);   

        // helper function
        void GetFLOPStats (double& dAS, double& dM, double& dD) const;

    private:
        // these are declared double to avoid integer overflow
        double m_dASOP; // # of floating point additions and subtractions
        double m_dMOP;  // # of floating point multiplications
        double m_dDOP;  // # of floating point divisions

    protected:

};

// ctor
template <class T>
CMatToolBox<T>::CMatToolBox ()
// ==================================================================
// Function: default constructor
//    Input: none
//   Output: none
// ==================================================================
{
    m_dASOP = m_dMOP = m_dDOP = 0.0;
}

// dtor
template <class T>
CMatToolBox<T>::~CMatToolBox ()
// ==================================================================
// Function: destructor
//    Input: none
//   Output: none
// ==================================================================
{
}

// ---------------------------------------------------------------
// ------------------------ vector functions ---------------------
// ---------------------------------------------------------------
template <class T>
void CMatToolBox<T>::Display (const std::string& szMessage,
                              const CVector<T>& A) const
// ==================================================================
// Function: displays a message and the elements of a vector
//    Input: message string and the vector 
//   Output: None
// ==================================================================
{
    std::cout << '\n' << szMessage << '\n';
    std::cout.setf(std::ios::left);
    for (int i=1; i <= A.GetSize(); i++)
    {
        std::cout << "(" << i << ") "
                  << std::setw(FW) << A(i) << " ";
        if ((i % ELEMENTSPERLINE) == 0)
			std::cout << "\n";
    }
	std::cout<<"\n..................................................................\n";
}

template <class T>
bool CMatToolBox<T>::Add (const CVector<T>& A, const CVector<T>& B, 
                          CVector<T>& C)
// ==================================================================
// Function: adds two vectors and stores the result in the
//           third vector C = A + B
//    Input: vectors A and B 
//   Output: vector C
// ==================================================================
{
    // check for incompatible vectors
    int n = A.GetSize();
    if (n != B.GetSize() || n != C.GetSize())
        return false;

    // add
    for (int i=1; i <= n; i++)
        C(i) = A(i) + B(i);
    m_dASOP += static_cast<double>(n);       //Remember GOOD STEP!!

    return true;
}

template <class T>
bool CMatToolBox<T>::Subtract (const CVector<T>& A, 
                               const CVector<T>& B, CVector<T>& C)
// ==================================================================
// Function: subtracts one vector from another and stores the result
//           in the third vector C = A - B
//    Input: vectors A and B 
//   Output: vector C
// ==================================================================
{
	int n=A.GetSize();
	if(n!=B.GetSize()||n!=C.GetSize())
	return false;
	for(int i=1;i<=n;i++)
		C(i)=A(i)-B(i);
	m_dASOP+=static_cast<double>(n);
	return true;
}

template <class T>
bool CMatToolBox<T>::DotProduct (const CVector<T>& A,const CVector<T>& B, T& product)
// ==================================================================
// Function: computes the dot product of two vectors such that
//           product = A dot B
//    Input: vectors A and B 
//   Output: product 
// ==================================================================
{
		int n=A.GetSize();
		 product=0;
		if(n!=B.GetSize())
			return false;
		for(int i=1;i<=n;i++)
			product+=A(i)*B(i);
				m_dMOP+=static_cast<double>(n);
				m_dASOP+=static_cast<double>(n-1);
					return true;
}

template <class T>
bool CMatToolBox<T>::Normalize (CVector<T>& A)
// ==================================================================
// Function: normalizes a vector
//    Input: vector A 
//   Output: normalized vector A 
// ==================================================================
{
    T sum=0;
    int n = A.GetSize();
    for (int i=1; i <= n; i++)
        sum += A(i)*A(i);

    m_dASOP += static_cast<double>(n);
    m_dMOP += static_cast<double>(n-1);

    if (sum ==0)
        return false;
    else
    {
        sum = sqrt(sum);
        for (int i=1; i <= n; i++)
            A(i) /= sum;
    }

    m_dDOP += static_cast<double>(n);

    return true;
}

template <class T>
void CMatToolBox<T>::Scale (CVector<T>& A, T c)
// ==================================================================
// Function: scales a vector by a constant c such that A = c A
//    Input: vector A and constant c 
//   Output: scaled vector A
// ==================================================================
{
	int n=A.GetSize();
	for (int i=1;i<=n;i++)
		A(i)=c*A(i);
	m_dMOP+=static_cast<double>(n);
}

template <class T>
T CMatToolBox<T>::MaxValue (const CVector<T>& A)
// ==================================================================
// Function: finds the largest value among all the elements in A
//    Input: vector A 
//   Output: return value is the largest element in A
// ==================================================================
{   int n=A.GetSize();
T max=A(1);
for(int i=1;i<=n;i++)
{if (A(i)>=max)
   max=A(i);
}
return max;    
}

template <class T>
T CMatToolBox<T>::MinValue (const CVector<T>& A)
// ==================================================================
// Function: finds the smallest value among all the elements in A
//    Input: vector A 
//   Output: return value is the smallest element in A
// ==================================================================
{   int n=A.GetSize();
T min=A(1);
for(int i=1;i<=n;i++)
{if (A(i)<=min)
   min=A(i);
}
    return min;
}

template <class T>
T CMatToolBox<T>::TwoNorm (const CVector<T>& A)
// ==================================================================
// Function: computes the two norm of vector A
//    Input: vector A 
//   Output: return value is the two-norm
// ==================================================================
{ double tn,sum=0;
int n=A.GetSize();
for (int i=1;i<=n;i++)
{sum+=A(i)*A(i);}
	tn=pow(sum,0.5);
    m_dMOP+=static_cast<double>(n);
	m_dASOP+=static_cast<double>(n-1);
	return tn;
}

template <class T>
T CMatToolBox<T>::MaxNorm (const CVector<T>& A)
// ==================================================================
// Function: computes the two norm of vector A
//    Input: vector A 
//   Output: return value is the two-norm
// ==================================================================
{T max=abs(A(1));
	int n=A.GetSize();
	for(int i=1;i<=n;i++)
	{if(abs(A(i))>max)
	  max=abs(A(i));
	}
return max;
}
	


template <class T>
bool CMatToolBox<T>::CrossProduct (const CVector<T>& A,const CVector<T>& B,CVector<T>& C)
// ==================================================================
// Function: computes the cross-product of two vectors and stores the
//           result in the third vector such that C = A x B
//           (3-dimensional space)
//    Input: vectors A, B and C
//   Output: vector C
// ==================================================================
{int n1=A.GetSize();
 int n2=B.GetSize();
 if(n1!=3||n2!=3)
     return false;
C(3)=A(1)*B(2)-A(2)*B(1);
C(2)=A(3)*B(1)-B(3)*A(1);
C(1)=A(2)*B(3)-A(3)*B(2);
m_dMOP+=static_cast<double>(6);
m_dASOP+=static_cast<double>(3);
return true;
}

// ---------------------------------------------------------------
// ------------------------ matrix functions ---------------------
// ---------------------------------------------------------------
template <class T>
void CMatToolBox<T>::Display (const std::string& szMessage,
                              const CMatrix<T>& A) const
// ==================================================================
// Function: displays a message and the elements of a matrix
//           rowwise
//    Input: message string and the matrix
//   Output: None
// ==================================================================
{
    std::cout << '\n' << szMessage << '\n';
    std::cout.setf(std::ios::left);
    for (int i=1; i <= A.GetRows(); i++)
    {
        int nC = 0;
        for (int j=1; j <= A.GetColumns(); j++)
        {
            ++nC;
            std::cout << "(" << i << "," << j << ") "
                      << std::setw(FW) << A(i,j) << " ";
            if ((nC % ELEMENTSPERLINE) == 0)
                std::cout << '\n';
        }
        
    }
	std::cout << "\n..................................................................\n";
}

template <class T>
bool CMatToolBox<T>::Add (const CMatrix<T>& A, const CMatrix<T>& B,CMatrix<T>& C)
// ==================================================================
// Function: adds two matrices and stores the result in the
// third matrix C = A + B
// Input: matrices A and B 
// Output: matrix C
// ==================================================================
{
	if(A.GetColumns()!=B.GetColumns()||A.GetRows()!=B.GetRows()||C.GetColumns()!=B.GetColumns()||C.GetRows()!=B.GetRows()||A.GetColumns()!=C.GetColumns()||A.GetRows()!=C.GetRows())
    return false;
for(int i=1;i<=A.GetRows();i++)
{for(int j=1;j<=A.GetColumns();j++)
		C(i,j)=A(i,j)+B(i,j);
}
int n= A.GetRows()*A.GetColumns();
m_dASOP+=static_cast<double>(n);
	return true;
}

template <class T>
bool CMatToolBox<T>::Subtract (const CMatrix<T>& A,const CMatrix<T>& B, CMatrix<T>& C)
// ==================================================================
// Function: subtracts one matrix from another and stores the result
//           in the third matrix C = A - B
//    Input: matrices A and B 
//   Output: matrix C
// ==================================================================

{
	int m=A.GetRows();
	int n=A.GetColumns();
	int o=B.GetRows();
    int p=B.GetColumns();
	
	if(m!=o||n!=p)
    return false;
   
	for(int i=1;i<=A.GetRows();i++)
    {
		for(int j=1;j<=A.GetColumns();j++)
		C(i,j)=A(i,j)-B(i,j);
}
	int q=m*o;
	m_dASOP+=static_cast<double>(q);
	return true;
}

template <class T>
bool CMatToolBox<T>::Multiply (const CMatrix<T>& A,const CMatrix<T>& B, CMatrix<T>& C)
// ==================================================================
// Function: multiplies two matrices and stores the result
//           in the third matrix C = A * B
//    Input: matrices A and B 
//   Output: matrix C
// ==================================================================
{
	if(A.GetColumns()!=B.GetRows())
	return false;
	
		for ( int i=1;i<=A.GetRows();i++)
	{
		
		for( int j=1;j<=B.GetColumns();j++)
		{C(i,j)=0;
		 for ( int k=1;k<=B.GetRows();k++)
		 {
			 C(i,j)+=A(i,k)*B(k,j);
			 m_dMOP++;
		 }
	 }
		m_dASOP+=static_cast<double>(A.GetRows()*B.GetColumns()*(B.GetRows()-1));
	}
	return true;
}


template <class T>
bool CMatToolBox<T>::Determinant (const CMatrix<T>& A, T& det)
// ==================================================================
// Function: computes the determinant of matrix A
//    Input: matrix A and variable to hold the determinant
//   Output: determinant
// ==================================================================
{
    return false;
}

template <class T>
void CMatToolBox<T>::Scale (CMatrix<T>& A, T c)
// ==================================================================
// Function: scales all the elements of a matrix by a constant c
//           such that A = c A
//    Input: matrix A and constant c
//   Output: scaled matrix A
// ==================================================================
{ int n=A.GetColumns();
  int m=A.GetRows();
for(int i=1;i<=m;i++)
{
	for (int j=1;j<=n;j++)
{		A(i,j)=c*A(i,j);
	m_dMOP++;
}
}

}

template <class T>
T CMatToolBox<T>::MaxNorm (const CMatrix<T>& A)
// ==================================================================
// Function: computes the max norm of matrix A
//    Input: matrix A 
//   Output: return value is the max norm
// ==================================================================
{
	int n=A.GetColumns();
    int m=A.GetRows();
	T max=0;
	for (int i =1;i<=m;i++)
	{
		for(int j=1;j<=n;j++)
          {
			  if(abs(A(i,j))>max)
			   max=abs(A(i,j));
	       }
	}
	return max;
}


template <class T>
bool CMatToolBox<T>::Transpose (const CMatrix<T>& A,CMatrix<T>& B)
// ==================================================================
// Function: computes the transpose of a matrix and stores the result
//           in another matrix B = A(T)
//    Input: matrices A and B
//   Output: matrix B
// ==================================================================
{
	if(A.GetRows()!=B.GetColumns()||A.GetColumns()!=B.GetRows())
	{
		return false;
	}

	for(int i=1 ; i<=B.GetRows() ; i++)
	{ 
	for(int j=1;j<=B.GetColumns();j++)
	{
		B(i,j)=A(j,i);
	}
	
	}

	return true;

}


template <class T>
bool CMatToolBox<T>::MatMultVec (const CMatrix<T>& A,const CVector<T>& x,CVector<T>& b)
// ==================================================================
// Function: multiplies a matrix and a vector and stores the result
//           in a vector b = A * x
//    Input: vectors A and x 
//   Output: vector b
// ==================================================================
{
    int n=A.GetColumns();
	int m=A.GetRows();
	int p=x.GetSize();
	int q=b.GetSize();
	float sum=0;
	if(n!=p||m!=q||m!=p||n!=q)
		return false;
	for(int i=1;i<=m;i++)
		{
			for(int j=1;j<=n;j++)
			{
				sum+=A(i,j)*x(j);
			     m_dMOP++;
			} 
			m_dASOP+=(n-1);
			
			b(i)=sum;
	}
				return true;
}

template <class T>
bool CMatToolBox<T>::LUFactorization (CMatrix<T>& A, T TOL)
// ==================================================================
// Function: carries out LU factorization of matrix A
//           A is replaced with L and U
//    Input: matrix A and tolerance value to detect singular A
//   Output: matrix A 
// ==================================================================
{
    return false;
}

template <class T>
bool CMatToolBox<T>::LUSolve (const CMatrix<T>& A,
                              CVector<T>& x,
                              const CVector<T>& b)
// ==================================================================
// Function: carries out forward and backward substitution so as to
//           solve A x = b. A contains L and U terms.
//    Input: matrix A, vectors x and b
//   Output: vector x 
// ==================================================================
{
    return false;
}

template <class T>
bool CMatToolBox<T>::AxEqb (CMatrix<T>& A,
                            CVector<T>& x,
                            CVector<T>& b,
                            T TOL)
// ==================================================================
// Function: solves A x = b using Gaussian Elimination Technique
//           (this version only for one rhs vector)
//    Input: Matrices A and b
//   Output: Matrix x
//           A false return value indicates a singular A matrix.
// ==================================================================
{
    int i, j, k;   // loop indices
    T c;           // multiplier (Step 4)

    // number of equations to solve
    int n = A.GetRows();

    // x initially contains b
    x = b;

    // forward elimination
    for (k=1; k <= n-1; k++)              // Step 1
    {
        if (fabs(A(k,k)) <= TOL)          // Step 2
            return false;
        for (i=k+1; i <= n; i++)          // Step 3
        {
            c = A(i,k)/A(k,k);            // Step 4
            for (j=k+1; j <= n; j++)      // Step 5
                A(i,j) -= c * A(k,j);     // Step 6
            x(i) -= c * x(k);             // Step 8
        }                                 // Step 9
        int nC = n-k;
        if (nC > 0)
        {
            m_dDOP += static_cast<double>(nC);
            m_dMOP += static_cast<double>(nC*nC);
            m_dASOP += static_cast<double>(nC+nC*nC);
        }
    }                                     // Step 10 

    // back substitution
    if (fabs(A(n,n)) <= TOL)              
        return false;
    x(n) /= A(n,n);                       // Step 11

    for (i=n-1; i >= 1; i--)              // Step 12
    {
        T sum = T(0);
        for (j=i+1; j <= n; j++)
            sum += A(i,j) * x(j);         // Step 13
        if ((n-i) > 0)
        {
            m_dASOP += static_cast<double>(n-i);
            m_dMOP += static_cast<double>(n-i);
        }
        x(i) = (x(i) - sum)/A(i,i);       // Step 14
    }                                     // Step 15
    m_dASOP += static_cast<double>(n);
    m_dDOP += static_cast<double>(n+1);

    return true;
}

template <class T>
bool CMatToolBox<T>::LDLTFactorization (CMatrix<T>& A,
                                        T TOL) 
// ==================================================================
// Function: carries out LDL(T) factorization of matrix A
//           A is replaced with L and D. A is a symmetric matrix.
//    Input: matrix A and tolerance value to detect singular A
//   Output: matrix A 
// ==================================================================
{   	

int n=A.GetColumns();   
int m=A.GetRows();
if (m!=n)
	return false;
double sum2=0;
for (int i=1;i<=n;i++)
       
	{
		if (A(i,i)<=TOL)
			return false;

		for(int j=1;j<=i-1;j++)
		{
		A(i,i)=A(i,i)-(A(i,j)*A(i,j)*A(j,j));
		m_dASOP++;
		m_dMOP+=2;
		}
		
		for(int j=i+1;j<=n;j++)
		{
		 for(int k=1;k<=i-1;k++)
		{
			sum2+=A(j,k)*A(k,k)*A(i,k);
		    m_dASOP++; 
			m_dMOP+=2;
		 }
	   A(j,i)=(A(j,i)-sum2)/A(i,i);
	   m_dDOP++;
	   m_dASOP++;
		}	

	}													//CHANGE START
for(int i=1;i<=n;i++)								 //CHANGED THIS ...To Fill up the upper triangular portion
{
	for(int j=i+1;j<=n;j++)		
	{
		A(i,j)=A(j,i);
	}
}
return true ;
}




template <class T>
bool CMatToolBox<T>::LDLTSolve (const CMatrix<T>& A,					//CHANGED THIS..changed  the vectors to  matrix
                                CMatrix<T>& x,
                                const CMatrix<T>& b)
// ==================================================================
// Function: carries out forward and backward substitution so as to
//           solve A x = b. A contains L and D terms.
//    Input: matrix A, vectors x and b
//   Output: vector x 
// ==================================================================
{   double sum1 =0.0;
	int n = A.GetColumns();
	int m = A.GetRows();
	int q=x.GetRows();
	int w=b.GetRows();
	CMatrix<double> y(n,1);
	
	if(n!=m||n!=q||m!=w||q!=w)
		return false;

	
		y(1,1)=b(1,1);
	for (int i=2;i<=n;i++)
	{    
		for (int j=1;j<=i-1;j++)
		{
			sum1+=A(i,j)*y(j,1);
		m_dASOP++;
		m_dMOP++;
		}
		y(i,1)=(b(i,1)-sum1);
	m_dASOP++;
	}
	
	 
	x(n,1)=y(n,1)/A(n,n);
	m_dDOP++;
	for (int i=n-1;i>=1;i--)
	{  
		sum1=0;
		for(int j=i+1;j<=n;j++)
		{
			sum1-=x(j,1)*A(j,i);
		m_dASOP++;
		m_dMOP++;
		}
	x(i,1)=(y(i,1)/A(i,i))+sum1;
	m_dDOP++;
	}
return true;
}

template <class T>
void CMatToolBox<T>::GetFLOPStats (double& dAS, double& dM, 
                                   double& dD) const
// ==================================================================
// Function: retrieves floating point operations
//    Input: variables to store +-, * and / operations
//   Output: variables with their values
// ==================================================================
{
    dAS = m_dASOP;
    dM  = m_dMOP;
    dD  = m_dDOP;
}

template <class T>
void CMatToolBox<T>::PrintMatrixRowWise (CMatrix<T>& A, const std::string& heading,   //CHANGED THIS...Added these functions from matlib.h
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
void CMatToolBox<T>::PrintMatrixColumnWise (CMatrix<T>& A, const std::string& heading,
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
}																			//Change end
	
#endif			
