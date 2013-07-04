#pragma once
/*
The EdgePreservingDecomposition files contain standard C++ (standard except the first line) code for creating and, to a
limited extent (create your own uses!), messing with multi scale edge preserving decompositions of a 32 bit single channel
image. As a byproduct it contains a lot of linear algebra which can be useful for optimization problems that
you want to solve in rectangles on rectangular grids.

Anyway. Basically, this is an implementation of what's presented in the following papers:
	Edge-Preserving Decompositions for Multi-Scale Tone and Detail Manipulation
	An Iterative Solution Method for Linear Systems of Which the Coefficient Matrix is a Symetric M-Matrix
	Color correction for tone mapping
	Wikipedia, the free encyclopedia

First one is most of what matters, next two are details, last everything else. I did a few things differently, especially:
	Reformulated the minimization with finite elements instead of finite differences. This results in better conditioning,
	slightly better accuracy (less artifacts), the possibility of a better picked edge stopping function, but more memory consumption.

	A single rotationally invariant edge stopping function is used instead of two non-invariant ones.

	Incomplete Cholseky factorization instead of Szeliski's LAHBF. Slower, but not subject to any patents.

	For tone mapping, original images are decomposed instead of their logarithms, and just one decomposition is made;
	I find that this way works plenty good (theirs isn't better or worse... just different) and is simpler.

Written by ben_pcc in Portland, Oregon, USA. Some history:
	Late April 2010, I develop interest in this stuff because photos of my ceramics lack local contrast.
	Mid 2010, it works but is too slow to be useful.
	Fall 2010, various unsuccessful attempts at speeding up are tried.
	Early December 2010, I get off the path of least resistance and write a matrix storage class with incomplete Cholesky decomposition.
	31 December 2010, the FEM reformulation works very well.
	1 January 2011, I'm cleaning up this file and readying it for initial release.

This release: 1 January 2011.

It's likely that I'll take apart and rerelease contents of this file (in the distant future) as most of it isn't edge preserving decomposition
and rather supporting material. SparseConjugateGradient alone is a workhorse I and a few others have been exploiting for a few years.

EdgePreservingDecomposition.h and EdgePreservingDecomposition.cpp are released under the following licence:
	• It's free.
	• You may not incorporate this code as part of proprietary or commercial software, but via freeware you may use its output for profit.
	• You may modify and redistribute, but keep this big comment block intact and not for profit in any way unless I give specific permission.
	• If you're unsure about anything else, treat as public domain.
	• Don't be a dick.

My email address is my screen name followed by @yahoo.com. I'm also known as ben_s or nonbasketless. Enjoy!

http://www.rawtherapee.com/forum/viewtopic.php?t=2553&start=0&postdays=0&postorder=asc&highlight=


http://web.cecs.pdx.edu/~ben_s/Linked/EPD/EdgePreservingDecomposition.rar

To put this into RT, you should basically:
¥ Add the two files to the RT project.
¥ Write a function that converts CIELab L into 32 bit (0 to 1) continuous memory.
¥ Feed that into EdgePreservingDecomposition::CompressDynamicRange.
¥ Hope that untested ModifyCIELabLuminance will replace the 16 bit source L you started with 
  (my testbed didn't involve 16 bit Lab... actual results in RT would look a bit less cartoony I think).
¥ Link the function arguments of EdgePreservingDecomposition::CompressDynamicRange to the GUI.


Compression, 0
DetailBoost, 0
EdgeStopping, 1.4
Scale, 1.0
SaturationControl, 0.8
Iterates, 20

These are all arguments to the EdgePreservingDecomposition::CompressDynamicRange function,
except SaturationControl which is fed to ModifyCIELabLuminance. 
Iterates is somewhat optional, probably, more is slower but more accurate. 
"More accurate" means "better detail boosting", fewer iterates (even 1) don't seem to result in more artifact.

*/



#include <cmath>
#include <iostream>

//This is for solving big symmetric positive definite linear problems.
float *SparseConjugateGradient(void Ax(float *Product, float *x, void *Pass), float *b, unsigned int n, bool OkToModify_b = true, float *x = NULL, float RMSResidual = 0.0f, void *Pass = NULL, unsigned int MaximumIterates = 0, void Preconditioner(float *Product, float *x, void *Pass) = NULL);

//Storage and use class for symmetric matrices, the nonzero contents of which are confined to diagonals.
class MultiDiagonalSymmetricMatrix{
public:
	MultiDiagonalSymmetricMatrix(unsigned int Dimension, unsigned int NumberOfDiagonalsInLowerTriangle);
	~MultiDiagonalSymmetricMatrix();

	/* Storage of matrix data, and a function to create memory for Diagonals[index].
	Here's the storage scheme, designed to be simple both on paper and in C++:

	Let's say you have some diagonal. The StartRows is the row on which, at the left edge of the matrix, the diagonal "starts",
	and StartRows must strictly increase with its index. The main diagonal for example has start row 0, its subdiagonal has 1, etc.
	Then, Diagonal[j] is the matrix entry on the diagonal at column j. For efficiency, you're expected to learn this and fill in
	public Diagonals manually. Symmetric matrices are represented by this class, and all symmetry is handled internally, you
	only every worry or think about the lower trianglular (including main diagonal) part of the matrix.
	*/
	float **Diagonals;
	unsigned int *StartRows;
	bool CreateDiagonal(unsigned int index, unsigned int StartRow);
	unsigned int n, m;	//The matrix is n x n, with m diagonals on the lower triangle. Don't change these. They should be private but aren't for convenience.
	inline unsigned int DiagonalLength(unsigned int StartRow){	//Gives number of elements in a diagonal.
		return n - StartRow;
	};

	//Not efficient, but you can use it if you're lazy, or for early tests. Returns false if the row + column falls on no loaded diagonal, true otherwise.
	bool LazySetEntry(float value, unsigned int row, unsigned int column);

	//Calculates the matrix-vector product of the matrix represented by this class onto the vector x.
	void VectorProduct(float *Product, float *x);

	//Given the start row, attempts to find the corresponding index, or -1 if the StartRow doesn't exist.
	int FindIndex(unsigned int StartRow);

	//This is the same as above, but designed to take this class as a pass through variable. By this way you can feed
	//the meat of this class into an independent function, such as SparseConjugateGradient.
	static void PassThroughVectorProduct(float *Product, float *x, void *Pass){
		((MultiDiagonalSymmetricMatrix *)Pass)->VectorProduct(Product, x);
	};

	/* CreateIncompleteCholeskyFactorization creates another matrix which is an incomplete (or complete if MaxFillAbove is big enough)
	LDLt factorization of this matrix. Storage is like this: the first diagonal is the diagonal matrix D and the remaining diagonals
	describe all of L except its main diagonal,	which is a bunch of ones. Read up on the LDLt Cholesky factorization for what all this means.
	Note that VectorProduct is nonsense. More useful to you is CholeskyBackSolve which fills x, where LDLt x = b. */
	bool CreateIncompleteCholeskyFactorization(unsigned int MaxFillAbove = 0);
	void KillIncompleteCholeskyFactorization(void);
	void CholeskyBackSolve(float *x, float *b);
	MultiDiagonalSymmetricMatrix *IncompleteCholeskyFactorization;

	static void PassThroughCholeskyBackSolve(float *Product, float *x, void *Pass){
		((MultiDiagonalSymmetricMatrix *)Pass)->CholeskyBackSolve(Product, x);
	};

};

class EdgePreservingDecomposition{
public:
	EdgePreservingDecomposition(unsigned int width, unsigned int height);
	~EdgePreservingDecomposition();

	//Create an edge preserving blur of Source. Will create and return, or fill into Blur if not NULL. In place not ok.
	float *CreateBlur(float *Source, float Scale = 1.0f, float EdgeStopping = 1.4f, unsigned int Iterates = 20, float *Blur = NULL);

	/*Lowers global contrast while preserving or boosting local contrast. Can fill into Compressed. The smaller Compression
	the more compression is applied, with Compression = 1 giving no effect and above 1 the opposite effect. You can totally
	use Compression = 1 and play with DetailBoost for some really sweet unsharp masking.
	In place calculation to save memory (Source == Compressed) is totally ok. */
	float *CompressDynamicRange(float *Source, float Scale = 1.0f, float EdgeStopping = 1.4f, float Compression = 0.1f, float DetailBoost = 0.1f, unsigned int Iterates = 20, float *Compressed = NULL);

private:
	MultiDiagonalSymmetricMatrix *A;	//The equations are simple enough to not mandate a matrix class, but fast solution NEEDS a complicated preconditioner.
	unsigned int w, h, n;

	//Convenient access to the data in A.
	float *a0, *a_1, *a_w, *a_w_1, *a_w1;
};

//This function is for replacing 16 bit luminance in Lab (or Luv, or Yuv, or similar) with saturation control per Mantiuk's paper.
void ModifyCIELabLuminance(float *NewL, unsigned short *L, unsigned short *a, unsigned short *b, unsigned int npx, float SaturationControl = 1.0f);

