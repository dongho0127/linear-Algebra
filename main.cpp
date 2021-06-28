#include <iostream>
#include <cmath>
#include <cstdlib>

#define MAT_SIZE 100
using namespace std;

class LeastError
{
public:
	LeastError(int POLYrank);

	void getCoor(void);
	//좌표 입력받음
	void MainCalculator(void);
	//전체 계산
	void Graph(void);
	//그래프화 함

private:
	int POLYRANK, CoorS;
	double x[MAT_SIZE], y[MAT_SIZE];
	double polyCOEF[100];
	double A[MAT_SIZE][MAT_SIZE], Q[MAT_SIZE][MAT_SIZE], R[MAT_SIZE][MAT_SIZE];
	double MatCash[MAT_SIZE][MAT_SIZE], MatCash2[MAT_SIZE][MAT_SIZE], VecCash[MAT_SIZE], VecCash2[MAT_SIZE];
	double Error[MAT_SIZE];


	void MatrixMul(double A[MAT_SIZE][MAT_SIZE], double B[MAT_SIZE][MAT_SIZE], int Arow, int Acol, int Brow, int Bcol);
	//두개의 행렬을 곱해서(AB) Matcash에 저장하는 함수

	void MatrixMul(double A[][MAT_SIZE], double vector[], int Arow, int Acol);
	//벡터와 행렬의 곱을 구한 후 값을 Vecash에 저장한다

	void MatrixINV(double A[][MAT_SIZE], int rank);
	//역행렬을 계산하여 Matcash에 저장하는 함수

	void matCpy(double newmat[][MAT_SIZE], double oldmat[][MAT_SIZE], int x, int y, int r1, int r2, int c1, int c2);
	//주어진 행, 열의 범위에 해당하는 행렬을 새로운 행렬의(x,y)에 복사한다

	void Swab(double mat1[][MAT_SIZE], int r1, int r2, int n, int m);
	//받은 행렬의 두 행을 바꾼다

	void Add(double mat2[][MAT_SIZE], int r1, double k, int r2, int n, int m);
	//받은 행렬에 Ri+k*Rj를 시행한다

	void showMat(double mat3[][MAT_SIZE], int r1, int r2, int c1, int c2);
	//받은 행렬의 i~j행, x~y행의 행렬만 보여준다

	void Echelon(double mat4[][MAT_SIZE], int& rowcash, int& colcash, int n, int m);
	//행렬을 행줄임 한다.

	void simple(double mat5[][MAT_SIZE], int n, int m);
	//행줄임된 행렬의 기약형태를 계산한다.

	void QRfactorization(double A[][MAT_SIZE], int  rowS, int colS);
	//A에 대해 QR분해를 한 후 각각의 행렬에 저장한다

	void normalization(double A[][MAT_SIZE], int rowS, int whatcol);
	//행렬을 정규화한다

	void Transpose(double A[][MAT_SIZE], int rowS, int colS);
	//A의 전치화된 행렬을 Matcash 에 저장한다

};

int main(void)
{
	cout.setf(ios::fixed);
	cout.setf(ios::showpoint);
	cout.precision(4);

	int polyrank;
	while (1)
	{
		cout << "what is the polynomial's rank? : ";
		cin >> polyrank;

		LeastError LE(polyrank);
		LE.getCoor();
		LE.MainCalculator();
		cout << endl;
	}
	return 0;
}

void LeastError:: Swab(double mat1[][MAT_SIZE], int r1, int r2, int n, int m)
{
	double cash;
	for (int i = 1; i <= m; i++)
	{
		cash = mat1[r1][i];
		mat1[r1][i] = mat1[r2][i];
		mat1[r2][i] = cash;
	}
	return;
}

void LeastError:: Add(double mat2[][MAT_SIZE], int r1, double k, int r2, int n, int m)

{
	for (int i = 1; i <= m; i++)
	{
		mat2[r1][i] += (mat2[r2][i] * k);
	}
	return;
}

void LeastError:: showMat(double mat3[][MAT_SIZE], int r1, int r2, int c1, int c2)
{
	for (int i = r1; i <= r2; i++)
	{
		for (int j = c1; j <= c2; j++)
		{
			cout << mat3[i][j] << "    ";
		}
		cout << "\n" << endl;
	}
	cout << "\n";
	return;
}

void LeastError::Echelon(double mat4[][MAT_SIZE], int& rowcash, int& colcash, int n, int m)
{
	double firstNum;
	int i, j, firstRow;
	if ((rowcash == n) || (colcash == m + 1))
	{
		return;
	}

	for (i = rowcash; i <= n; i++)
	{
		if (mat4[i][colcash] != 0)
		{
			firstRow = i;
			firstNum = mat4[i][colcash];
			Swab(mat4, firstRow, rowcash, n, m);
			for (j = rowcash + 1; j <= n; j++)
			{
				Add(mat4, j, static_cast<double>(-1 * mat4[j][colcash] / firstNum), rowcash, n, m);
			}
			rowcash++;
			colcash++;
			return Echelon(mat4, rowcash, colcash, n, m);

		}
	}
}

void LeastError::MatrixINV(double A[][MAT_SIZE], int rank)
	{
		int rowcash = 1, colcash = 1;
		for (int i = 1; i <= rank; i++)
		{
			for (int j = rank + 1; j <= 2 * rank; j++)
			{
				A[i][j] = 0;
			}
			A[i][i + rank] = 1;
		}
		Echelon(A, rowcash, colcash, rank, 2 * rank);
		simple(A, rank, 2 * rank);
		matCpy(MatCash, A, 1, 1, 1, rank, rank + 1, 2 * rank);
	}

void LeastError::matCpy(double newmat[][MAT_SIZE], double oldmat[][MAT_SIZE], int x, int y, int r1, int r2, int c1, int c2)
{
	for (int i = r1; i <= r2; i++)
	{
		for (int j = c1; j <= c2; j++)
		{
			newmat[x + i - r1][y + j - c1] = oldmat[i][j];
		}
	}
}

void LeastError::MatrixMul(double A[][MAT_SIZE], double vector[], int Arow, int Acol)
{
	double sumcash, veccash[MAT_SIZE];
	for (int i = 1; i <= Arow; i++)
	{
		sumcash = 0;
		for (int j = 1; j <= Acol; j++)
		{
			sumcash = sumcash + A[i][j] * vector[j];
		}
		veccash[i] = sumcash;
	}
	for (int i = 1; i <= Arow; i++)
	{
		VecCash[i] = veccash[i];
	}
}

void LeastError::simple(double mat5[][MAT_SIZE], int n, int m)
{
	double firstcol;
	for (int i = 1; i <= n; i++)
	{
		firstcol = mat5[i][i];
		for (int j = 1; j <= m; j++)
		{
			mat5[i][j] = mat5[i][j] / firstcol;
		}
	}
	for (int i = 2; i <= n; i++)
	{
		for (int j = 1; j < i; j++)
		{
			Add(mat5, j, -1 * mat5[j][i], i, n, m);
		}
	}
	return;
}

void LeastError::MatrixMul(double A1[MAT_SIZE][MAT_SIZE], double B[MAT_SIZE][MAT_SIZE], int Arow, int Acol, int Brow, int Bcol)
{
	double cash;
	for (int i = 1; i <= Arow; i++)
	{
	
		for (int j = 1; j <= Bcol; j++)
		{
			cash = 0;
			for (int k = 1; k <= Acol; k++)
			{
				cash += A1[i][k] * B[k][j];
			}
			MatCash[i][j] = cash;
		}
	}
}

void LeastError::normalization(double A[][MAT_SIZE], int rowS, int whatcol)
{
	double colNormCASH=0;
	for (int j = 1; j <= rowS; j++)
	{
		colNormCASH += A[j][whatcol] * A[j][whatcol];
	}
	colNormCASH = sqrt(colNormCASH);
	for (int i = 1; i <= rowS; i++)
	{
		A[i][whatcol] /= colNormCASH;
	}
	
}

void LeastError::Transpose(double A[][MAT_SIZE], int rowS, int colS)
{
	for (int i = 1; i <= rowS; i++)
	{
		for (int j = 1; j <= colS; j++)
		{
			MatCash[j][i] = A[i][j];
		}
	}
	
}

void LeastError::QRfactorization(double A[][MAT_SIZE], int  rowS, int colS)
{
	matCpy(Q, A, 1, 1, 1, rowS, 1, colS);
	double innerPr;
	for (int i = 1; i <= colS; i++)
	{
		for (int k = 1; k <= i - 1; k++)
		{	
			innerPr = 0;
			for (int j = 1; j <= rowS; j++)
			{
				innerPr += Q[j][k] * A[j][i];
			}
			for (int j = 1; j <= rowS; j++)
			{
				Q[j][i] -= innerPr * Q[j][k];
			}
		}
		normalization(Q, rowS, i);
	}
	Transpose(Q, rowS, colS);
	matCpy(MatCash2, MatCash, 1, 1, 1, colS, 1, rowS);
	MatrixMul(MatCash2, A, colS, rowS, rowS, colS);
	matCpy(R, MatCash, 1, 1, 1, colS, 1, colS);
}

LeastError::LeastError(int POLYrank) : POLYRANK(POLYrank)
{
	if (POLYRANK <= 0)
	{
		cout << "the rank of the polynomial must be greather than 0!" << endl;
		exit(1);
	}
}

void LeastError::getCoor(void)
{
	cout << "How many coordinates do you want to input? : ";
	cin >> CoorS;
	cout << endl;
	for (int i = 1; i <= CoorS; i++)
	{
		cout << "the " << i << "th coordinate: ";
		cin >> x[i] >> y[i];
	}
}

void LeastError::MainCalculator(void)
{
	for (int i = 1; i <= POLYRANK + 1; i++)
	{
		for (int j = 1; j <= CoorS; j++)
		{
			A[j][i] = pow(x[j], i - 1);
		}
	}

	QRfactorization(A, CoorS, POLYRANK + 1);
	Transpose(Q, CoorS, POLYRANK + 1);
	MatrixMul(MatCash, y, POLYRANK+1, CoorS);
	MatrixINV(R, POLYRANK + 1);
	MatrixMul(MatCash, VecCash, POLYRANK + 1, POLYRANK + 1);
	cout << endl;

	cout << "y= ";
	for (int i = 1; i <= POLYRANK + 1; i++)
	{
		cout << "(" << VecCash[i] << ")";
		if (i != 1)
		{
			cout << "x^"<<i-1;
		}
		if (i != POLYRANK + 1)
		{
			cout << " + ";
		}
	}
	cout << endl<<endl;

	//이제 에러계산->Error[0]에 에러의 크기 들어감!

	for (int i = 1; i <= CoorS; i++)
	{
		VecCash2[i] = VecCash[i];
	}
	MatrixMul(A, VecCash2, CoorS, POLYRANK + 1);
	Error[0] = 0;
	for (int i = 1; i <= CoorS; i++)
	{
		Error[i] = y[i] - VecCash[i];
		Error[0] += Error[i] * Error[i];
	}
	Error[0] = sqrt(Error[0]);
	cout << "The norm of the error vector is " << Error[0] << endl << endl;
}

void LeastError::Graph(void)
{
	;
}