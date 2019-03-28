//#include <mpi.h>
//#include "Logger.h"
//#include "GaussMethod.h"
//#include "GaussMethod_MPI.h"
//#include <string>
//#include <fstream>
//#include <iostream>
//
//int ProcNum; // ����� ��������� �����������
//int ProcRank; // ���� �������� ����������
//int *pParallelPivotPos; // ������ �����, ������� ���� ������� ��������
//int *pProcPivotIter; // ������ ��������, �� ������� ������
// // ���������� �������������� � �������� �������
//void main(int argc, char* argv[]) {
//	double* pMatrix; // ������� �������� �������
//	double* pVector; // ������ ������ ������ �������� �������
//	double* pResult; // ������ �����������
//	double *pProcRows; // ������ ������� A
//	double *pProcVector; // ���� ������� b
//	double *pProcResult; // ���� ������� x
//	int Size; // ������ ������� � ��������
//	int RowNum; // ���������� ����� �������
//	double start, finish, duration;
//	setvbuf(stdout, 0, _IONBF, 0);
//	MPI_Init(&argc, &argv);
//	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
//	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
//
//		if (ProcRank == 0)
//			printf("������������ ����� ������ ��� ������� ������ �������� ���������\n");
//	// ��������� ������ � ������������� ������
//	ProcessInitialization(pMatrix, pVector, pResult,
//		pProcRows, pProcVector, pProcResult, Size, RowNum);
//	// ������������� �������� ������
//	DataDistribution(pMatrix, pProcRows, pVector, pProcVector,
//		Size, RowNum);
//	// ���������� ������������� ��������� ������
//	ParallelResultCalculation(pProcRows, pProcVector, pProcResult, Size,
//		RowNum);
//	// ���� ���������� ������� ����������� �� ������� ��������
//	ResultCollection(pProcResult, pResult);
//
//	// ���������� �������� ����������
//	ProcessTermination(pMatrix, pVector, pResult, pProcRows,
//		pProcVector, pProcResult);
//	MPI_Finalize();
//}////// ������� ��� ������������� ���������� ������ ������
//void ParallelResultCalculation(double* pProcRows,
//	double* pProcVector, double* pProcResult, int Size, int RowNum) {
//	ParallelGaussianElimination(pProcRows, pProcVector, Size,
//		RowNum);
//	ParallelBackSubstitution(pProcRows, pProcVector, pProcResult,
//		Size, RowNum);
//}//////// ������� ��� ������������� ���������� ������� ���� ������ ������
//void ParallelGaussianElimination(double* pProcRows,
//	double* pProcVector, int Size, int RowNum) {
//	double MaxValue; // �������� �������� �������� �� ����������
//	int PivotPos; // ��������� ������� ������ � ������ ��������
//	// ������� �������� ����������
//	// ��������� ��� ������ �������� ��������
//	struct { double MaxValue; int ProcRank; } ProcPivot, Pivot;
//	// pPivotRow ������������ ��� �������� ������� ������ ������� �
//	// ���������������� �������� ������� b
//	double* pPivotRow = new double[Size + 1];
//	// �������� ������� ���� ������ ������
//	for (int i = 0; i < Size; i++) {
//
//		// ���������� ������� ������ ����� ����� ��������
//		double MaxValue = 0;
//		for (int j = 0; j < RowNum; j++) {
//			if ((pProcPivotIter[j] == -1) &&
//				(MaxValue < fabs(pProcRows[j*Size + i]))) {
//				MaxValue = fabs(pProcRows[j*Size + i]);
//				PivotPos = j;
//			}
//		}
//		ProcPivot.MaxValue = MaxValue;
//		ProcPivot.ProcRank = ProcRank;
//		// ���������� �������� �������� (��������, ������� ��������
//		// ������������ �������� ���������� MaxValue)
//		MPI_Allreduce(&ProcPivot, &Pivot, 1, MPI_DOUBLE_INT,
//			MPI_MAXLOC, MPI_COMM_WORLD);
//		// �������� ������� ������
//		if (ProcRank == Pivot.ProcRank) {
//			pProcPivotIter[PivotPos] = i; // ����� ��������
//			pParallelPivotPos[i] = pProcInd[ProcRank] + PivotPos;
//		}
//		MPI_Bcast(&pParallelPivotPos[i], 1, MPI_INT, Pivot.ProcRank,
//			MPI_COMM_WORLD);
//		if (ProcRank == Pivot.ProcRank) {
//			// ���������� ������� ������
//				for (int j = 0; j < Size; j++) {
//					pPivotRow[j] = pProcRows[PivotPos*Size + j];
//				}
//			pPivotRow[Size] = pProcVector[PivotPos];
//		}
//		MPI_Bcast(pPivotRow, Size + 1, MPI_DOUBLE, Pivot.ProcRank,
//			MPI_COMM_WORLD);
//		ParallelEliminateColumns(pProcRows, pProcVector, pPivotRow,
//			Size, RowNum, i);
//	}
//}////// ������� ��� ������������� ���������� ��������� ���� ������ ������
//void ParallelBackSubstitution(double* pProcRows, double* pProcVector,
//	double* pProcResult, int Size, int RowNum) {
//	int IterProcRank; // ���� ����������, ������� �������� ������� ������
//	int IterPivotPos; // ��������� ������� ������ � ������ ����������
//	double IterResult; // ����������� �������� ��������� �����������
//	double val;
//	// �������� ��������� ���� ������ ������
//	for (int i = Size - 1; i >= 0; i--) {
//		// ���������� ����� ����������, ������� �������� ������� ������
//		FindBackPivotRow(pParallelPivotPos[i], Size, IterProcRank,
//			IterPivotPos);
//
//		// ���������� �������� �����������
//		if (ProcRank == IterProcRank) {
//			IterResult = pProcVector[IterPivotPos] / pProcRows[IterPivotPos*Size + i];
//			pProcResult[IterPivotPos] = IterResult;
//		}
//		// �������� �������� ��������� �����������
//		MPI_Bcast(&IterResult, 1, MPI_DOUBLE, IterProcRank, MPI_COMM_WORLD);
//		// ���������� ������� ������ ������
//		for (int j = 0; j < RowNum; j++)
//			if (pProcPivotIter[j] < i) {
//				val = pProcRows[j*Size + i] * IterResult;
//				pProcVector[j] = pProcVector[j] - val;
//			}
//	}
//}//////// ������� ��� ��������� ������ � ������������� �������� ������
//void ProcessInitialization(double* &pAMatrix, double* &pBMatrix,
//	double* &pCMatrix, double* &pAblock, double* &pBblock,
//	double* &pCblock, double* &pTemporaryAblock, int &Size,
//	int &BlockSize) {
//	if (ProcRank == 0) {
//		do {
//			printf("\n������� ������ ������: ");
//			scanf("%d", &Size);
//
//			if (Size%GridSize != 0) {
//				printf("������ ������ ������ ���� ������ ������� �����! \n");
//			}
//		} while (Size%GridSize != 0);
//	}
//	MPI_Bcast(&Size, 1, MPI_INT, 0, MPI_COMM_WORLD);
//	BlockSize = Size / GridSize;
//	pAblock = new double[BlockSize*BlockSize];
//	pBblock = new double[BlockSize*BlockSize];
//	pCblock = new double[BlockSize*BlockSize];
//	pTemporaryAblock = new double[BlockSize*BlockSize];
//	for (int i = 0; i < BlockSize*BlockSize; i++) {
//			pCblock[i] = 0;
//	}
//	if (ProcRank == 0) {
//		pAMatrix = new double[Size*Size];
//		pBMatrix = new double[Size*Size];
//		pCMatrix = new double[Size*Size];
//		RandomDataInitialization(pAMatrix, pBMatrix, Size);
//	}
//}