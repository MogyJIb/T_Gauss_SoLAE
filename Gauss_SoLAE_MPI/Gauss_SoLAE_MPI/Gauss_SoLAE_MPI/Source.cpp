//#include <mpi.h>
//#include "Logger.h"
//#include "GaussMethod.h"
//#include "GaussMethod_MPI.h"
//#include <string>
//#include <fstream>
//#include <iostream>
//
//int ProcNum; // Число доступных процессоров
//int ProcRank; // Ранг текущего процессора
//int *pParallelPivotPos; // Номера строк, которые были выбраны ведущими
//int *pProcPivotIter; // Номера итераций, на которых строки
// // процессора использовались в качестве ведущих
//void main(int argc, char* argv[]) {
//	double* pMatrix; // Матрица линейной системы
//	double* pVector; // Вектор правых частей линейной системы
//	double* pResult; // Вектор неизвестных
//	double *pProcRows; // Строки матрицы A
//	double *pProcVector; // Блок вектора b
//	double *pProcResult; // Блок вектора x
//	int Size; // Размер матрицы и векторов
//	int RowNum; // Количество строк матрицы
//	double start, finish, duration;
//	setvbuf(stdout, 0, _IONBF, 0);
//	MPI_Init(&argc, &argv);
//	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
//	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
//
//		if (ProcRank == 0)
//			printf("Параллельный метод Гаусса для решения систем линейных уравнений\n");
//	// Выделение памяти и инициализация данных
//	ProcessInitialization(pMatrix, pVector, pResult,
//		pProcRows, pProcVector, pProcResult, Size, RowNum);
//	// Распределение исходных данных
//	DataDistribution(pMatrix, pProcRows, pVector, pProcVector,
//		Size, RowNum);
//	// Выполнение параллельного алгоритма Гаусса
//	ParallelResultCalculation(pProcRows, pProcVector, pProcResult, Size,
//		RowNum);
//	// Сбор найденного вектора неизвестных на ведущем процессе
//	ResultCollection(pProcResult, pResult);
//
//	// Завершение процесса вычислений
//	ProcessTermination(pMatrix, pVector, pResult, pProcRows,
//		pProcVector, pProcResult);
//	MPI_Finalize();
//}////// Функция для параллельного выполнения метода Гаусса
//void ParallelResultCalculation(double* pProcRows,
//	double* pProcVector, double* pProcResult, int Size, int RowNum) {
//	ParallelGaussianElimination(pProcRows, pProcVector, Size,
//		RowNum);
//	ParallelBackSubstitution(pProcRows, pProcVector, pProcResult,
//		Size, RowNum);
//}//////// Функция для параллельного выполнения прямого хода метода Гаусса
//void ParallelGaussianElimination(double* pProcRows,
//	double* pProcVector, int Size, int RowNum) {
//	double MaxValue; // Значение ведущего элемента на процессоре
//	int PivotPos; // Положение ведущей строки в полосе линейной
//	// системы даннного процессора
//	// Структура для выбора ведущего элемента
//	struct { double MaxValue; int ProcRank; } ProcPivot, Pivot;
//	// pPivotRow используется для хранения ведущей строки матрицы и
//	// соответствующего элемента вектора b
//	double* pPivotRow = new double[Size + 1];
//	// Итерации прямого хода метода Гаусса
//	for (int i = 0; i < Size; i++) {
//
//		// Нахождение ведущей строки среди строк процесса
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
//		// Нахождение ведущего процесса (процесса, который содержит
//		// максимальное значение переменной MaxValue)
//		MPI_Allreduce(&ProcPivot, &Pivot, 1, MPI_DOUBLE_INT,
//			MPI_MAXLOC, MPI_COMM_WORLD);
//		// Рассылка ведущей строки
//		if (ProcRank == Pivot.ProcRank) {
//			pProcPivotIter[PivotPos] = i; // номер итерации
//			pParallelPivotPos[i] = pProcInd[ProcRank] + PivotPos;
//		}
//		MPI_Bcast(&pParallelPivotPos[i], 1, MPI_INT, Pivot.ProcRank,
//			MPI_COMM_WORLD);
//		if (ProcRank == Pivot.ProcRank) {
//			// заполнение ведущей строки
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
//}////// Функция для параллельного выполнения обратного хода метода Гаусса
//void ParallelBackSubstitution(double* pProcRows, double* pProcVector,
//	double* pProcResult, int Size, int RowNum) {
//	int IterProcRank; // Ранг процессора, который содержит ведущую строку
//	int IterPivotPos; // Положение ведущей строки в полосе процессора
//	double IterResult; // Вычисленное значение очередной неизвестной
//	double val;
//	// Итерации обратного хода метода Гаусса
//	for (int i = Size - 1; i >= 0; i--) {
//		// Вычисление ранга процессора, который содержит ведущую строку
//		FindBackPivotRow(pParallelPivotPos[i], Size, IterProcRank,
//			IterPivotPos);
//
//		// Вычисление значения неизвестной
//		if (ProcRank == IterProcRank) {
//			IterResult = pProcVector[IterPivotPos] / pProcRows[IterPivotPos*Size + i];
//			pProcResult[IterPivotPos] = IterResult;
//		}
//		// Рассылка значения очередной неизвестной
//		MPI_Bcast(&IterResult, 1, MPI_DOUBLE, IterProcRank, MPI_COMM_WORLD);
//		// Обновление вектора правых частей
//		for (int j = 0; j < RowNum; j++)
//			if (pProcPivotIter[j] < i) {
//				val = pProcRows[j*Size + i] * IterResult;
//				pProcVector[j] = pProcVector[j] - val;
//			}
//	}
//}//////// Функция для выделения памяти и инициализации исходных данных
//void ProcessInitialization(double* &pAMatrix, double* &pBMatrix,
//	double* &pCMatrix, double* &pAblock, double* &pBblock,
//	double* &pCblock, double* &pTemporaryAblock, int &Size,
//	int &BlockSize) {
//	if (ProcRank == 0) {
//		do {
//			printf("\nВведите размер матриц: ");
//			scanf("%d", &Size);
//
//			if (Size%GridSize != 0) {
//				printf("Размер матриц должен быть кратен размеру сетки! \n");
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