#include "pch.h"
#include "GaussMethod_MPI.h"
#include "mpi.h"
#include <math.h>


GaussMethod_MPI::GaussMethod_MPI(int processId, int processCount) {
	this->processId = processId;
	this->processCount = processCount;
}


GaussMethod_MPI::~GaussMethod_MPI()
{
}

double* GaussMethod_MPI::solve(double **data, int size) {
	this->columnCount = size + 1;
	this->rowCount = size;

	makeTriangle(data);
	return solveTriangle(data, UPPER);
}

void GaussMethod_MPI::makeTriangle(double **data) {
	int step = getThreadsStep(processCount); 

	for (int i = 0; i < rowCount - 1; i++) {
		int startInd = i + 1 + processId * step,
			endInd = i + 1 + (processId + 1) * step;
		if (endInd > rowCount) endInd = rowCount;
		int count = endInd - startInd;

		double *lineToAdd = data[i];
		MPI_Bcast(&(lineToAdd[0]), columnCount, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		for (int j = 0; j < rowCount; j++) {
			string to_log("");
			for (int k = 0; k < columnCount; k++) {
				to_log.append(to_string(data[j][k])).append(" ");
			}
			logger.logI(to_log
				.append(to_string(processId).append(" -processId data before ")));
		} //LOG

		logger.logI(string("scatter - ").append("startInd  ").append(to_string(startInd)));
		double **tempData = new double*[count];
		for (int j = 0; j < count; j++) {
			tempData[j] = new double[columnCount];
		}

		MPI_Scatter(data, count * columnCount, MPI_DOUBLE, tempData, count * columnCount, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		for (int j = 0; j < count; j++) {
			string to_log("");
			for (int k = 0; k < columnCount; k++) {
				to_log.append(to_string(tempData[j][k])).append(" ");
			}
			logger.logI(to_log
				.append(to_string(processId).append(" -processId tempdata ")));
		} //LOG


		for (int j = 0; j < count; j++) {
			double coefficient = -tempData[j][i] / lineToAdd[i];
			for (int k = 0; k < columnCount; k++)
				tempData[j][k] += lineToAdd[k] * coefficient;
		}

		for (int j = 0; j < count; j++) {
			string to_log("");
			for (int k = 0; k < columnCount; k++) {
				to_log.append(to_string(tempData[j][k])).append(" ");
			}
			logger.logI(to_log
				.append(to_string(processId).append(" -processId tempdata ")));
		} //LOG

		for (int j = 0; j < count; j++) {
			MPI_Gather(tempData, count * columnCount, MPI_DOUBLE, data, count * columnCount, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		}

		for (int j = 0; j < count; j++) 
			delete[] tempData[j];
		delete[] tempData;

	}
	
}

int GaussMethod_MPI::getThreadsStep(int threadCount) {
	return (threadCount > rowCount) 
				? 1
				: (int) floor(1.0 * rowCount / threadCount);
}


double* GaussMethod_MPI::solveTriangle(double **data, bool isLow) {
	if (processId != 0) return NULL;

	if (columnCount < 10) //LOG
		for (int i = 0; i < rowCount; i++) {
			string to_log("");
			for (int j = 0; j < columnCount; j++) {
				to_log.append(to_string(data[i][j])).append(" ");
			}
			logger.logI(to_log
				.append(to_string(processId).append(" -processId  ")));
		} //LOG

	double* result = new double[rowCount];
	for (int j = 0; j < rowCount; j++)
		result[j] = 0.0;

	double rowSum = 0.0;
	int i = !isLow ? rowCount - 1 : 0;

	while ((!isLow && i >= 0) || (isLow && i < rowCount)) {
		rowSum = 0.0;
		for (int j = 0; j < columnCount - 1; j++)
			rowSum += result[j] * data[i][j];
		rowSum = data[i][columnCount - 1] - rowSum;
		result[i] = data[i][i] == 0 ? 0 : rowSum / data[i][i];

		i = !isLow ? i - 1 : i + 1;
	}

	return result;
}