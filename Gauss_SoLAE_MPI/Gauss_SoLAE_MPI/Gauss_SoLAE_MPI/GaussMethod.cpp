#include "pch.h"
#include "GaussMethod.h"


GaussMethod::GaussMethod()
{
}


GaussMethod::~GaussMethod()
{
}

double* GaussMethod::solve(double **data, int size) {
 	this->columnCount = size + 1;
	this->rowCount = size;

	makeTriangle(data);
	return solveTriangle(data, UPPER);
}

void GaussMethod::makeTriangle(double **data) {
	for (int i = 0; i < rowCount - 1; i++) 
		for (int j = i + 1; j < rowCount; j++) {
			double coefficient = -data[j][i] / data[i][i];
			logger.logI(to_string(coefficient).append("   -coeff  ").append(to_string(data[i][i]).append("   -diagonal"))); //log
			for (int k = 0; k < columnCount; k++)
				data[j][k] += data[i][k] * coefficient;
		}
	
}


double* GaussMethod::solveTriangle(double **data, bool isLow) {
	if (columnCount < 10) //LOG
		for (int i = 0; i < rowCount; i++) {
			string to_log("");
			for (int j = 0; j < columnCount; j++) {
				to_log.append(to_string(data[i][j])).append(" ");
			}
			logger.logI(to_log); 
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