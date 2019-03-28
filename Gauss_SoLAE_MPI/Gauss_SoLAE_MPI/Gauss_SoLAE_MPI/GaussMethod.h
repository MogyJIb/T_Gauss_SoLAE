#pragma once
#include "Logger.h"
#include <string>

using namespace std;

class GaussMethod
{
private:
	Logger logger;
	const bool LOW = true;
	const bool UPPER = !LOW;
	void makeTriangle(double **data);
	double* solveTriangle(double **data, bool isLow);

public:
	int columnCount;
	int rowCount;

	double* solve(double **data, int size);

	GaussMethod();
	~GaussMethod();
};



//3.0 2.0 3.0 1.0
//4.0 4.0 3.0 1.0
//1.0 4.0 4.0 2.0

//-2f / 11.0, 1f / 11.0, 5f / 11.0
// -0.18181818   0.0909090909   04545454545