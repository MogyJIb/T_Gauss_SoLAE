#define _CRT_SECURE_NO_WARNINGS


#include <mpi.h>
#include "Logger.h"
#include "GaussMethod.h"
#include "GaussMethod_MPI.h"
#include <string>
#include <fstream>
#include <iostream>

using namespace std;

void checkArgs(int argc, char *argv[]);
void solve();
void solve_MPI(int argc, char *argv[]);


static bool USE_MPI = true;
static char inputFileName[256];
static char outputFileName[256];
static string PATH("");

int main(int argc, char *argv[])
{
	checkArgs(argc, argv);

	if (USE_MPI)
		solve_MPI(argc, argv);
	else
		solve();

	return 0;
}



void checkArgs(int argc, char *argv[]) {

	if (argc > 1) {
		PATH.append(argv[1]);
	}
	if (!PATH.empty()) {
		strcpy(inputFileName, string(PATH).append("input.txt").data());
		strcpy(outputFileName, string(PATH).append("output.txt").data());
	}
	else {
		strcpy(inputFileName, "input.txt");
		strcpy(outputFileName, "output.txt");
	}

	if (argc > 2) {
		//if (strcmp(argv[2], "1") == 0) USE_MPI = true;
	}
}


void solve() {

	Logger logger(PATH);

	logger.logI("\r\n\r\n\r\n\r\n");
	logger.logI(PATH); //log
	logger.logI(inputFileName); //log
	logger.logI(outputFileName); //log

	double **matrix;
	int size;

	logger.logI("read input"); //log

	ifstream fin;
	fin.open(inputFileName);
	fin >> size;

	matrix = new double*[size];
	for (int i = 0; i < size; i++)
		matrix[i] = new double[size + 1];

	for (int i = 0; i < size; i++) {
		string to_log("");
		for (int j = 0; j < size + 1; j++) {
			fin >> matrix[i][j];
			to_log.append(to_string(matrix[i][j])).append(" ");
		}
		if (size < 10) logger.logI(to_log); //log
	}

	fin.close();

	logger.logI("input read successful"); //log

	logger.logI(
		to_string(size).append("  -size")
	); //log


	double* result = new double[size];

	logger.logI("solve matrix"); //log

	try {
		GaussMethod gauss;
		result = gauss.solve(matrix, size);
	}
	catch (...) {
		logger.logI("exception in solving integral"); //log
	}

	logger.logI("solved"); //log


	logger.logI("write output"); //log

	ofstream fout;
	fout.open(outputFileName, ofstream::out);
	for (int i = 0; i < size - 1; i++)
		fout << result[i] << " ";
	fout << result[size - 1];
	fout.close();

	logger.logI("output write successful"); //log
}


void solve_MPI(int argc, char *argv[]) {
	int processId;
	int processCount;
	static int master = 0;
	int ierr;

	Logger logger(PATH);

	logger.logI("\r\n\r\n\r\n\r\n");
	logger.logI(PATH); //log
	logger.logI(inputFileName); //log
	logger.logI(outputFileName); //log

	ierr = MPI_Init(&argc, &argv);

	if (ierr != 0)
	{
		logger.logI("can't init mpi, fatal error");
		return;
	}

	MPI_Comm_size(MPI_COMM_WORLD, &processCount);
	MPI_Comm_rank(MPI_COMM_WORLD, &processId);

	logger.logI(to_string(processId).append("   processId ")); //log
	logger.logI(to_string(processCount).append("   processCount ")); //log

	double **matrix;
	int size;

	logger.logI("read input"); //log

	ifstream fin;
	fin.open(inputFileName);
	fin >> size;

	matrix = new double*[size];
	for (int i = 0; i < size; i++)
		matrix[i] = new double[size + 1];
	
	if (processId == 0) {
		for (int i = 0; i < size; i++) {
			string to_log("");
			for (int j = 0; j < size + 1; j++) {
				fin >> matrix[i][j];
				to_log.append(to_string(matrix[i][j])).append(" ");
			}
			if (size < 10) logger.logI(to_log); //log
		}

		logger.logI("input read successful"); //log

		
		logger.logI(
			to_string(size).append("  -size")
		); //log
	}
	fin.close();

	MPI_Bcast(&(matrix[0][0]), size * (size + 1), MPI_DOUBLE, 0, MPI_COMM_WORLD);

	double* result = new double[size];

	logger.logI("solve matrix"); //log

	try {
		GaussMethod_MPI gauss(processId, processCount);
		result = gauss.solve(matrix, size);
	}
	catch (...) {
		logger.logI("exception in solving integral"); //log
	}

	

	if (processId == 0) {
		logger.logI("solved"); //log
		logger.logI("write output"); //log

		ofstream fout;
		fout.open(outputFileName, ofstream::out);
		for (int i = 0; i < size - 1; i++)
			fout << result[i] << " ";
		fout << result[size - 1];
		fout.close();

		logger.logI("output write successful"); //log
	}
	

	MPI_Finalize();
}