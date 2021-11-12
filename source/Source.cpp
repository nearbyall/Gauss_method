#include <iostream>
#include <fstream>

#include "GaussMethod.h"

int main(int argc, char* argv[], char* envp[]) 
{
	setlocale(LC_ALL, "Russian");

	GaussMethod gaussMethod;

	int n = gaussMethod.dimensionOfMatrixFromFile("gauss_41.lb1");
	double** slae = gaussMethod.slaeFromFile("gauss_41.lb1");
	double** A = gaussMethod.aFromSlae(slae, n);
	double* B = gaussMethod.bFromSlae(slae, n);
	double** aTemp = gaussMethod.matrixClone(A, n);
	double* bTemp = gaussMethod.arrayClone(B, n);
	double* x = gaussMethod.gauss(A, B, n);
	double discrepancy = gaussMethod.discrepancy(aTemp, bTemp, x, n);
	//std::map<int, double> table = gaussMethod.checkAsymptotics(100);

	return 1;
}