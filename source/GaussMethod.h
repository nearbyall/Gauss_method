#pragma once

#include <iostream>
#include <fstream>
#include <map>
#include <ctime>
#include <chrono>
#include <iomanip>

class GaussMethod
{

public:

	//���������� ���-�� ��������� ���� �� �����
	int dimensionOfMatrixFromFile(std::string slaepath);

	//���������� �� ����� ������� ������������� �������
	double** slaeFromFile(std::string slaepath);

	//���������� �� ���� ������� ������������
	double** aFromSlae(double** slae, int n);

	//���������� �� ���� ��������� ������
	double* bFromSlae(double** slae, int n);

	//��������� ������� ������������ ����������� n
	double** randMatrixA(int n);

	//��������� ��������� ������ ����������� n
	double* randMatrixB(int n);

	//����� ������
	double* gauss(double** a, double* y, int n);

	//�������
	double discrepancy(double** a, double* b, double* x, int n);

	//�������� ������� ������������ nxn �� ������� nx1
	double* multMatrix(double** a, double* b, int n);

	//�������������� ������� � ��������.
	double** reverseMatrix(double** matrix, int n);

	//���������� ������������ �������
	double determinant(double** matrix, int n, double epsilon);

	//����� ��������
	void sweepMethod(double* v1, double* v2, double* v3, double* d, double* result, int n);

	//������������ �������
	double** matrixClone(double** matrix, int n);

	//������������ �������
	double* arrayClone(double* array, int n);

	//���������� ����������� ������� ���������� �� ���-�� ���������
	std::map<int, double> checkAsymptotics(int count);
};