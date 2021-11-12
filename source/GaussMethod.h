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

	//Считывание кол-ва уравнений слау из файла
	int dimensionOfMatrixFromFile(std::string slaepath);

	//Считывание из файла матрицы произвольного размера
	double** slaeFromFile(std::string slaepath);

	//Считывание из слау матрицы коэффицентов
	double** aFromSlae(double** slae, int n);

	//Считывание из слау свободных членов
	double* bFromSlae(double** slae, int n);

	//Генерация матрицы коэффицентов размерности n
	double** randMatrixA(int n);

	//Генерация свободных членов размерности n
	double* randMatrixB(int n);

	//Метод Гаусса
	double* gauss(double** a, double* y, int n);

	//Невязка
	double discrepancy(double** a, double* b, double* x, int n);

	//Умножаем матрицу размерностью nxn на матрицу nx1
	double* multMatrix(double** a, double* b, int n);

	//Преобразование матрицы в обратную.
	double** reverseMatrix(double** matrix, int n);

	//Нахождение определителя матрицы
	double determinant(double** matrix, int n, double epsilon);

	//Метод прогонки
	void sweepMethod(double* v1, double* v2, double* v3, double* d, double* result, int n);

	//Клонирование матрицы
	double** matrixClone(double** matrix, int n);

	//Клонирование массива
	double* arrayClone(double* array, int n);

	//Нахождение зависимости времени выполнения от кол-ва уравнений
	std::map<int, double> checkAsymptotics(int count);
};