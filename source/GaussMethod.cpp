#include "GaussMethod.h"

int GaussMethod::dimensionOfMatrixFromFile(std::string slaepath)
{
	std::ifstream file;
	file.open(slaepath);
	if (file.is_open())
	{
		//Если открытие файла прошло успешно

		//Подсчитаем кол-во чисел в файле
		int count = 0;
		int temp;

		while (!file.eof())
		{
			file >> temp;
			count++;
		}

		//Вначале переведем каретку в потоке в начало файла
		file.seekg(0, std::ios::beg);
		file.clear();

		//Подсчитаем кол-во отступов в строке
		int count_space = 0;
		char symbol;

		while (!file.eof())
		{
			file.get(symbol);
			if (symbol == '\t') count_space++;
			if (symbol == '\n') break;
		}

		//Теперь мы знаем сколько чисел в файле и сколько пробелов в первой строке.
		//Теперь можем считать матрицу.
		int n = count / (count_space + 1);
		return n;
		file.close();
	}
	else
	{
		//Если открытие файла прошло не успешно
		std::cout << "Файл не найден.";
		return 0;
	}
}

double** GaussMethod::slaeFromFile(std::string slaepath)
{
	//Считаем СЛАУ из файла в двумерный динамический массив
	std::ifstream file;
	file.open(slaepath);
	if (file.is_open())
	{
		//Если открытие файла прошло успешно

		//Подсчитаем кол-во чисел в файле
		int count = 0;
		int temp;

		while (!file.eof())
		{
			file >> temp;
			count++;
		}

		//Вначале переведем каретку в потоке в начало файла
		file.seekg(0, std::ios::beg);
		file.clear();

		//Подсчитаем кол-во отступов в строке
		int count_space = 0;
		char symbol;

		while (!file.eof())
		{
			file.get(symbol);
			if (symbol == '\t') count_space++;
			if (symbol == '\n') break;
		}

		//Опять переходим в потоке в начало файл
		file.seekg(0, std::ios::beg);
		file.clear();

		//Теперь мы знаем сколько чисел в файле и сколько пробелов в первой строке.
		//Теперь можем считать матрицу.
		int n = count / (count_space + 1);
		int m = count_space + 1;
		double** slae = new double* [n];
		for (int i = 0; i < n; i++) slae[i] = new double[m];

		//Считаем матрицу из файла
		for (int i = 0; i < n; i++)
			for (int j = 0; j < m; j++)
				file >> slae[i][j];

		return slae;

		file.close();
	}
	else
	{
		//Если открытие файла прошло не успешно
		std::cout << "Файл не найден.";
		return nullptr;
	}

}

double** GaussMethod::aFromSlae(double** slae, int n) 
{
	double** A = new double* [n];
	for (int i = 0; i < n; i++) A[i] = new double[n];
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
			A[i][j] = slae[i][j];
	}
	return A;
}

double* GaussMethod::bFromSlae(double** slae, int n)
{
	double* B = new double[n];
	for (int i = 0; i < n; i++) {
		B[i] = slae[i][n];
	}
	return B;
}

double** GaussMethod::randMatrixA(int n)
{

	double** a = new double* [n];
	for (int i = 0; i < n; i++)
		a[i] = new double[n];

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			a[i][j] = (rand() % 100) + 1;

	return a;

}

double* GaussMethod::randMatrixB(int n)
{
	double* b = new double[n];
	for (int i = 0; i < n; i++)
		b[i] = (rand() % 100) + 1;
	return b;
}

double* GaussMethod::gauss(double** a, double* y, int n)
{
	double* x, max;
	int k, index;
	const double eps = 0.000001;  // точность
	x = new double[n];
	k = 0;
	std::ofstream fout("slae.txt");

	//Вывод в файл матрицы до преобразования
	fout << "Матрица до преобразования:" << std::endl;
	for (int i = 0; i < n; i++) {
		fout << "|";
		for (int j = 0; j < n; j++) {
			fout << std::setw(8) << a[i][j] << "\t";
		}
		fout << "|" << "\t" << std::setw(3) << y[i] << "|";
		fout << std::endl;
	}
	fout.close();

	while (k < n)
	{
		// Поиск строки с максимальным a[i][k]
		max = abs(a[k][k]);
		index = k;
		for (int i = k + 1; i < n; i++)
		{
			if (abs(a[i][k]) > max)
			{
				max = abs(a[i][k]);
				index = i;
			}
		}
		// Перестановка строк
		if (max < eps)
		{
			// нет ненулевых диагональных элементов
			std::cout << "Решение получить невозможно из-за нулевого столбца ";
			std::cout << index << " матрицы A" << std::endl;
			return 0;
		}
		for (int j = 0; j < n; j++)
		{
			double temp = a[k][j];
			a[k][j] = a[index][j];
			a[index][j] = temp;
		}
		double temp = y[k];
		y[k] = y[index];
		y[index] = temp;
		// Нормализация уравнений
		for (int i = k; i < n; i++)
		{
			double temp = a[i][k];
			if (abs(temp) < eps) continue; // для нулевого коэффициента пропустить
			for (int j = 0; j < n; j++)
				a[i][j] = a[i][j] / temp;
			y[i] = y[i] / temp;
			if (i == k)  continue; // уравнение не вычитать само из себя
			for (int j = 0; j < n; j++)
				a[i][j] = a[i][j] - a[k][j];
			y[i] = y[i] - y[k];
		}
		k++;
	}

	//Вывод в файл матрицы после преобразования
	fout.open("upper_triangle_matrix.txt");
	fout << "Матрица после преобразования:" << std::endl;
	for (int i = 0; i < n; i++) {
		fout << "|";
		for (int j = 0; j < n; j++) {
			fout << std::setw(8) << a[i][j] << "\t";
		}
		fout << "|" << "\t" << std::setw(8) << y[i] << "|";
		fout << std::endl;
	}
	fout.close();
	
	// Обратная подстановка и вывод в файл найденных неизвестных
	fout.open("result.txt");
	for (k = n - 1; k >= 0; k--)
	{
		x[k] = y[k];
		for (int i = 0; i < k; i++)
			y[i] = y[i] - a[i][k] * x[k];
		fout << "X" << k + 1 << "= " << x[k] << "\t";
	}
	fout.close();
	return x;

}

double GaussMethod::discrepancy(double** a, double* b, double* x, int n)
{
	double* dx = new double[n];

	dx = multMatrix(a, x, n);
	for (int i = 0; i < n; i++) {        //находим dx как  Ax-b=dx
		dx[i] -= b[i];
		//cout << dx[i] << endl;
	}
	double norma;
	double norma1 = 0;
	double norma2 = 0;
	for (int i = 0; i < n; i++) {           // Находим норму 1
		if (abs(dx[i]) >= norma1) {
			norma1 = abs(dx[i]);
		}
	}
	for (int i = 0; i < n; i++) {          //Находим норму 2
		norma2 += abs(dx[i]);
	}
	if (norma1 >= norma2) {
		norma = norma1;
	}
	else {
		norma = norma2;
	}
	double p = 0.00000000000001;
	norma = norma * p;

	//Вывод в файл невязки
	std::ofstream fout("discrepancy.txt");
	fout << "||Ax-B||=";                        //Находим максимальную норму для оценки погрешности 
	fout << norma << std::endl;
	fout.close();

	delete[] dx;

	return norma;
}

double* GaussMethod::multMatrix(double** a, double* b, int n)
{
	double* c = new double[n];

	for (int j = 0; j < n; j++)
	{
		c[j] = 0;
		for (int k = 0; k < n; k++)
			c[j] += a[j][k] * b[k];                  //Умножаем матрицу n x n на матрицу n x 1 

	}
	return c;
}

double** GaussMethod::reverseMatrix(double** matrix, int n)
{
	int i, j, k;
	//создание единичной матрицы 
	double** mob = new double* [n];
	for (i = 0; i < n; i++)
	{
		mob[i] = new double[n];
		for (j = 0; j < n; j++)mob[i][j] = 0;
		mob[i][i] = 1;
	}
	//прямой ход методом Гаусса
	double a, b;
	for (i = 0; i < n; i++)
	{
		a = matrix[i][i];
		for (j = i + 1; j < n; j++)
		{
			b = matrix[j][i];
			for (k = 0; k < n; k++)
			{
				matrix[j][k] = matrix[i][k] * b - matrix[j][k] * a;
				mob[j][k] = mob[i][k] * b - mob[j][k] * a;
			}
		}
	}
	return mob;
}

double GaussMethod::determinant(double** matrix, int n, double epsilon)
{
	//Предполагается, что матрица квадратная
	int pivot_index = -1;
	double pivot_value = 0;
	double determinant = 1;

	for (int i = 0; i < n; i++)
	{
		for (int j = i; j < n; j++)
		{
			if (abs(matrix[j][i]) > pivot_value)
			{
				pivot_index = j;
				pivot_value = abs(matrix[j][i]);
			}
		}

		//Если опорный элемент равен нулю (эпсилон для сброса погрешности)
		if (pivot_value < epsilon)
		{
			//Матрица вырождена
			return 0;
		}

		if (pivot_index != i)
		{
			//Обменяем строки местами
			double* temp = matrix[i];
			matrix[i] = matrix[pivot_index];
			matrix[pivot_index] = temp;
			determinant *= -1;
		}

		for (int j = i + 1; j < n; j++)
		{
			if (matrix[j][i] != 0)
			{
				double multiplier = 1 / matrix[i][i] * matrix[j][i];

				for (int k = i; k < n; k++)
				{
					matrix[j][k] -= matrix[i][k] * multiplier;
				}
			}
		}

		determinant *= matrix[i][i];
	}

	return determinant;
}

void GaussMethod::sweepMethod(double* v1, double* v2, double* v3, double* d, double* result, int n)
{
	double* y = new double[n];
	double* alpha = new double[n];
	double* beta = new double[n];

	y[0] = v2[0]; //Нулевые коэффиценты
	alpha[0] = -1 * v3[0] / y[0];
	beta[0] = d[0] / y[0];

	for (int i = 1; i < n - 1; i++)
	{
		y[i] = v2[i] + (v1[i] * alpha[i - 1]);
		alpha[i] = -1 * v3[i] / y[i];
		beta[i] = (d[i] - (v1[i] * beta[i - 1])) / y[i];
	}
	y[n - 1] = v2[n - 1] + (v1[n - 1] * alpha[n - 2]);
	beta[n - 1] = (d[n - 1] - (v2[n - 1] * beta[n - 2])) / y[n - 1];


	result[n - 1] = beta[n - 1];
	for (int i = n - 2; i > -1; i--)
	{
		result[i] = (alpha[i] * result[i + 1]) + beta[i];
	}

	delete[] alpha;
	delete[] beta;
	delete[] y;
}

double** GaussMethod::matrixClone(double** matrix, int n)
{
	double** newMatrix = new double* [n];
	for (int i = 0; i < n; i++) newMatrix[i] = new double[n];
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			newMatrix[i][j] = matrix[i][j];
	return newMatrix;
}

double* GaussMethod::arrayClone(double* array, int n)
{
	double* newArray = new double[n];
	for (int i = 0; i < n; i++) newArray[i] = array[i];
	return newArray;
}

std::map<int, double> GaussMethod::checkAsymptotics(int count)
{
	std::map<int, double> table;

	std::ofstream fout1("dimension.txt");
	std::ofstream fout2("complexity.txt");

	for (int i = 5; i < count + 5; i++) {
		double** a = randMatrixA(i);
		double* b = randMatrixB(i);

		auto t1 = std::chrono::high_resolution_clock::now();	

		double* result = gauss(a, b, i);

		auto t2 = std::chrono::high_resolution_clock::now();
		auto duration = (std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1)).count();

		table.insert(std::pair<int, int>(i, duration));

		fout1 << i << "\n";
		fout2 << duration << "\n";

		for (int j = 0; j < i; j++) {
			delete[] a[j];
		}
		delete[] a;
		delete[] b;
	}

	fout1.close();
	fout2.close();

	return table;
}