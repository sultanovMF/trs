/// TODO:
//* Метод SOR
//* Точный метод (Гаусс, LU, метод сопряженных градиентов)
#pragma once
#include <span>
#include <vector>
#include <algorithm>
#include <stdexcept>

const double EPSILON = 0.000001;

/*
	Solve the matrix in the following form:
	| -B1  C1 0  ... | D1 // A1 = 0
	|  A2 -B2 C2 ... | D2
	| .............. | 
	| ....... An -Bn | Dn // Cn = 0
*/

std::vector<double> solve_tridiagonal(
	const std::span<double> A,
	const std::span<double> B,
	const std::span<double> C,
	const std::span<double> D)
{
	if (!(A.size() == B.size()  && B.size()  == C.size()  && C.size()  == D.size())) {
		throw std::invalid_argument("Tapes must be the same size!");
	}

	const int n = D.size();

	std::vector<double> P(n - 1);
	std::vector<double> Q(n);
	std::vector<double> X(n);

	// forward
	P[0] = C[0] / B[0];
	Q[0] = D[0] / B[0];
	for (int i = 1; i < n; ++i)
	{
		if (i < n - 1)
			P[i] = C[i] / (B[i] - A[i] * P[i - 1]);
		Q[i] = (D[i] - A[i] * Q[i - 1]) / (B[i] - A[i] * P[i - 1]);
	}

	// backward
	X[n - 1] = Q[n - 1];
	for (int i = n - 2; i >= 0; --i)
	{
		X[i] = Q[i] - P[i] * X[i + 1];
	}

	return X;
}

static void diag3_method(std::vector<double> a, std::vector<double> b, std::vector<double> c, std::vector<double> d, int n, std::vector<double>& x) {
	/*
	* a_i*x_{i-1} - b_i*x_i + c_i*x_{i+1} = d_i
	* a_0 = c_{n-1} = 0
	* i = 0..(n-1)
	* n equations
	*/
	std::vector<double> xi = std::vector<double>(n + 1);
	std::vector<double> theta = std::vector<double>(n + 1);
	xi[0] = 0;
	theta[0] = 0;
	for (int i = 0; i < n; i++) {
		double denominator = b[i] - a[i] * xi[i];
		xi[i + 1] = c[i] / denominator;
		theta[i + 1] = (a[i] * theta[i] - d[i]) / denominator;
	}
	x = std::vector<double>(n);
	x[n - 1] = theta[n];
	for (int i = n - 2; i >= 0; i--) {
		x[i] = xi[i + 1] * x[i + 1] + theta[i + 1];
	}
}

// int solve_SOR(int N, const BandMatrix& A, const std::vector<double>& b, std::vector<double>& x, double w, double eps = 1e-6) {
// 	std::vector<double> x_old;
// 	double error;
// 	int steps = 0;
// 	do {
// 		x_old = x;
// 		error = 0;
// 		for (int i = 0; i < N; ++i) {
// 			x[i] = b[i];
// 			for (int j = 0; j < i; ++j) {
// 				x[i] -= A(i, j) * x[j];
// 			}
// 			for (int j = i + 1; j < N; ++j) {
// 				x[i] -= A(i, j) * x_old[j];
// 			}
// 			x[i] = x[i] * w / A(i, i) + (1 - w) * x_old[i];

// 			error = std::max(error, abs(x[i] - x_old[i]));
// 		}
// 		++steps;
// 	} while (error > eps);

// 	return steps;
// }


// double dot_product(std::vector<double> vect_A, std::vector<double> vect_B)
// {
// 	double product = 0;

// 	for (int i = 0; i < vect_A.size(); i++)
// 		product = product + vect_A[i] * vect_B[i];
// 	return product;
// }

// int conj_grad_matrix_solver(int N, const BandMatrix& A, const std::vector<double>& b, std::vector<double>& x, double eps = 1e-6) {
// 	// init r0
// 	std::vector<double> r_old(N);
// 	for (int i = 0; i < N; ++i) {
// 		r_old[i] = b[i];
// 		for (int j = 0; j < N; ++j) {
// 			r_old[i] -= A(i, j) * x[j];
// 		}
// 	}
// 	// init z0
// 	std::vector<double> z = r_old;
	
// 	std::vector<double> r(N);
// 	double alpha;
// 	double beta;
// 	double error = 0;
// 	int steps = 0;
// 	do {
// 		error = 0;

// 		std::vector<double> Az(N);

// 		for (int j = 0; j < N; ++j) {
// 			Az[j] = 0;
// 			for (int k = 0; k < N; ++k) {
// 				Az[j] += A(j, k) * z[k];
// 			}
// 		}
// 		alpha = dot_product(r_old, r_old) / dot_product(Az, z);

// 		for (int i = 0; i < N; ++i) {
// 			double x_old = x[i];
// 			x[i] = x[i] + alpha * z[i];
// 			r[i] = r_old[i] - alpha * Az[i];
// 			error = std::max(error, abs(x[i] - x_old));
// 		}
// 		beta = dot_product(r, r) / dot_product(r_old, r_old);


// 		for (int i = 0; i < N; ++i) {
// 			z[i] = r[i] + beta * z[i];
// 		}

// 		//error = *std::max_element(r.begin(), r.end()) / *std::max_element(b.begin(), b.end());
// 		//error = eps + 1;
// 		r_old = r;
// 		++steps;
// 	} while (error > eps);
// 	return steps;
// }


auto solve_gauss(std::span<double> A, std::span<double> b) {
	const size_t n = b.size();
	std::vector<double> x (n);
	for (int i = 0; i < n; ++i) {
		double pivot = A[i + i * n];
		if (abs(pivot) < EPSILON) {
			// Если пилотный элемент равен нулю
			double max = 0.;
			int max_index = i;
			for (int j = i + 1; j < n; ++j) {
				if (abs(A[i + j * n]) > max) {
					max_index = j;
					max = A[j + j * n];
				}
			}

			if (max < EPSILON) {
				throw std::invalid_argument("det(A) = 0");
			}
			for (int k = i; k < n; ++k) {
				std::swap(A[k + max_index * n], A[k + i * n]);
			}
			std::swap(b[i], b[max_index]);
			i--;
			continue;
		}
		for (int j = i; j < n; ++j) {
			A[j + i * n] /= pivot;
		}
		b[i] /= pivot;

		for (int j = i + 1; j < n; ++j) {
			pivot = A[i + j * n];
			for (int k = 0; k < n; ++k) {
				A[k + j * n] -= pivot * A[k + n * i];
			}
			b[j] -= pivot * b[i];
		}
	}

	//backward
	for (int i = n - 1; i >= 0; --i) {
		x[i] = b[i];
		for (int j = i + 1; j < n; ++j) {
			x[i] -= A[j + i * n] * x[j];
		}
	}

	return x;
}