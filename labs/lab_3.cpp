#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <vector>
#include <span>
using namespace std;


const double delta = 1e-8;

std::vector<double> solve_tridiagonal(
	const std::span<double> A,
	const std::span<double> B,
	const std::span<double> C,
	const std::span<double> D)
{
	if (!(A.size() == B.size() && B.size() == C.size() && C.size() == D.size())) {
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

struct EllipticPDEMatrixGen {
	double alpha;
	double beta;
	double gamma;
																																																													
	int N;
	int M;
	EllipticPDEMatrixGen(double alpha, double beta, double gamma, int N, int M) : alpha(alpha), beta(beta), gamma(gamma), N(N), M(M) {}

	const double& operator()(int i, int j) const {
		if (i == j) {
			return alpha;
		}
		else if (i - j == -1 && (i+1) % M != 0) {
			return beta;
		}
		else if (i - j == 1 && (i % M != 0)) {
			return beta;
		}
		else if (abs(i - j) == M) {
			return gamma;
		}
		else {
			return 0.;
		}
	}
};

void lu_decompostion(const unsigned int n, EllipticPDEMatrixGen A, vector<double>& L, vector<double>& U) {
	for (int i = 0; i < n; ++i) {
		L[i * n] = A(i, 0);
		U[i] = A(i, 0) / L[0];
	}

	for (int i = 1; i < n; i++)
	{
		for (int j = i; j < n; j++)
		{
			U[i * n + j] = A(i,j);

			for (int k = 0; k < i; k++)
			{
				U[i * n + j] -= L[i * n + k] * U[k * n + j];
			}

			L[j * n + i] = A(i, j);

			for (int k = 0; k < i; k++)
			{
				L[j * n + i] -= L[j * n + k] * U[k * n + i];
			}

			L[j * n + i] /= U[i * n + i];
		}
	}
	//for (int i = 0; i < n; ++i) {
	//	L[i * n] = A[i * n];
	//	U[i] = A[i] / L[0];
	//}

	//for (int i = 1; i < n; i++)
	//{
	//	for (int j = i; j < n; j++)
	//	{
	//		U[i * n + j] = A[i * n + j];

	//		for (int k = 0; k < i; k++)
	//		{
	//			U[i * n + j] -= L[i * n + k] * U[k * n + j];
	//		}

	//		L[j * n + i] = A[j * n + i];

	//		for (int k = 0; k < i; k++)
	//		{
	//			L[j * n + i] -= L[j * n + k] * U[k * n + i];
	//		}

	//		L[j * n + i] /= U[i * n + i];
	//	}
	//}

}


void backward_up(const unsigned int n, vector<double> U, vector<double> b, vector<double>& x) {
	for (int i = n - 1; i >= 0; --i) {
		x[i] = b[i];
		for (int j = i + 1; j < n; ++j) {
			x[i] -= U[j + i * n] * x[j];
		}
		x[i] /= U[i + i * n];
	}
}

void backward_low(const unsigned int n, vector<double>  L, vector<double>  b, vector<double> & y) {
	for (int i = 0; i < n; ++i) {
		y[i] = b[i];
		for (int j = 0; j < i; ++j) {
			y[i] -= L[i * n + j] * y[j];
		}
		y[i] /= L[i + i * n];
	}
}


void solve_lu(const unsigned int n, vector<double>  L, vector<double>  U, vector<double>  b, vector<double>& x) {
	vector<double> y(n);
	backward_low(n, L, b, y);
	backward_up(n, U, y, x);
}

vector<double> Ax(double alpha, double beta, double gamma, vector<double> U_prev, double N, double M) {
	vector<double> U_mult(N * M, 0.);

	// умножаем левые гаммы
	for (int j = 1; j < N; ++j) {
		for (int i = 0; i < M; ++i) {
			U_mult[j * M + i] += - gamma * U_prev[(j - 1) * M + i] /alpha;
		}
	}
	// умножаем правые гаммы
	for (int j = 0; j < N - 1; ++j) {
		for (int i = 0; i < M; ++i) {
			U_mult[j * M + i] += - gamma * U_prev[(j + 1) * M + i] / alpha;
		}
	}

	// умножаем левые беты
	for (int j = 0; j < N; ++j) {
		for (int i = 1; i < M; ++i) {
			U_mult[j * M + i] += - beta * U_prev[j * M + (i - 1)] / alpha;
		}
	}

	// умножаем правые беты
	for (int j = 0; j < N; ++j) {
		for (int i = 0; i < M - 1; ++i) {
			U_mult[j * M + i] += - beta * U_prev[j * M + (i + 1)] / alpha;
		}
	}

	return U_mult;
}

int main() {
	// Метод простых итераций с правильным параметром tau https://zaochnik.com/spravochnik/matematika/issledovanie-slau/iteratsionnye-metody-reshenija-slau/
	vector<int> N_arr = { 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100 };

	{
		// Part I
		const double lx = 1.;
		const double ly = 2.;

		auto f = [](double x, double y) {
			return -0.2e1 * 0.3141592654e1 * cos(0.3141592654e1 * x * x) * sin(0.3141592654e1 * y / 0.2e1) + 0.4e1 * 0.3141592654e1 * 0.3141592654e1 * x * x * sin(0.3141592654e1 * x * x) * sin(0.3141592654e1 * y / 0.2e1) + sin(0.3141592654e1 * x * x) * 0.3141592654e1 * 0.3141592654e1 * sin(0.3141592654e1 * y / 0.2e1) / 0.4e1;;
		};

		auto phi1 = [](double x = 0, double y = 0) {
			return 0;
		};

		auto phi2 = [](double x = 0, double y = 0) {
			return 0;
		};

		auto phi3 = [](double x = 0, double y = 0) {
			return 0;
		};

		auto phi4 = [](double x = 0, double y = 0) {
			return 0;
		};


		auto exact_sol = [lx, ly](double x, double y) {
			return sin(M_PI * x * x / lx / lx) * sin(M_PI * y / ly);
		};


		auto getF = [phi1, phi2, phi3, phi4, f](int N, int M, double hx, double hy) {
			vector<double> F(N * M);
			{
				int i, j;
				// F для левой границы
				i = 0;
				j = 0;
				F[j * M + i] = -f(i * hx, j * hy) - phi1() / hx / hx - phi3() / hy / hy;
				for (j = 1; j < N - 1; ++j) {
					F[j * M + i] = -f(i * hx, j * hy) - phi1() / hx / hx;
				}
				F[j * M + i] = -f(i * hx, j * hy) - phi1() / hx / hx - phi4() / hy / hy;

				// F для верхней границы
				for (i = 1; i < M - 1; ++i) {
					F[j * M + i] = -f(i * hx, j * hy) - phi4() / hy / hy;
				}

				F[j * M + i] = -f(i * hx, j * hy) - phi2() / hx / hx - phi4() / hy / hy;

				// F для правой границы
				for (j = N - 2; j > 0; --j) {
					F[j * M + i] = -f(i * hx, j * hy) - phi2() / hx / hx;
				}
				F[j * M + i] = -f(i * hx, j * hy) - phi2() / hx / hx - phi3() / hy / hy;

				// F  для нижней границы
				for (i = M - 2; i > 0; --i) {
					F[j * M + i] = -f(i * hx, j * hy) - phi3() / hy / hy;
				}

				// внутренняя область
				for (j = 1; j < N - 1; ++j) {
					for (i = 1; i < M - 1; ++i) {
						F[j * M + i] = -f(i * hx, j * hy);
					}
				}
			}
			return F;
		};

		// Метод простых итераций
		for (auto N : N_arr) { // N число по y
			double hy = ly / N;
			double hx = hy * lx / ly;

			int M = long long(lx / hx); // M число по x

			vector<double> U(N * M, 0.); // Начальное приближение нулевой вектор
			vector<double> U_prev(N * M);

			double alpha =- 2. * (1 / hx / hx + 1 / hy / hy);
			double beta = 1 / hx / hx;
			double gamma = 1 / hy / hy;

			vector<double> F = getF(N, M, hx, hy);

			// Метод простых итераций u^(k+1) = b + A * x^(k) https://ru.wikipedia.org/wiki/%D0%9C%D0%B5%D1%82%D0%BE%D0%B4_%D0%B8%D1%82%D0%B5%D1%80%D0%B0%D1%86%D0%B8%D0%B8
			
			vector<double> b(N * M);

			for (int k = 0; k < N * M; ++k) {
				b[k] = F[k] / alpha;
			}

			auto Ax = [alpha, beta, gamma, N, M](vector<double> U_prev) {
				vector<double> U_mult(N * M, 0.);

				// умножаем левые гаммы
				for (int j = 1; j < N; ++j) {
					for (int i = 0; i < M; ++i) {
						U_mult[j * M + i] += -gamma * U_prev[(j - 1) * M + i] / alpha;
					}
				}
				// умножаем правые гаммы
				for (int j = 0; j < N - 1; ++j) {
					for (int i = 0; i < M; ++i) {
						U_mult[j * M + i] += -gamma * U_prev[(j + 1) * M + i] / alpha;
					}
				}

				// умножаем левые беты
				for (int j = 0; j < N; ++j) {
					for (int i = 1; i < M; ++i) {
						U_mult[j * M + i] += -beta * U_prev[j * M + (i - 1)] / alpha;
					}
				}

				// умножаем правые беты
				for (int j = 0; j < N; ++j) {
					for (int i = 0; i < M - 1; ++i) {
						U_mult[j * M + i] += -beta * U_prev[j * M + (i + 1)] / alpha;
					}
				}

				return U_mult;
			};


			double error = 0.;
			long long iter_count = 0;
			do {
				iter_count++;
				copy(U.begin(), U.end(), U_prev.begin()); // копируем предыдущее значение вектора

				auto ax = Ax(U_prev);

				for (int j = 0; j < N; ++j) {
					for (int i = 0; i < M; ++i) {
						U[j * M + i] = b[j * M + i] + ax[j * M + i];
					}
				}

				error = 0.;
				for (int k = 0; k < N * M; ++k) {
					error = max(abs(U_prev[k] - U[k]), error);
				}

			} while (error > delta);

			error = 0.;

			for (int j = 0; j < N; ++j) {
				for (int i = 0; i < M; ++i) {
					double x = i * hx;
					double y = j * hy;
					error = max(error, abs(U[j * M + i] - sin(M_PI * x * x / lx / lx) * sin(M_PI * y / ly)));
				}
			}

			cout << "Simple iterations N = " << N << " M  = " << M << " " << " error: " << error  << " iter count: " << iter_count << endl;

		}

		// Метод SOR
		//for (auto N : N_arr) {
		//	double hy = ly / N;
		//	double hx = hy * lx / ly;

		//	int M = long long(lx / hx); // M число по x

		//	vector<double> U(N * M, 0.); // Начальное приближение нулевой вектор
		//	vector<double> U_prev(N * M);

		//	double alpha = - 2. * (1 / hx / hx + 1 / hy / hy);
		//	double beta  = 1 / hx / hx;
		//	double gamma = 1 / hy / hy;

		//	double omega = 1.;

		//	vector<double> F = getF(N, M, hx, hy);

		//		double error = 0.;
		//		long long iter_count = 0;
		//		do {
		//			iter_count++;
		//			copy(U.begin(), U.end(), U_prev.begin()); // копируем предыдущее значение вектора
		//			int i, j;

		//			// Первый блок

		//			i = 0;
		//			j = 0;

		//			U[j * M + i] = (1 - omega) * U_prev[j * M + i] 
		//				+ omega / alpha * (F[j * M + i] - beta * U_prev[j * M + i + 1] - gamma * U_prev[(j+1) * M + i]);

		//			for (i = 1; i < M - 1; ++i) {
		//				U[j * M + i] = (1 - omega) * U_prev[j * M + i] + 
		//					omega / alpha * (F[j * M + i] - beta * U[j * M + (i - 1)] - beta * U_prev[j * M + i + 1] - gamma * U_prev[(j + 1) * M + i]);
		//			}

		//			U[j * M + i] = (1 - omega) * U_prev[j * M + i] +
		//				omega / alpha * (F[j * M + i] - beta * U[j * M + (i - 1)] - gamma * U_prev[(j + 1) * M + i]);

//			// Внутренние блоки
//			for (j = 1; j < N - 1; ++j) {
//				i = 0;
//				U[j * M + i] = (1 - omega) * U_prev[j * M + i] + 
//					omega / alpha * (F[j * M + i] - gamma * U[(j-1) * M + i] - beta * U_prev[j * M + i + 1] - gamma * U_prev[(j + 1) * M + i]);

//				for (i = 1; i < M - 1; ++i) {
//					U[j * M + i] = (1 - omega) * U_prev[j * M + i] +
//						omega / alpha * (F[j * M + i] - gamma * U[(j - 1) * M + i] - beta * U[j * M + (i - 1)] - beta * U_prev[j * M + i + 1] - gamma * U_prev[(j + 1) * M + i]);
//				}

//				U[j * M + i] = (1 - omega) * U_prev[j * M + i] +
//					omega / alpha * (F[j * M + i] - gamma * U[(j - 1) * M + i] - beta * U[j * M + (i - 1)] - gamma * U_prev[(j + 1) * M + i]);
//			}


//			// Последний блок
//			i = 0;
//			U[j * M + i] = (1 - omega) * U_prev[j * M + i] +
//				omega / alpha * (F[j * M + i] - gamma * U[(j - 1) * M + i] - beta * U_prev[j * M + i + 1]);

//			for (i = 1; i < M - 1; ++i) {
//				U[j * M + i] = (1 - omega) * U_prev[j * M + i] +
//					omega / alpha * (F[j * M + i] - gamma * U[(j - 1) * M + i] - beta * U[j * M + (i - 1)] - beta * U_prev[j * M + i + 1]);
//			}

//			U[j * M + i] = (1 - omega) * U_prev[j * M + i] +
//				omega / alpha * (F[j * M + i] - gamma * U[(j - 1) * M + i] - beta * U[j * M + (i - 1)]);


//			error = 0.;
//			for (int k = 0; k < N * M; ++k) {
//				error = max(abs(U_prev[k] - U[k]), error);
//			}

//		} while (error > delta);

//		error = 0.;

//		for (int j = 0; j < N; ++j) {
//			for (int i = 0; i < M; ++i) {
//				double x = i * hx;
//				double y = j * hy;
//				error = max(error, abs(U[j * M + i] - sin(M_PI * x * x / lx / lx) * sin(M_PI * y / ly)));
//			}
//		}

//		cout << "SOR method N = " << N << " M  = " << M << " " << " error: " << error  << " iter count: " << iter_count << endl;
//}
//int N = 4;
//int M = 5;
//EllipticPDEMatrixGen pm(1, 2, 3, N, M);
//for (int i = 0; i < N * M; ++i) {
//	for (int j = 0; j < M * N; ++j) {
//		cout << pm(i, j) << " ";
//	}
//	cout << endl;
//}

//// Прямой метод:
//{
//	for (auto N : N_arr) {
//		double hy = ly / N;
//		double hx = hy * lx / ly;
//
//		int M = long long(lx / hx); // M число по x
//
//		vector<double> U(N * M, 0.); // Начальное приближение нулевой вектор
//
//		double alpha = -2. * (1 / hx / hx + 1 / hy / hy);
//		double beta = 1 / hx / hx;
//		double gamma = 1 / hy / hy;
//
//		vector<double> F = getF(N, M, hx, hy);
//
//		EllipticPDEMatrixGen pm(alpha, beta, gamma, N, M);
//
//		vector<double> Lower(M * N * M * N);
//		vector<double> Upper(M * N * M * N);
//
//		lu_decompostion(N * M, pm, Lower, Upper);
//
//		solve_lu(N * M, Lower, Upper, F, U);
//
//
//		double error = 0.;
//
//		for (int j = 0; j < N; ++j) {
//			for (int i = 0; i < M; ++i) {
//				double x = i * hx;
//				double y = j * hy;
//				error = max(error, abs(U[j * M + i] - sin(M_PI * x * x / lx / lx) * sin(M_PI * y / ly)));
//			}
//		}
//
//		cout << "Gauss N = " << N << " M  = " << M << " " << " error: " << error << endl;
//	}
//		}

	}
	{
		// Part 2
	}
}