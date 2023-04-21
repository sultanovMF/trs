#pragma once 

#include <vector>
#include <functional>
#include <random>
#include "matrix_solvers.hpp"
std::random_device rseed;
std::mt19937 rng(rseed());
std::uniform_int_distribution<int> dist(-10000, 10000);

struct CauchyProblem
{
	std::function<std::vector<double>(double, std::vector<double>)> foo;
	std::vector<double> y;
	double t0;
	double tf;
};

auto euler_method(const CauchyProblem& problem, double h, int n) {
	auto t = problem.t0;
	auto y = problem.y;

	std::vector<std::pair<double, std::vector<double>>> result;
	result.push_back(std::pair(t, y)); // handler
	for (int i = 1; i < n; i++)
	{
		auto y_prev = problem.foo(t, y);
		for (int j = 0; j < y.size(); ++j)
		{
			y[j] += h * y_prev[j];
		}
		t += h;

		result.push_back(std::pair(t, y)); // handler
	}
	return result;
}

auto adams_method(const CauchyProblem& problem, double h, int n) {
	auto t = problem.t0;
	auto y = problem.y;

	std::vector<std::pair<double, std::vector<double>>> result;
	result.push_back(std::pair(t, y)); // handler T0 Y0

	// Set Y0
	auto y_prev_prev = y;
	auto foo_prev_prev = problem.foo(t, y_prev_prev); // f(T0, Y0)

	// Calculate Y1 using Euler method
	auto y_prev = y_prev_prev;
	for (int j = 0; j < y.size(); ++j)
	{
		y_prev[j] += h * foo_prev_prev[j];
	}

	t += h;
	auto foo_prev = problem.foo(t, y_prev); // f(T1, Y1)

	result.push_back(std::pair(t, y_prev)); // handler T1 Y1

	for (int i = 2; i < n; i++)
	{
		// Calculate Yn+1
		y = y_prev;
		for (int j = 0; j < y.size(); ++j)
		{
			y[j] += h / 2. * (3 * foo_prev[j] - foo_prev_prev[j]);
		}

		t += h;

		result.push_back(std::pair(t, y)); // handler

		y_prev_prev = y_prev;
		y_prev = y;

		foo_prev_prev = foo_prev; // f(Tn, Yn)
		foo_prev = problem.foo(t, y);	  // f(Tn+1, Yn+1)
	}
	return result;
}

auto rk4_method(const CauchyProblem& problem, double h, int n) {
	auto t = problem.t0;
	auto y = problem.y;
	std::vector<std::pair<double, std::vector<double>>> result;

	result.push_back(std::pair(t, y)); // handler T0 Y0

	auto shift = [](double new_h, std::vector<double> k, std::vector<double> y)
	{
		std::vector<double> temp = y;
		for (int i = 0; i < y.size(); ++i)
		{
			temp[i] += new_h * k[i];
		}
		return temp;
	};

	for (int i = 1; i < n; ++i)
	{
		auto k1 = problem.foo(t, y);
		auto k2 = problem.foo(t + h / 2., shift(h / 2., k1, y));
		auto k3 = problem.foo(t + h / 2., shift(h / 2., k2, y));
		auto k4 = problem.foo(t + h,      shift(h,      k3, y));

		for (int j = 0; j < y.size(); ++j)
		{
			y[j] += h / 6. * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]);
		}

		t += h;

		result.push_back(std::pair(t, y)); // handler
	}

	return result;
}


struct BVP1d {
	// Second Order BVP problem
	// u'' + p(x) u' + q(x) u = f(x)
	// with mixed boundary conditions

	std::function<double(double)> p;
	std::function<double(double)> q;
	std::function<double(double)> f;

	// left boundary condition coeffs.
	double alpha[2];
	double A;

	// right boundary condition coeffs.
	double beta[2];
	double B;

	double x0;
	double xl;
};

auto fd1d(const BVP1d& problem, int n) {
	double h = (problem.xl - problem.x0) / (n - 1.);

	std::vector<std::pair<double, double>> result(n);

	std::vector<double> A(n);
	std::vector<double> B(n);
	std::vector<double> C(n);
	std::vector<double> D(n);

	// left points approximation x0 
	A[0] = 0; // за границей
	B[0] = h * problem.alpha[0] - problem.alpha[1];
	C[0] = problem.alpha[1];
	D[0] = problem.A * h;

	for (int i = 1; i < n - 1; ++i) {
		// inner points approximation
		double xi = problem.x0 + h * i;

		A[i] = 2 - h * problem.p(xi);     // (2 - p(Xi) * h) Ui-1
		B[i] = 2 * h * h * problem.q(xi) - 4; // (2 q(Xi) h^2 - 4) Ui
		C[i] = 2 + h * problem.p(xi);     // (2 + h p(Xi) ) Ui+1
		D[i] = 2 * h * h * problem.f(xi);     //  2 h^2 f(Xi)
	}

	// right points approximation xl
	A[n - 1] = problem.beta[0] * h - problem.beta[1]; // (2 - p(Xi) * h) Ui-1
	B[n - 1] = problem.beta[1]; // (2 q(Xi) h^2 - 4) Ui
	C[n - 1] = 0; // за границей
	D[n - 1] = problem.B * h;

	auto y = solve_tridiagonal(A, B, C, D);
	double error = 0;
	for (int i = 0; i < n; ++i) {
		double x = problem.x0 + h * i;
		result[i] = std::pair(x, y[i]);
		error = std::max(error, std::abs(y[i] - (exp((double)(2 * x)) - exp((double)x) * (double)x + exp((double)x) - exp((double)x) * (double)(x * x) - 0.3e1 * exp(0.2e1) + 0.5e1 * exp(0.1e1) + 0.2e1)));
	}
	std::cout << "ex 1 n = " << n << " error " << error << std::endl;

	return result;
}

auto shooting(const BVP1d& bvp_problem, int n) {
	// for system of 2 ode
	double h = (bvp_problem.xl - bvp_problem.x0) / (n - 1.);
	std::vector<std::pair<double, double>> result(n);

	CauchyProblem ivp_problem;
	ivp_problem.foo = [&](double t, std::vector<double> y) {
		auto u_t = y[1];
		auto v_t = bvp_problem.f(t) - bvp_problem.p(t) * y[1] - bvp_problem.q(t) * y[0];
		return std::vector<double>({ u_t, v_t });
	};
	ivp_problem.t0 = bvp_problem.x0;
	ivp_problem.tf = bvp_problem.xl;

	double mu;
	double nu = dist(rng);

	do {
		// Пристреливаемся
		mu = nu;
		nu = dist(rng);

		ivp_problem.y = { mu, (bvp_problem.A - bvp_problem.alpha[0] * mu) / bvp_problem.alpha[1] };
		auto mu_sol = adams_method(ivp_problem, h, n);
		ivp_problem.y = { nu, (bvp_problem.A - bvp_problem.alpha[0] * nu) / bvp_problem.alpha[1] };
		auto nu_sol = adams_method(ivp_problem, h, n);

		double mu_psi = mu_sol[n - 1].second[0] * bvp_problem.beta[0] + mu_sol[n - 1].second[1] * bvp_problem.beta[1] - bvp_problem.B;
		double nu_psi = nu_sol[n - 1].second[0] * bvp_problem.beta[0] + nu_sol[n - 1].second[1] * bvp_problem.beta[1] - bvp_problem.B;
		if (mu_psi * nu_psi < 0) {
			if (mu > nu) std::swap(mu, nu);
			break;
		}
	} while (true);

	// Метод Дихотомии

	while (nu - mu > 0.00000001) {
		double mid = (mu + nu) / 2;

		ivp_problem.y = { nu, (bvp_problem.A - bvp_problem.alpha[0] * nu) / bvp_problem.alpha[1] };
		auto nu_sol = adams_method(ivp_problem, h, n);
		ivp_problem.y = { mid, (bvp_problem.A - bvp_problem.alpha[0] * mid) / bvp_problem.alpha[1] };
		auto mid_sol = adams_method(ivp_problem, h, n);

		double mid_psi = mid_sol[n - 1].second[0] * bvp_problem.beta[0] + mid_sol[n - 1].second[1] * bvp_problem.beta[1] - bvp_problem.B;
		double nu_psi = nu_sol[n - 1].second[0] * bvp_problem.beta[0] + nu_sol[n - 1].second[1] * bvp_problem.beta[1] - bvp_problem.B;

		if (mid_psi * nu_psi < 0) {
			mu = mid;
		}
		else {
			nu = mid;
		}
	}

	return adams_method(ivp_problem, h, n);
}