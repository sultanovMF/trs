#include <iostream>
#include <span>
#include <vector>
#include <functional>
#include <fstream>
#include <tuple>
#include <format> // works in MCVS 17.1.1 and above, for linux wait GCC 13 
#include <array>
#include "ode.hpp"
#include "matrix_solvers.hpp"
#include <span>



// http://www.ees.nmt.edu/outside/courses/hyd510/PDFs/Lecture%20notes/Lectures%20Part%202.6%20FDMs.pdf

void write_results(std::fstream &output, std::vector<std::pair<double, std::vector<double>>> res, int exnum, int n) {
	std::string path = std::format("../../../../results/lab_1/ex_{}_{}.txt", exnum, n);
	output.open(path, std::ios::out);
	if (!output)
	{
		std::cout << "File not created!";
	}
	else
	{
		for (auto [t, y] : res)
		{
			output << t << "\t";
			for (int i = 0; i < y.size(); ++i) {
				output << y[i];
				if (i != y.size() - 1) {
					output << "\t";
				}
			}
			output << "\n";
		}
		output.close();
	}
}

void write_results(std::fstream &output, std::vector<std::pair<double, double>> res, int exnum, int n) {
	std::string path = std::format("../../../../results/lab_1/ex_{}_{}.txt", exnum, n);
	output.open(path, std::ios::out);
	if (!output)
	{
		std::cout << "File not created!";
	}
	else
	{
		for (auto [t, y] : res) {
			output << t << "\t" << y << "\n";
		}
		output.close();
	}
}



int main()
{
	std::fstream output;
	{
		// Part I
		CauchyProblem problem;
		problem.foo = [] (double t, std::vector<double> y) {
			auto x_t = y[1]; // dx/dt = v
			auto v_t = -2 * y[0] * sinh(y[0] * y[0]) + 2 * y[0] * y[0] * y[0] * cosh(y[0] * y[0]);
			return std::vector<double>({ x_t, v_t });
		};

		problem.t0 = 0;
		problem.tf = 1;
		problem.y = {0, 1};

		double t_len = (problem.tf - problem.t0);

		for (int n : {10, 50, 100, 150, 200, 300, 450, 500, 550, 600, 650, 750, 800, 850, 900, 950, 1000})
		{
			double h = t_len / n;
			write_results(output, euler_method(problem, h, n), 1, n);
			write_results(output, adams_method(problem, h, n), 2, n);
			write_results(output, rk4_method  (problem, h, n), 3, n);
		}
	}
	{
		// Part II
		BVP1d problem;
		problem.p = [] (double x) {
			return -2.;
		};

		problem.q = [] (double x) {
			return 0.;
		};
		
		problem.f = [] (double x) {
			return std::exp(x) * (x*x + x - 3);
		};

		problem.A = 2;
		problem.alpha[0] = 0;
		problem.alpha[1] = 1;

		problem.B = 2;
		problem.beta[0] = 1;
		problem.beta[1] = 1;

		problem.x0 = 0;
		problem.xl = 1;

		for (auto n : {10, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000}) {
			//write_results(output, fd1d(problem, n), 4, n);

			/*auto result = shooting(problem, n);
			double error = 0.;
			double h = 1 / n;

			for (int i = 0; i < n; ++i) {
				double x = result[i].first;
				error = std::max(error, std::abs(result[i].second[0] - (exp((double)(2 * x)) - exp((double)x) * (double)x + exp((double)x) - exp((double)x) * (double)(x * x) - 0.3e1 * exp(0.2e1) + 0.5e1 * exp(0.1e1) + 0.2e1)));
			}
			std::cout << "n = " << n << " error " << error << std::endl;*/

			//write_results(output, shooting(problem, n), 5, n);
		}
	}
	return 0;
}
