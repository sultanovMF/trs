#include <iostream>
#include <span>
#include <vector>
#include <functional>
#include <fstream>
#include <cmath>
#include <fstream>
#include <tuple>
#include <format>
#include "matrix_solvers.hpp"
// http://math.phys.msu.ru/data/785/uravneniya_teploprovodnosti.pdf

void write_results(std::fstream& output, std::vector<std::tuple<double, double, double>> res, int exnum, int n) {
    std::string path = std::format("../../../../results/lab_2/ex_{}_{}.txt", exnum, n);
    output.open(path, std::ios::out);
    if (!output)
    {
        std::cout << "File not created!";
    }
    else
    {
        for (auto [t, x, u] : res) {
            output << t << "\t" << x << "\t" << u << "\n";
        }
        output.close();
    }
}

class Handler {
    std::fstream output_res;
    std::fstream output_error;
    double ex_num;
    int n;
    std::string res_path;
    std::string error_path;
    std::function<double(double, double)> target;
    double max_error = 0;
    bool key_result;
    bool key_error;
public:
    Handler(double ex_num, int n, bool key_result = true, bool key_error = false, std::function<double(double, double)> target = 0)
        : ex_num(ex_num), n(n), target(target), key_error(key_error), key_result(key_result) {
        if (key_result) {
            res_path   = std::format("../../../../results/lab_2/ex_{}_{}.txt", ex_num, n);
            output_res.open(res_path, std::ios::out);

            if (!output_res)
            {
                std::cout << "Cannot open files!!";
            }
        }

        if (key_error) {
            error_path = std::format("../../../../results/lab_2/ex_{}_{}_error.txt", ex_num, n);
            output_error.open(error_path, std::ios::out);

            if (!output_error)
            {
                std::cout << "Cannot open files!!";
            }
        }
    }
    
    void Handle(double t, double x, double u) {
        
        if(key_error) {
            double error = std::max(error, std::abs(u - target(t, x)));
            max_error = std::max(max_error, error);
            if (key_result) {
                output_error << n << "\t" << error << "\n";
            }
        }
        if (key_result) {
            output_res << t << "\t" << x << "\t" << u << "\n";
        }
    }

    ~Handler() {
        output_res.close();
        if (key_error) {
            std::cout << ex_num << " " << n << " max error: " << max_error << std::endl;
        }
        output_error.close();
    }

};

int main() {
    const double tf = 1.;
    std::vector<int> N_arr = { 10, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000 };
    {
        double k = 1;
        auto f = [](double x, double t) {
            return t * x * x - t * t;
        };

        double alpha[2] = { 0., 1. };
        double beta[2] = { 1., 1. };

        auto phi = [](double x) {
            return x;
        };
        auto psi_left = [](double t) {
            return t;
        };

        auto psi_right = [](double t) {
            return 3. / 2. * t * t + 2;
        };
        // ex 1
        for (int N : N_arr) {
            Handler handler(1, N, false, true, [](double t, double x) {return t * t * x * x / 2. + x; });

            double xi = 1. / N;
            N = N + 1;
            double tau = xi;

            int K = 100000;
            // IC
            std::vector<double> u(N);
            std::vector<double> u_prev(N);

            for (int i = 0; i < N; ++i) {
                u[i] = phi(i * xi);
            }

            for (int j = 1; j < K; ++j) {
                std::copy(u.begin(), u.end(), u_prev.begin());

                for (int i = 1; i < N - 1; ++i) {
                    u[i] = tau * k / xi / xi * (u_prev[i + 1] - 2 * u_prev[i] + u_prev[i - 1]) + u_prev[i] + tau * f(xi * i, tau * j);
                }

                u[0] = (psi_left(j * tau) - beta[0] / xi * u[1]) / (alpha[0] - (beta[0] / xi));
                u[N - 1] = (psi_right(j * tau) + beta[1] / xi * u[N - 2]) / (alpha[1] + beta[1] / xi);


                for (int i = 0; i < N; ++i) {
                    handler.Handle(j * tau, i * xi, u[i]);
                }
            }

        }
        //    //    }

        //        //    // ex 2.1
        //        //    for (int N : N_arr) {
        //        //        Handler handler(2.1, N, false, true, [](double t, double x) {return t * t * x * x / 2. + x; });

        //        //        double xi = 1. / N;
        //        //        double tau = xi;

        //        //        int K = tf / tau + 1;
        //        //        std::vector<double> u(N);
        //        //        std::vector<double> u_prev(N);

        //        //        std::vector<double> A(N);
        //        //        std::vector<double> B(N);
        //        //        std::vector<double> C(N);
        //        //        std::vector<double> D(N);

        //        //        // ic
        //        //        for (int i = 0; i < N; ++i) {
        //        //            u[i] = phi(xi * i);
        //        //        }

        //        //        for (int j = 1; j < K; ++j) {
        //        //            // система составлена неправильно
        //        //            std::copy(u.begin(), u.end(), u_prev.begin());

        //        //            A[0] = 0;
        //        //            B[0] = (alpha[0] - beta[0] / xi);
        //        //            C[0] = beta[0] / xi;
        //        //            D[0] = psi_left(j * tau);

        //        //            for (int i = 1; i < N - 1; ++i) {
        //        //                A[i] = -k / xi / xi;
        //        //                B[i] = (1 / tau + 2 * k / xi / xi);
        //        //                C[i] = -k / xi / xi;
        //        //                D[i] = f(i * xi, j * tau) + u_prev[i] / tau;
        //        //            }

        //        //            A[N - 1] = -beta[1] / xi;
        //        //            B[N - 1] = (alpha[1] + beta[1] / xi);
        //        //            C[N - 1] = 0;
        //        //            D[N - 1] = psi_right(j * tau);

        //        //            u = solve_tridiagonal(A, B, C, D);

        //        //            for (int i = 0; i < N; ++i) {
        //        //                handler.Handle(j * tau, i * xi, u[i]);
        //        //            }
        //        //        }

        //        //    }

        //        //    // ex 2.2
        //        //    for (int N : N_arr) {
        //        //        Handler handler(2.2, N, false, true, [](double t, double x) {return t * t * x * x / 2. + x; });

        //        //        double xi = 1. / N;
        //        //        double tau = xi;

        //        //        int K = tf / tau + 1;
        //        //        std::vector<double> u(N);
        //        //        std::vector<double> u_prev(N);

        //        //        std::vector<double> A(N);
        //        //        std::vector<double> B(N);
        //        //        std::vector<double> C(N);
        //        //        std::vector<double> D(N);

        //        //        // ic
        //        //        for (int i = 0; i < N; ++i) {
        //        //            u[i] = phi(xi * i);
        //        //        }

        //        //        for (int j = 1; j < K; ++j) {
        //        //            // система составлена неправильно
        //        //            std::copy(u.begin(), u.end(), u_prev.begin());

        //        //            A[0] = 0;
        //        //            B[0] = (alpha[0] - beta[0] / xi);
        //        //            C[0] = beta[0] / xi;
        //        //            D[0] = psi_left(j * tau);

        //        //            for (int i = 1; i < N - 1; ++i) {
        //        //                A[i] = -k / xi / xi / 2.;
        //        //                B[i] = (1 / tau + k / xi / xi);
        //        //                C[i] = -k / xi / xi / 2.;
        //        //                D[i] = f(i * xi, j * tau) + u_prev[i] / tau
        //        //                    + 0.5 * (u_prev[i + 1] - 2 * u_prev[i] + u_prev[i - 1]) / xi / xi;
        //        //            }

        //        //            A[N - 1] = -beta[1] / xi;
        //        //            B[N - 1] = (alpha[1] + beta[1] / xi);
        //        //            C[N - 1] = 0;
        //        //            D[N - 1] = psi_right(j * tau);

        //        //            u = solve_tridiagonal(A, B, C, D);

        //        //            for (int i = 0; i < N; ++i) {
        //        //                handler.Handle(j * tau, i * xi, u[i]);
        //        //            }
        //        //            //double t = j * tau;
        //        //            //double x = i * xi;
        //        //            //error = std::max(error, std::abs(u[i] - (t * t * x * x / 2. + x)));
        //        //         // к правому краю ошибка нарастает
        //        //        }

        //        //    }
        //        //}

        //        // PART II
        //    {
        //        //{
        //        //    // for linear
        //        //    auto k = [](double u) {
        //        //        return 1;
        //        //    };
        //        //    auto F = [](double u) {
        //        //        return 1;
        //        //    };
        //        //    auto f = [](double x, double t) {
        //        //        return (t * x * x - t * t);
        //        //    };

        //        //    double alpha[2] = { 0., 1. };
        //        //    double beta[2] = { 1., 1. };

        //        //    auto phi = [](double x) {
        //        //        return x;
        //        //    };
        //        //    auto psi_left = [](double t) {
        //        //        return 1.;
        //        //    };

        //        //    auto psi_right = [](double t) {
        //        //        return 3. / 2. * t * t + 2;
        //        //    };

        //        //    for (int N : N_arr) {
        //        //        Handler handler(3.1, N, false, true, [](double t, double x) {return t * t * x * x / 2. + x; });

        //        //        double xi = 1. / N;
        //        //        double tau = xi;

        //        //        int K = tf / tau + 1;
        //        //        std::vector<double> u(N);
        //        //        std::vector<double> u_prev(N);

        //        //        std::vector<double> A(N);
        //        //        std::vector<double> B(N);
        //        //        std::vector<double> C(N);
        //        //        std::vector<double> D(N);

        //        //        // ic
        //        //        for (int i = 0; i < N; ++i) {
        //        //            u[i] = phi(xi * i);
        //        //        }

        //        //        for (int j = 1; j < K; ++j) {
        //        //            // система составлена неправильно
        //        //            std::copy(u.begin(), u.end(), u_prev.begin());

        //        //            A[0] = 0;
        //        //            B[0] = (alpha[0] - beta[0] / xi);
        //        //            C[0] = beta[0] / xi;
        //        //            D[0] = psi_left(j * tau);

        //        //            for (int i = 1; i < N - 1; ++i) {
        //        //                double ke = (k((i + 1) * xi) + k(i * xi)) / 2.;
        //        //                double kw = (k(i * xi) + k((i - 1) * xi)) / 2.;

        //        //                A[i] = -kw / xi / xi;
        //        //                B[i] = (1 / tau + kw / xi / xi + ke / xi / xi);
        //        //                C[i] = -ke / xi / xi;
        //        //                D[i] = F(u_prev[i]) * f(i * xi, (j - 1) * tau) + u_prev[i] / tau;
        //        //            }

        //        //            A[N - 1] = -beta[1] / xi;
        //        //            B[N - 1] = (alpha[1] + beta[1] / xi);
        //        //            C[N - 1] = 0;
        //        //            D[N - 1] = psi_right(j * tau);

        //        //            u = solve_tridiagonal(A, B, C, D);

        //        //            for (int i = 0; i < N; ++i) {
        //        //                handler.Handle(j * tau, i * xi, u[i]);
        //        //            }
        //        //        }

        //        //    }
        //        //}

        //        //{
        //        //    // for nonlinear

        //        //    auto k = [](double u) {
        //        //        return sin(u);
        //        //    };

        //        //    auto F = [](double u) {
        //        //        return u;
        //        //    };
        //        //    auto f = [](double x, double t) {
        //        //        return (t * x * x - t * t);
        //        //    };

        //        //    double alpha[2] = { 0., 1. };
        //        //    double beta[2] = { 1., 1. };

        //        //    auto phi = [](double x) {
        //        //        return x;
        //        //    };
        //        //    auto psi_left = [](double t) {
        //        //        return 1.;
        //        //    };

        //        //    auto psi_right = [](double t) {
        //        //        return 3. / 2. * t * t + 2;
        //        //    };

        //        //    for (int N : N_arr) {
        //        //        Handler handler(3.2, N, true);

        //        //        double xi = 1. / N;
        //        //        double tau = xi;

        //        //        int K = tf / tau + 1;
        //        //        std::vector<double> u(N);
        //        //        std::vector<double> u_prev(N);

        //        //        std::vector<double> A(N);
        //        //        std::vector<double> B(N);
        //        //        std::vector<double> C(N);
        //        //        std::vector<double> D(N);

        //        //        // ic
        //        //        for (int i = 0; i < N; ++i) {
        //        //            u[i] = phi(xi * i);
        //        //        }

        //        //        for (int j = 1; j < K; ++j) {
        //        //            // система составлена неправильно
        //        //            std::copy(u.begin(), u.end(), u_prev.begin());

        //        //            A[0] = 0;
        //        //            B[0] = (alpha[0] - beta[0] / xi);
        //        //            C[0] = beta[0] / xi;
        //        //            D[0] = psi_left(j * tau);

        //        //            for (int i = 1; i < N - 1; ++i) {
        //        //                double ke = (k((i + 1) * xi) + k(i * xi)) / 2.;
        //        //                double kw = (k(i * xi) + k((i - 1) * xi)) / 2.;

        //        //                A[i] = -kw / xi / xi;
        //        //                B[i] = (1 / tau + kw / xi / xi + ke / xi / xi);
        //        //                C[i] = -ke / xi / xi;
        //        //                D[i] = F(u_prev[i]) * f(i * xi, (j - 1) * tau) + u_prev[i] / tau;
        //        //            }

        //        //            A[N - 1] = -beta[1] / xi;
        //        //            B[N - 1] = (alpha[1] + beta[1] / xi);
        //        //            C[N - 1] = 0;
        //        //            D[N - 1] = psi_right(j * tau);

        //        //            u = solve_tridiagonal(A, B, C, D);

        //        //            for (int i = 0; i < N; ++i) {
        //        //                handler.Handle(j * tau, i * xi, u[i]);
        //        //            }
        //        //        }

        //        //    }
        //        //}


        //        //{
        //        //    // for nonlinear with optimization

        //        //    auto k = [](double u) {
        //        //        return sin(u);
        //        //    };

        //        //    auto F = [](double u) {
        //        //        return u;
        //        //    };
        //        //    auto f = [](double x, double t) {
        //        //        return (t * x * x - t * t);
        //        //    };

        //        //    double alpha[2] = { 0., 1. };
        //        //    double beta[2] = { 1., 1. };

        //        //    auto phi = [](double x) {
        //        //        return x;
        //        //    };
        //        //    auto psi_left = [](double t) {
        //        //        return 1.;
        //        //    };

        //        //    auto psi_right = [](double t) {
        //        //        return 3. / 2. * t * t + 2;
        //        //    };

        //        //    for (int N : N_arr) {
        //        //        Handler handler(4, N, true);

        //        //        double xi = 1. / N;
        //        //        double tau = xi;

        //        //        int K = tf / tau + 1;
        //        //        std::vector<double> u(N);
        //        //        std::vector<double> uu(N);
        //        //        std::vector<double> u_prev(N);

        //        //        std::vector<double> A(N);
        //        //        std::vector<double> B(N);
        //        //        std::vector<double> C(N);
        //        //        std::vector<double> D(N);

        //        //        // ic
        //        //        for (int i = 0; i < N; ++i) {
        //        //            u[i] = phi(xi * i);
        //        //        }
        //        //        double delta;

        //        //        double ke;
        //        //        double kw;

        //        //        for (int j = 1; j < K; ++j) {
        //        //                // система составлена неправильно
        //        //                std::copy(u.begin(), u.end(), u_prev.begin());

        //        //                A[0] = 0;
        //        //                B[0] = (alpha[0] - beta[0] / xi);
        //        //                C[0] = beta[0] / xi;
        //        //                D[0] = psi_left(j * tau);

        //        //                A[N - 1] = -beta[1] / xi;
        //        //                B[N - 1] = (alpha[1] + beta[1] / xi);
        //        //                C[N - 1] = 0;
        //        //                D[N - 1] = psi_right(j * tau);
        //        //                std::copy(u_prev.begin(), u_prev.end(), uu.begin());
        //        //                do {
        //        //                    delta = 0;
        //        //                    for (int i = 1; i < N - 1; ++i) {
        //        //                        ke = (k(uu[i + 1]) + k(uu[i])) / 2.;
        //        //                        kw = (k(uu[i]) + k(uu[i - 1])) / 2.;

        //        //                        A[i] = -kw / xi / xi;
        //        //                        B[i] = (1 / tau + kw / xi / xi + ke / xi / xi);
        //        //                        C[i] = -ke / xi / xi;
        //        //                        D[i] = F(u_prev[i]) * f(i * xi, (j - 1) * tau) + u_prev[i] / tau;
        //        //                    }

        //        //                    u = solve_tridiagonal(A, B, C, D);

        //        //                    for (int i = 0; i < N; ++i) {
        //        //                        delta = std::max(delta, std::abs(u[i] - uu[i]));
        //        //                    }

        //        //                    std::copy(u.begin(), u.end(), uu.begin());

        //        //                } while (delta > 0.00001);


        //        //                for (int i = 0; i < N; ++i) {
        //        //                    handler.Handle(j * tau, i * xi, u[i]);
        //        //                }
        //        //            }
        //        //    }
        //        //}
        //    }


        //}
    }
    return 0;
}