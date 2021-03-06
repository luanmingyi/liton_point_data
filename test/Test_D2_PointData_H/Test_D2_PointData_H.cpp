#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
using namespace std;
#include "../../scr/liton_cpp_snippets/lion_snippets.hpp"
#include "../../scr/liton_point_data/PointData.hpp"

using namespace liton_pd;

int main(int argc, char** argv)
{
	string name("test");
	cout << name << endl;
	ofstream out((name + "_out.txt").c_str());
	ofstream err((name + "_err.txt").c_str());

	liton_sp::env::disp_env(out);
	out << endl;

	D2::PointData<float, 1, LO::center, LO::center> x0;
	out << x0.disp() << endl;
	D2::PointData<double, 1, LO::center, LO::half> x1;
	out << x1.disp() << endl;
	D2::PointData<double, 1, LO::half, LO::center> x2;
	out << x2.disp() << endl;
	D2::PointData<double, 2, LO::half, LO::half> x3;
	out << x3.disp() << endl;
	//D2::PointData<double, 0, LO::half, LO::half> x4;
	//D2::PointData<double, -1, LO::half, LO::half> x5;
	//D2::PointData<double, 1, LO::half, LO::half> x6(x3);
	//x3 = x3;
	out << endl;

	liton_sp::debug::exec_except([&]() {x1.alloc(0, 0, 0, 0, 0, 0); }, out, err);
	liton_sp::debug::exec_except([&]() {x1.alloc(0, 0, 0, 1, 0, 0); }, out, err);
	const unsigned big = 0.1 * 1024 / sizeof(float) / x1.N * 1024 * 1024;
	out << big* big << endl;
	liton_sp::debug::exec_except([&]() {x1.alloc(0, big, 0, 1, big, 0); }, out, err);
	liton_sp::debug::exec_except([&]() {x1.alloc(0, 10, 0, 1, 5, 0); }, out, err);
	liton_sp::debug::exec_except([&]() {x1.alloc(0, 10, 0, 1, 5, 0); }, out, err);
	out << endl;

	liton_sp::debug::exec_except([&]() {x1.realloc(1, 5, 2, 1, 10, 3); }, out, err);
	liton_sp::debug::exec_except([&]() {x2.realloc(1, 5, 2, 1, 10, 3); }, out, err);
	out << endl;

	liton_sp::debug::exec_except([&]() {x2.clear(); }, out, err);
	liton_sp::debug::exec_except([&]() {x3.clear(); x3.alloc(1, 5, 2, 1, 10, 3); }, out, err);
	liton_sp::debug::exec_except([&]() {x2.alloc(1, 5, 2, 1, 10, 3); x3.alloc(1, 5, 2, 1, 10, 3); }, out, err);
	out << endl;

	out << x3.size().disp() << endl;
	out << endl;

	out << x1.disp() << endl;
	D2::PD_For_2D(x1.size().range(RA::IN, RA::IN), [&x1]PD_F_ij(i, j) { x1(0, i, j, FL::C, FL::N) = static_cast<double>(i) + static_cast<double>(j); });
	D2::PD_For_2D(x1.size().range(RA::N, RA::N), [&x1]PD_F_ij(i, j) { x1(0, i, j, FL::C, FL::N) = 1; });
	D2::PD_For_2D(x1.size().range(RA::N, RA::P), [&x1]PD_F_ij(i, j) { x1(0, i, j, FL::C, FL::N) = 2; });
	D2::PD_For_2D(x1.size().range(RA::P, RA::N), [&x1]PD_F_ij(i, j) { x1(0, i, j, FL::C, FL::N) = 3; });
	D2::PD_For_2D(x1.size().range(RA::P, RA::P), [&x1]PD_F_ij(i, j) { x1(0, i, j, FL::C, FL::N) = 4; });
	D2::PD_For_2D(x1.size().range(RA::N, RA::IN), [&x1]PD_F_ij(i, j) { x1(0, i, j, FL::C, FL::N) = 100 + static_cast<double>(j); });
	D2::PD_For_2D(x1.size().range(RA::P, RA::IN), [&x1]PD_F_ij(i, j) { x1(0, i, j, FL::C, FL::N) = 200 + static_cast<double>(j); });
	D2::PD_For_2D(x1.size().range(RA::IN, RA::N), [&x1]PD_F_ij(i, j) { x1(0, i, j, FL::C, FL::N) = static_cast<double>(i) + 300; });
	D2::PD_For_2D(x1.size().range(RA::IN, RA::P), [&x1]PD_F_ij(i, j) { x1(0, i, j, FL::C, FL::N) = static_cast<double>(i) + 400; });
	D2::PD_For_1D(0, x1.size().range(RA::ALL, RA::ALL), [&x1]PD_F_i(i) { x1(0, i, x1.size().last(1, RA::ALL), FL::C, FL::P) = 11; });
	out << x1.disp_data() << endl;
	out << endl;

	out << x2.disp() << endl;
	D2::PD_For_2D(x2.size().range(RA::IN, RA::IN), [&x2]PD_F_ij(i, j) { x2(0, i, j, FL::N, FL::C) = static_cast<double>(i) + static_cast<double>(j); });
	D2::PD_For_2D(x2.size().range(RA::N, RA::N), [&x2]PD_F_ij(i, j) { x2(0, i, j, FL::N, FL::C) = 1; });
	D2::PD_For_2D(x2.size().range(RA::N, RA::P), [&x2]PD_F_ij(i, j) { x2(0, i, j, FL::N, FL::C) = 2; });
	D2::PD_For_2D(x2.size().range(RA::P, RA::N), [&x2]PD_F_ij(i, j) { x2(0, i, j, FL::N, FL::C) = 3; });
	D2::PD_For_2D(x2.size().range(RA::P, RA::P), [&x2]PD_F_ij(i, j) { x2(0, i, j, FL::N, FL::C) = 4; });
	D2::PD_For_2D(x2.size().range(RA::N, RA::IN), [&x2]PD_F_ij(i, j) { x2(0, i, j, FL::N, FL::C) = 100 + static_cast<double>(j); });
	D2::PD_For_2D(x2.size().range(RA::P, RA::IN), [&x2]PD_F_ij(i, j) { x2(0, i, j, FL::N, FL::C) = 200 + static_cast<double>(j); });
	D2::PD_For_2D(x2.size().range(RA::IN, RA::N), [&x2]PD_F_ij(i, j) { x2(0, i, j, FL::N, FL::C) = static_cast<double>(i) + 300; });
	D2::PD_For_2D(x2.size().range(RA::IN, RA::P), [&x2]PD_F_ij(i, j) { x2(0, i, j, FL::N, FL::C) = static_cast<double>(i) + 400; });
	D2::PD_For_1D(1, x2.size().range(RA::ALL, RA::ALL), [&x2]PD_F_i(i) { x2(0, x2.size().last(0, RA::ALL), i, FL::P, FL::C) = 11; });
	out << x2.disp_data() << endl;
	out << endl;

	out << x3.disp() << endl;
	D2::PD_For_N_2D(0, x3.N, x3.size().range(RA::IN, RA::IN), [&x3]PD_F_n_ij(n, i, j) { x3(n, i, j, FL::N, FL::N) = static_cast<double>(i) + static_cast<double>(j) + n; });
	D2::PD_For_N_2D(0, x3.N, x3.size().range(RA::N, RA::N), [&x3]PD_F_n_ij(n, i, j) { x3(n, i, j, FL::N, FL::N) = 1 + n; });
	D2::PD_For_N_2D(0, x3.N, x3.size().range(RA::N, RA::P), [&x3]PD_F_n_ij(n, i, j) { x3(n, i, j, FL::N, FL::N) = 2 + n; });
	D2::PD_For_N_2D(0, x3.N, x3.size().range(RA::P, RA::N), [&x3]PD_F_n_ij(n, i, j) { x3(n, i, j, FL::N, FL::N) = 3 + n; });
	D2::PD_For_N_2D(0, x3.N, x3.size().range(RA::P, RA::P), [&x3]PD_F_n_ij(n, i, j) { x3(n, i, j, FL::N, FL::N) = 4 + n; });
	D2::PD_For_N_2D(0, x3.N, x3.size().range(RA::N, RA::IN), [&x3]PD_F_n_ij(n, i, j) { x3(n, i, j, FL::N, FL::N) = 100 + static_cast<double>(j) + n; });
	D2::PD_For_N_2D(0, x3.N, x3.size().range(RA::P, RA::IN), [&x3]PD_F_n_ij(n, i, j) { x3(n, i, j, FL::N, FL::N) = 200 + static_cast<double>(j) + n; });
	D2::PD_For_N_2D(0, x3.N, x3.size().range(RA::IN, RA::N), [&x3]PD_F_n_ij(n, i, j) { x3(n, i, j, FL::N, FL::N) = static_cast<double>(i) + 300 + n; });
	D2::PD_For_N_2D(0, x3.N, x3.size().range(RA::IN, RA::P), [&x3]PD_F_n_ij(n, i, j) { x3(n, i, j, FL::N, FL::N) = static_cast<double>(i) + 400 + n; });
	D2::PD_For_N_1D(0, x3.N, 0, x3.size().range(RA::ALL, RA::ALL), [&x3]PD_F_n_i(n, i) { x3(n, i, x3.size().last(1, RA::ALL), FL::N, FL::P) = 11 + n; });
	D2::PD_For_N_1D(0, x3.N, 1, x3.size().range(RA::ALL, RA::ALL), [&x3]PD_F_n_i(n, i) { x3(n, x3.size().last(0, RA::ALL), i, FL::P, FL::N) = 22 + n; });
	out << x3.disp_data() << endl;
	out << endl;

	liton_sp::debug::exec_except([&]() {out << x1(0, -1, 0, FL::C, FL::P) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x1(0, -2, 0, FL::C, FL::P) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x1(0, 7, 0, FL::C, FL::P) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x1(0, 0, -2, FL::C, FL::P) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x1(0, 0, 13, FL::C, FL::P) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x1(0, 7, 13, FL::C, FL::P) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x1(0, -1, 0, FL::N, FL::N) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x1(0, -1, 0, FL::C, FL::C) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x2(0, -1, 0, FL::N, FL::N) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x2(0, -1, 0, FL::C, FL::C) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x3(0, -1, 0, FL::C, FL::N) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x3(0, -1, 0, FL::N, FL::C) << endl; }, out, err);
	const D2::PointData<double, 1, LO::center, LO::half> xc(0, 5, 0, 0, 10, 0);
	liton_sp::debug::exec_except([&]() {out << xc(0, 0, 0, FL::C, FL::N) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << xc(0, 0, 0) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x3(0, -2, 0, FL::N, FL::N) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x3(0, 7, 0, FL::P, FL::N) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x3(0, 0, -2, FL::N, FL::N) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x3(0, 0, 13, FL::N, FL::P) << endl; }, out, err);
	out << endl;

	double min = 10000;
	D2::PD_Reduce_2D(x1.size().range(RA::IN, RA::IN), min,
	[]PD_RF(double, x, xx) { xx = x < xx ? x : xx; },
	[&x1]PD_F_ij(i, j)->double { return x1(0, i, j, FL::C, FL::P) - 10; }
	                );
	out << "min " << min << endl;

	double max = 0;
	D2::PD_Reduce_N_2D(0, x3.N, x3.size().range(RA::IN, RA::IN), max,
	[]PD_RF(double, x, xx) { xx = x > xx ? x : xx; },
	[&x3]PD_F_n_ij(n, i, j)->double { return x3(n, i, j, FL::P, FL::P); }
	                  );
	out << "max " << max << endl;

	double sum_1 = 0;
	D2::PD_Reduce_1D(1, x2.size().range(RA::IN, RA::IN), sum_1,
	[]PD_RF(double, x, xx) { xx += x; },
	[&x2]PD_F_i(i)->double { return x2(0, 0, i, FL::N, FL::C) + x2(0, 0, i, FL::P, FL::C); }
	                );
	out << "sum_1 " << sum_1 << endl;

	double sum_2 = 0;
	D2::PD_Reduce_N_1D(0, x3.N, 1, x3.size().range(RA::ALL, RA::ALL), sum_2,
	[]PD_RF(double, x, xx) { xx += x; },
	[&x3]PD_F_n_i(n, i)->double { return x3(n, -1, i, FL::N, FL::N) + x3(n, -1, i, FL::N, FL::P); }
	                  );
	out << "sum_2 " << sum_2 << endl;

	out << x2.disp() << endl;
	out << x2.disp_data() << endl;
	liton_sp::debug::exec_except([&]() {x2.copy_from(x2, RA::IN, RA::IN, RA::P, RA::IN, FL::P, FL::P);}, out, err);
	liton_sp::debug::exec_except([&]() {x2.copy_from(x2, RA::IN, RA::IN, RA::P, RA::IN, FL::N, FL::P);}, out, err);
	out << x2.disp_data() << endl;

	out.close();
	err.close();

	return 0;
}