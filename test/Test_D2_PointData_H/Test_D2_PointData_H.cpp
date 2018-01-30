#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
using namespace std;
#include "../dep/liton_cpp_snippets/lion_snippets.hpp"

#ifdef _DEBUG
#define _CHECK_POINTDATA_RANGE
#endif
#include "../../scr/PointData.hpp"

using namespace liton_pd;

int main(int argc, char** argv)
{
	string name(__FILE__);
	name.erase(name.find_last_of('.'));
	cout << name << endl;
	ofstream out((name + "_out.txt").c_str());
	ofstream err((name + "_err.txt").c_str());

	liton_sp::env::disp_env(out);
	out << endl;

	D2::PointData<float, LO::center, LO::center, 1> x0;
	out << x0.disp() << endl;
	D2::PointData<double, LO::center, LO::half, 1> x1;
	out << x1.disp() << endl;
	D2::PointData<double, LO::half, LO::center, 1> x2;
	out << x2.disp() << endl;
	D2::PointData<double, LO::half, LO::half, 2> x3;
	out << x3.disp() << endl;
	//D2::PointData<double, LO::half, LO::half, 0> x4;
	//D2::PointData<double, LO::half, LO::half, -1> x5;
	//D2::PointData<double, LO::half, LO::half, 1> x6(x3);
	//x3 = x3;
	out << endl;

	liton_sp::debug::exec_except([&]() {x1.alloc(0, 0, 0, 0, 0, 0); }, out, err);
	liton_sp::debug::exec_except([&]() {x1.alloc(0, 0, 0, 1, 0, 0); }, out, err);
	const unsigned big = 0.1 * 1024 / sizeof(float) / x1.N * 1024 * 1024;
	out << big* big << endl;
	liton_sp::debug::exec_except([&]() {x1.alloc(big, big, 0, 1, 0, 0); }, out, err);
	liton_sp::debug::exec_except([&]() {x1.alloc(10, 5, 0, 1, 0, 0); }, out, err);
	liton_sp::debug::exec_except([&]() {x1.alloc(10, 5, 0, 1, 0, 0); }, out, err);
	out << endl;

	liton_sp::debug::exec_except([&]() {x1.realloc(5, 10, 1, 1, 2, 3); }, out, err);
	liton_sp::debug::exec_except([&]() {x2.realloc(5, 10, 1, 1, 2, 3); }, out, err);
	out << endl;

	liton_sp::debug::exec_except([&]() {x2.clear(); }, out, err);
	liton_sp::debug::exec_except([&]() {x3.clear(); x3.alloc(5, 10, 1, 1, 2, 3); }, out, err);
	liton_sp::debug::exec_except([&]() {x2.alloc(5, 10, 1, 1, 2, 3); x3.alloc(5, 10, 1, 1, 2, 3); }, out, err);
	out << endl;

	out << x3.size().disp() << endl;
	out << endl;

	out << x1.disp() << endl;
	D2::For_PD_2D(x1.size().range(RA::IN, RA::IN), [&x1]PD_F_ij(i, j) { x1(i, j, 0, FL::C, FL::N) = static_cast<double>(i) + static_cast<double>(j); });
	D2::For_PD_2D(x1.size().range(RA::N, RA::N), [&x1]PD_F_ij(i, j) { x1(i, j, 0, FL::C, FL::N) = 1; });
	D2::For_PD_2D(x1.size().range(RA::N, RA::P), [&x1]PD_F_ij(i, j) { x1(i, j, 0, FL::C, FL::N) = 2; });
	D2::For_PD_2D(x1.size().range(RA::P, RA::N), [&x1]PD_F_ij(i, j) { x1(i, j, 0, FL::C, FL::N) = 3; });
	D2::For_PD_2D(x1.size().range(RA::P, RA::P), [&x1]PD_F_ij(i, j) { x1(i, j, 0, FL::C, FL::N) = 4; });
	D2::For_PD_2D(x1.size().range(RA::N, RA::IN), [&x1]PD_F_ij(i, j) { x1(i, j, 0, FL::C, FL::N) = 100 + static_cast<double>(j); });
	D2::For_PD_2D(x1.size().range(RA::P, RA::IN), [&x1]PD_F_ij(i, j) { x1(i, j, 0, FL::C, FL::N) = 200 + static_cast<double>(j); });
	D2::For_PD_2D(x1.size().range(RA::IN, RA::N), [&x1]PD_F_ij(i, j) { x1(i, j, 0, FL::C, FL::N) = static_cast<double>(i) + 300; });
	D2::For_PD_2D(x1.size().range(RA::IN, RA::P), [&x1]PD_F_ij(i, j) { x1(i, j, 0, FL::C, FL::N) = static_cast<double>(i) + 400; });
	D2::For_PD_1D(x1.size().range(RA::ALL, RA::ALL), 0, [&x1]PD_F_i(i) { x1(i, x1.size().last(1, RA::ALL), 0, FL::C, FL::P) = 11; });
	out << x1.disp_data() << endl;
	out << endl;

	out << x2.disp() << endl;
	D2::For_PD_2D(x2.size().range(RA::IN, RA::IN), [&x2]PD_F_ij(i, j) { x2(i, j, 0, FL::N, FL::C) = static_cast<double>(i) + static_cast<double>(j); });
	D2::For_PD_2D(x2.size().range(RA::N, RA::N), [&x2]PD_F_ij(i, j) { x2(i, j, 0, FL::N, FL::C) = 1; });
	D2::For_PD_2D(x2.size().range(RA::N, RA::P), [&x2]PD_F_ij(i, j) { x2(i, j, 0, FL::N, FL::C) = 2; });
	D2::For_PD_2D(x2.size().range(RA::P, RA::N), [&x2]PD_F_ij(i, j) { x2(i, j, 0, FL::N, FL::C) = 3; });
	D2::For_PD_2D(x2.size().range(RA::P, RA::P), [&x2]PD_F_ij(i, j) { x2(i, j, 0, FL::N, FL::C) = 4; });
	D2::For_PD_2D(x2.size().range(RA::N, RA::IN), [&x2]PD_F_ij(i, j) { x2(i, j, 0, FL::N, FL::C) = 100 + static_cast<double>(j); });
	D2::For_PD_2D(x2.size().range(RA::P, RA::IN), [&x2]PD_F_ij(i, j) { x2(i, j, 0, FL::N, FL::C) = 200 + static_cast<double>(j); });
	D2::For_PD_2D(x2.size().range(RA::IN, RA::N), [&x2]PD_F_ij(i, j) { x2(i, j, 0, FL::N, FL::C) = static_cast<double>(i) + 300; });
	D2::For_PD_2D(x2.size().range(RA::IN, RA::P), [&x2]PD_F_ij(i, j) { x2(i, j, 0, FL::N, FL::C) = static_cast<double>(i) + 400; });
	D2::For_PD_1D(x2.size().range(RA::ALL, RA::ALL), 1, [&x2]PD_F_i(i) { x2(x2.size().last(0, RA::ALL), i, 0, FL::P, FL::C) = 11; });
	out << x2.disp_data() << endl;
	out << endl;

	out << x3.disp() << endl;
	D2::For_PD_2D_N(x3.size().range(RA::IN, RA::IN), 0, x3.N, [&x3]PD_F_ij_n(i, j, n) { x3(i, j, n, FL::N, FL::N) = static_cast<double>(i) + static_cast<double>(j) + n; });
	D2::For_PD_2D_N(x3.size().range(RA::N, RA::N), 0, x3.N, [&x3]PD_F_ij_n(i, j, n) { x3(i, j, n, FL::N, FL::N) = 1 + n; });
	D2::For_PD_2D_N(x3.size().range(RA::N, RA::P), 0, x3.N, [&x3]PD_F_ij_n(i, j, n) { x3(i, j, n, FL::N, FL::N) = 2 + n; });
	D2::For_PD_2D_N(x3.size().range(RA::P, RA::N), 0, x3.N, [&x3]PD_F_ij_n(i, j, n) { x3(i, j, n, FL::N, FL::N) = 3 + n; });
	D2::For_PD_2D_N(x3.size().range(RA::P, RA::P), 0, x3.N, [&x3]PD_F_ij_n(i, j, n) { x3(i, j, n, FL::N, FL::N) = 4 + n; });
	D2::For_PD_2D_N(x3.size().range(RA::N, RA::IN), 0, x3.N, [&x3]PD_F_ij_n(i, j, n) { x3(i, j, n, FL::N, FL::N) = 100 + static_cast<double>(j) + n; });
	D2::For_PD_2D_N(x3.size().range(RA::P, RA::IN), 0, x3.N, [&x3]PD_F_ij_n(i, j, n) { x3(i, j, n, FL::N, FL::N) = 200 + static_cast<double>(j) + n; });
	D2::For_PD_2D_N(x3.size().range(RA::IN, RA::N), 0, x3.N, [&x3]PD_F_ij_n(i, j, n) { x3(i, j, n, FL::N, FL::N) = static_cast<double>(i) + 300 + n; });
	D2::For_PD_2D_N(x3.size().range(RA::IN, RA::P), 0, x3.N, [&x3]PD_F_ij_n(i, j, n) { x3(i, j, n, FL::N, FL::N) = static_cast<double>(i) + 400 + n; });
	D2::For_PD_1D_N(x3.size().range(RA::ALL, RA::ALL), 0, 0, x3.N, [&x3]PD_F_i_n(i, n) { x3(i, x3.size().last(1, RA::ALL), n, FL::N, FL::P) = 11 + n; });
	D2::For_PD_1D_N(x3.size().range(RA::ALL, RA::ALL), 1, 0, x3.N, [&x3]PD_F_i_n(i, n) { x3(x3.size().last(0, RA::ALL), i, n, FL::P, FL::N) = 22 + n; });
	out << x3.disp_data() << endl;
	out << endl;

	liton_sp::debug::exec_except([&]() {out << x1(-1, 0, 0, FL::C, FL::P) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x1(-2, 0, 0, FL::C, FL::P) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x1(7, 0, 0, FL::C, FL::P) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x1(0, -2, 0, FL::C, FL::P) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x1(0, 13, 0, FL::C, FL::P) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x1(7, 13, 0, FL::C, FL::P) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x1(-1, 0, 0, FL::N, FL::N) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x1(-1, 0, 0, FL::C, FL::C) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x2(-1, 0, 0, FL::N, FL::N) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x2(-1, 0, 0, FL::C, FL::C) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x3(-1, 0, 0, FL::C, FL::N) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x3(-1, 0, 0, FL::N, FL::C) << endl; }, out, err);
	out << endl;

	double min = 10000;
	D2::Reduce_PD_2D(x1.size().range(RA::IN, RA::IN), min,
		[]PD_RF(double, x, xx) { xx = x < xx ? x : xx; },
		[&x1]PD_F_ij(i, j)->double { return x1(i, j, 0, FL::C, FL::P) - 10; }
	);
	out << "min " << min << endl;

	double max = 0;
	D2::Reduce_PD_2D_N(x3.size().range(RA::IN, RA::IN), 0, x3.N, max,
		[]PD_RF(double, x, xx) { xx = x > xx ? x : xx; },
		[&x3]PD_F_ij_n(i, j, n)->double { return x3(i, j, n, FL::P, FL::P); }
	);
	out << "max " << max << endl;

	double sum_1 = 0;
	D2::Reduce_PD_1D(x2.size().range(RA::IN, RA::IN), 1, sum_1,
		[]PD_RF(double, x, xx) { xx += x; },
		[&x2]PD_F_i(i)->double { return x2(0, i, 0, FL::N, FL::C) + x2(0, i, 0, FL::P, FL::C); }
	);
	out << "sum_1 " << sum_1 << endl;

	double sum_2 = 0;
	D2::Reduce_PD_1D_N(x3.size().range(RA::ALL, RA::ALL), 1, 0, x3.N, sum_2,
		[]PD_RF(double, x, xx) { xx += x; },
		[&x3]PD_F_i_n(i, n)->double { return x3(-1, i, n, FL::N, FL::N) + x3(-1, i, n, FL::N, FL::P); }
	);
	out << "sum_2 " << sum_2 << endl;

	out.close();
	err.close();

	return 0;
}