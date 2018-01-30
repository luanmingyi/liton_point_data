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

	D2::PointData<float, LO::center, LO::center, 1> x1;
	out << x1.disp() << endl;
	out << x1.disp_data() << endl;
	D2::PointData<double, LO::center, LO::center, 2> x2;
	out << x2.disp() << endl;
	D2::PointData<double, LO::center, LO::center, 2> x3(5, 10, 1, 1, 2, 3);
	out << x3.disp() << endl;
	//D2::PointData_C<double, 0> x4;
	//D2::PointData_C<double, -1> x5;
	//D2::PointData_C<float> x6(x1);
	//x3 = x2;
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

	D2::For_PD_2D(x1.size().range(RA::IN, RA::IN), [&x1]PD_F_ij(i, j) { x1(i, j, 0, FL::C, FL::C) = static_cast<float>(i) + static_cast<float>(j); });
	D2::For_PD_2D(x1.size().range(RA::N, RA::N), [&x1]PD_F_ij(i, j) { x1(i, j, 0, FL::C, FL::C) = 1; });
	D2::For_PD_2D(x1.size().range(RA::N, RA::P), [&x1]PD_F_ij(i, j) { x1(i, j, 0, FL::C, FL::C) = 2; });
	D2::For_PD_2D(x1.size().range(RA::P, RA::N), [&x1]PD_F_ij(i, j) { x1(i, j, 0, FL::C, FL::C) = 3; });
	D2::For_PD_2D(x1.size().range(RA::P, RA::P), [&x1]PD_F_ij(i, j) { x1(i, j, 0, FL::C, FL::C) = 4; });
	D2::For_PD_2D(x1.size().range(RA::N, RA::IN), [&x1]PD_F_ij(i, j) { x1(i, j, 0, FL::C, FL::C) = 100 + static_cast<float>(j); });
	D2::For_PD_2D(x1.size().range(RA::P, RA::IN), [&x1]PD_F_ij(i, j) { x1(i, j, 0, FL::C, FL::C) = 200 + static_cast<float>(j); });
	D2::For_PD_2D(x1.size().range(RA::IN, RA::N), [&x1]PD_F_ij(i, j) { x1(i, j, 0, FL::C, FL::C) = static_cast<float>(i) + 300; });
	D2::For_PD_2D(x1.size().range(RA::IN, RA::P), [&x1]PD_F_ij(i, j) { x1(i, j, 0, FL::C, FL::C) = static_cast<float>(i) + 400; });
	out << x1.disp_data() << endl;
	out << endl;

	liton_sp::debug::exec_except([&]() {out << x1(-1, 0, 0, FL::C, FL::C) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x1(-2, 0, 0, FL::C, FL::C) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x1(7, 0, 0, FL::C, FL::C) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x1(0, -2, 0, FL::C, FL::C) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x1(0, 13, 0, FL::C, FL::C) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x1(7, 13, 0, FL::C, FL::C) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x1(-1, 0, 0, FL::N, FL::C) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x1(-1, 0, 0, FL::C, FL::P) << endl; }, out, err);
	out << endl;

	D2::For_PD_2D_N(x2.size().range(RA::ALL, RA::ALL), 0, x2.N, [&x1, &x2]PD_F_ij_n(i, j, n)
	{
		x2(i, j, n, FL::C, FL::C) = x1(i, j, 0, FL::C, FL::C) * static_cast<double>(n + 1) + 5;
	});
	out << x2.disp_data() << endl;
	out << endl;

	D2::For_PD_1D(x1.size().range(RA::N, RA::ALL), 1, [&x1]PD_F_i(i)
	{
		x1(-1, i, 0, FL::C, FL::C) = 1;
		x1(x1.size().end(0, RA::ALL) - 1, i, 0, FL::C, FL::C) = 2;
	});
	D2::For_PD_1D(x1.size().range(RA::IN, RA::N), 0, [&x1]PD_F_i(i)
	{
		x1(i, -1, 0, FL::C, FL::C) = 3;
		x1(i, x1.size().end(1, RA::ALL) - 1, 0, FL::C, FL::C) = 4;
	});
	out << x1.disp_data() << endl;
	out << endl;

	D2::For_PD_1D_N(x2.size().range(RA::N, RA::ALL), 1, 0, x2.N, [&x2]PD_F_i_n(i, n)
	{
		x2(-1, i, n, FL::C, FL::C) = 10 + n + i;
		x2(x2.size().end(0, RA::ALL) - 1, i, n, FL::C, FL::C) = 20 + n + i;
	});
	D2::For_PD_1D_N(x2.size().range(RA::IN, RA::ALL), 0, 0, x2.N, [&x2]PD_F_i_n(i, n)
	{
		x2(i, -1, n, FL::C, FL::C) = 30 + n + i;
		x2(i, x2.size().end(1, RA::ALL) - 1, n, FL::C, FL::C) = 40 + n + i;
	});
	out << x2.disp_data() << endl;
	out << endl;

	float min = 10000;
	D2::Reduce_PD_2D(x1.size().range(RA::IN, RA::IN), min,
	[]PD_RF(float, x, xx) { xx = x < xx ? x : xx; },
	[&x1]PD_F_ij(i, j)->float { return x1(i, j, 0, FL::C, FL::C) - 10; }
	                );
	out << "min " << min << endl;

	double max = 0;
	D2::Reduce_PD_2D_N(x2.size().range(RA::IN, RA::IN), 0, x2.N, max,
	[]PD_RF(double, x, xx) { xx = x > xx ? x : xx; },
	[&x2]PD_F_ij_n(i, j, n)->double { return x2(i, j, n, FL::C, FL::C); }
	                  );
	out << "max " << max << endl;

	float sum_1 = 0;
	D2::Reduce_PD_1D(x1.size().range(RA::IN, RA::IN), 1, sum_1,
	[]PD_RF(float, x, xx) { xx += x; },
	[&x1]PD_F_i(i)->float { return x1(-1, i, 0, FL::C, FL::C) + x1(0, i, 0, FL::C, FL::C); }
	                );
	out << "sum_1 " << sum_1 << endl;

	double sum_2 = 0;
	D2::Reduce_PD_1D_N(x2.size().range(RA::IN, RA::IN), 1, 0, x2.N, sum_2,
	[]PD_RF(double, x, xx) { xx += x; },
	[&x2]PD_F_i_n(i, n)->double { return x2(-1, i, n, FL::C, FL::C) + x2(0, i, n, FL::C, FL::C); }
	                  );
	out << "sum_2 " << sum_2 << endl;

	out.close();
	err.close();

	return 0;
}