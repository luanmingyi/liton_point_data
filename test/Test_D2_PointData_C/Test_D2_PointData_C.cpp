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
	string name("tset");
	cout << name << endl;
	ofstream out((name + "_out.txt").c_str());
	ofstream err((name + "_err.txt").c_str());

	liton_sp::env::disp_env(out);
	out << endl;

	out << "line: " << __LINE__ << endl;
	D2::PointData<float, 1, LO::center, LO::center> x1;
	out << x1.disp() << endl;
	out << x1.disp_data() << endl;
	D2::PointData<double, 2, LO::center, LO::center> x2;
	out << x2.disp() << endl;
	D2::PointData<double, 2, LO::center, LO::center> x3(1, 5, 2, 1, 10, 3);
	out << x3.disp() << endl;
	//D2::PointData_C<double, 0, LO::center, LO::center> x4;
	//D2::PointData_C<double, -1, LO::center, LO::center> x5;
	//D2::PointData_C<float, LO::center, LO::center> x6(x1);
	//x3 = x2;
	out << endl;

	out << "line: " << __LINE__ << endl;
	liton_sp::debug::exec_except([&]() {out << x1(0, 0, 0) << endl; }, out, err);
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

	out << "line: " << __LINE__ << endl;
	out << x3.size().disp() << endl;
	out << endl;

	out << "line: " << __LINE__ << endl;
	D2::PD_For_2D(x1.size().range(RA::IN, RA::IN), [&x1]PD_F_ij(i, j) { x1(0, i, j) = static_cast<float>(i) + static_cast<float>(j); });
	D2::PD_For_2D(x1.size().range(RA::N, RA::N), [&x1]PD_F_ij(i, j) { x1(0, i, j) = 1; });
	D2::PD_For_2D(x1.size().range(RA::N, RA::P), [&x1]PD_F_ij(i, j) { x1(0, i, j) = 2; });
	D2::PD_For_2D(x1.size().range(RA::P, RA::N), [&x1]PD_F_ij(i, j) { x1(0, i, j) = 3; });
	D2::PD_For_2D(x1.size().range(RA::P, RA::P), [&x1]PD_F_ij(i, j) { x1(0, i, j) = 4; });
	D2::PD_For_2D(x1.size().range(RA::N, RA::IN), [&x1]PD_F_ij(i, j) { x1(0, i, j) = 100 + static_cast<float>(j); });
	D2::PD_For_2D(x1.size().range(RA::P, RA::IN), [&x1]PD_F_ij(i, j) { x1(0, i, j) = 200 + static_cast<float>(j); });
	D2::PD_For_2D(x1.size().range(RA::IN, RA::N), [&x1]PD_F_ij(i, j) { x1(0, i, j) = static_cast<float>(i) + 300; });
	D2::PD_For_2D(x1.size().range(RA::IN, RA::P), [&x1]PD_F_ij(i, j) { x1(0, i, j) = static_cast<float>(i) + 400; });
	out << x1.disp_data() << endl;
	out << endl;

	out << "line: " << __LINE__ << endl;
	D2::PD_For_2D(x1.size().range(RA::IN, RA::N), [&x1]PD_F_ij(i, j) { x1(0, i, j) = -x1(0, i, x1.size().mirror(1, FL::N, j)); });
	D2::PD_For_2D(x1.size().range(RA::IN, RA::P), [&x1]PD_F_ij(i, j) { x1(0, i, j) = -x1(0, i, x1.size().mirror(1, FL::P, j)); });
	D2::PD_For_2D(x1.size().range(RA::N, RA::IN), [&x1]PD_F_ij(i, j) { x1(0, i, j) = -x1(0, x1.size().mirror(0, FL::N, i), j); });
	D2::PD_For_2D(x1.size().range(RA::P, RA::IN), [&x1]PD_F_ij(i, j) { x1(0, i, j) = -x1(0, x1.size().mirror(0, FL::P, i), j); });
	out << x1.disp_data() << endl;
	out << endl;

	out << "line: " << __LINE__ << endl;
	liton_sp::debug::exec_except([&]() {out << x1(0, -1, 0, FL::C, FL::C) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x1(0, -1, 0) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x1(0, -2, 0) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x1(0, 7, 0) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x1(0, 0, -2) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x1(0, 0, 13) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x1(0, 7, 13) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x1(0, -1, 0, FL::N, FL::C) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x1(0, -1, 0, FL::C, FL::P) << endl; }, out, err);
	const D2::PointData<double, 1, LO::center, LO::center> xc(0, 5, 0, 0, 10, 0);
	liton_sp::debug::exec_except([&]() {out << xc(0, 0, 0, FL::C, FL::C) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << xc(0, 0, 0) << endl; }, out, err);
	out << endl;

	out << "line: " << __LINE__ << endl;
	D2::PD_For_N_2D(0, x2.N, x2.size().range(RA::ALL, RA::ALL), [&x1, &x2]PD_F_n_ij(n, i, j)
	{
		x2(n, i, j) = x1(0, i, j) * static_cast<double>(n + 1) + 5;
	});
	out << x2.disp_data() << endl;
	out << endl;

	out << "line: " << __LINE__ << endl;
	D2::PD_For_1D(1, x1.size().range(RA::N, RA::ALL), [&x1]PD_F_i(i)
	{
		x1(0, -1, i) = 1;
		x1(0, x1.size().last(0, RA::ALL), i) = 2;
	});
	D2::PD_For_1D(0, x1.size().range(RA::IN, RA::N), [&x1]PD_F_i(i)
	{
		x1(0, i, -1) = 3;
		x1(0, i, x1.size().last(1, RA::ALL)) = 4;
	});
	out << x1.disp_data() << endl;
	out << endl;

	out << "line: " << __LINE__ << endl;
	D2::PD_For_N_1D(0, x2.N, 1, x2.size().range(RA::N, RA::ALL), [&x2]PD_F_n_i(n, i)
	{
		x2(n, -1, i) = 10 + n + i;
		x2(n, x2.size().last(0, RA::ALL), i) = 20 + n + i;
	});
	D2::PD_For_N_1D(0, x2.N, 0, x2.size().range(RA::IN, RA::ALL), [&x2]PD_F_n_i(n, i)
	{
		x2(n, i, -1) = 30 + n + i;
		x2(n, i, x2.size().last(1, RA::ALL)) = 40 + n + i;
	});
	out << x2.disp_data() << endl;
	out << endl;

	out << "line: " << __LINE__ << endl;
	float min = 10000;
	D2::PD_Reduce_2D(x1.size().range(RA::IN, RA::IN), min,
	[]PD_RF(float, x, xx) { xx = x < xx ? x : xx; },
	[&x1]PD_F_ij(i, j)->float { return x1(0, i, j) - 10; }
	                );
	out << "min " << min << endl;

	out << "line: " << __LINE__ << endl;
	double max = 0;
	D2::PD_Reduce_N_2D(0, x2.N, x2.size().range(RA::IN, RA::IN), max,
	[]PD_RF(double, x, xx) { xx = x > xx ? x : xx; },
	[&x2]PD_F_n_ij(n, i, j)->double { return x2(n, i, j); }
	                  );
	out << "max " << max << endl;

	out << "line: " << __LINE__ << endl;
	float sum_1 = 0;
	D2::PD_Reduce_1D(1, x1.size().range(RA::IN, RA::IN), sum_1,
	[]PD_RF(float, x, xx) { xx += x; },
	[&x1]PD_F_i(i)->float { return x1(0, -1, i) + x1(0, 0, i); }
	                );
	out << "sum_1 " << sum_1 << endl;

	out << "line: " << __LINE__ << endl;
	double sum_2 = 0;
	D2::PD_Reduce_N_1D(0, x2.N, 1, x2.size().range(RA::IN, RA::IN), sum_2,
	[]PD_RF(double, x, xx) { xx += x; },
	[&x2]PD_F_n_i(n, i)->double { return x2(n, -1, i) + x2(n, 0, i); }
	                  );
	out << "sum_2 " << sum_2 << endl;

	out << endl;
	out << x2.disp() << endl;
	out << x2.disp_data() << endl;
	liton_sp::debug::exec_except([&]() {x2.copy_from(x2, RA::IN, RA::IN, RA::IN, RA::P, FL::P, FL::P);}, out, err);
	out << x2.disp_data() << endl;

	out.close();
	err.close();

	return 0;
}