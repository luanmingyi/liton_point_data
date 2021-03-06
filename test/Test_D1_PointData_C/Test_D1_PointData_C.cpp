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

	D1::PointData<float, 1, LO::center> x1;
	out << x1.disp() << endl;
	D1::PointData<double, 2, LO::center> x2;
	out << x2.disp() << endl;
	D1::PointData<double, 2, LO::center> x3(2, 10, 3);
	out << x3.disp() << endl;
	//D1::PointData<double, 0, LO::center> x4;
	//D1::PointData<double, -1, LO::center> x5;
	//D1::PointData<float, 1, LO::center> x6(x1);
	//x3 = x2;
	out << endl;

	liton_sp::debug::exec_except([&]() {out << x1(0, 0) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {x1.alloc(0, 0, 0); }, out, err);
	liton_sp::debug::exec_except([&]() {x1.alloc(1, 0, 0); }, out, err);
	out << 3 * 1024 / sizeof(float) / x1.N * 1024 * 1024  << endl;
	liton_sp::debug::exec_except([&]() {x1.alloc(2, 3 * 1024 / sizeof(float) / x1.N * 1024 * 1024, 3); }, out, err);
	liton_sp::debug::exec_except([&]() {x1.alloc(5, 10, 5); }, out, err);
	liton_sp::debug::exec_except([&]() {x1.alloc(2, 10, 3); }, out, err);
	out << endl;

	liton_sp::debug::exec_except([&]() {x1.realloc(2, 10, 3); }, out, err);
	liton_sp::debug::exec_except([&]() {x2.realloc(2, 10, 3); }, out, err);
	out << endl;

	liton_sp::debug::exec_except([&]() {x2.clear(); }, out, err);
	liton_sp::debug::exec_except([&]() {x3.clear(); x3.alloc(2, 10, 3); }, out, err);
	liton_sp::debug::exec_except([&]() {x2.alloc(2, 10, 3); x3.alloc(2, 10, 3); }, out, err);
	out << endl;

	out << x3.size().disp() << endl;
	out << endl;

	D1::PD_For_1D(x1.size().range(RA::IN), [&x1]PD_F_i(i) { x1(0, i) = static_cast<float>(i); });
	D1::PD_For_1D(x1.size().range(RA::N), [&x1]PD_F_i(i) { x1(0, i) = static_cast<float>(i - 5); });
	D1::PD_For_1D(x1.size().range(RA::P), [&x1]PD_F_i(i) { x1(0, i) = static_cast<float>(i + 6); });
	out << x1.disp_data() << endl;

	D1::PD_For_1D(x1.size().range(RA::N), [&x1]PD_F_i(i) { x1(0, i) = -x1(0, x1.size().mirror(0, FL::N, i)); });
	D1::PD_For_1D(x1.size().range(RA::P), [&x1]PD_F_i(i) { x1(0, i) = -x1(0, x1.size().mirror(0, FL::P, i)); });
	out << x1.disp_data() << endl;

	liton_sp::debug::exec_except([&]() {out << x1(0, -1, FL::C) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x1(0, -1) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x1(0, -5) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x1(0, 13) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x1(1, -1) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x1(0, -1) << endl; }, out, err);
	const D1::PointData<double, 1, LO::center> xc(0, 2, 0);
	liton_sp::debug::exec_except([&]() {out << xc(0, 0, FL::C) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << xc(0, 0) << endl; }, out, err);
	out << endl;

	D1::PD_For_N_1D(0, x2.N, x2.size().range(RA::ALL), [&x1, &x2, &x3]PD_F_n_i(n, i)
	{
		x2(n, i) = static_cast<double>(i) * static_cast<double>(n + 1);
		x3(n, i) = x2(n, i) + x1(0, i);
	});
	out << x2.disp_data() << endl;
	out << x3.disp_data() << endl;
	out << endl;

	double sum = 0;
	double max = -100;
	D1::PD_Reduce_1D(x3.size().range(RA::IN), sum,
	[]PD_RF(double, x, xx) { xx += x; },
	[&x3]PD_F_i(i)->double { return x3(0, i) / 2; });
	D1::PD_Reduce_N_1D(0, x3.N, x3.size().range(RA::ALL), max,
	[]PD_RF(double, x, xx) { xx = x > xx ? x : xx; },
	[&x3]PD_F_n_i(n, i)->double { return x3(n, i); });
	out << sum << endl;
	out << max << endl;

	out << endl;
	liton_sp::debug::exec_except([&]() {x3.copy_from(x2, RA::ALL, RA::ALL, 10, 0);}, out, err);
	out << x3.disp_data() << endl;
	liton_sp::debug::exec_except([&]() {x3.copy_from(x2, RA::ALL, RA::ALL, 0, 10);}, out, err);
	out << x3.disp_data() << endl;
	liton_sp::debug::exec_except([&]() {x3.copy_from(x2, RA::IN, RA::N, 0, -2);}, out, err);
	out << x3.disp_data() << endl;
	liton_sp::debug::exec_except([&]() {x3.copy_from(x2, RA::IN, RA::P, 9, 12);}, out, err);
	out << x3.disp_data() << endl;
	liton_sp::debug::exec_except([&]() {x3.copy_from(x3, RA::N, RA::P, 0, 12);}, out, err);
	out << x3.disp_data() << endl;
	liton_sp::debug::exec_except([&]() {x3.copy_from(x3, RA::ALL, RA::IN, 0, 0);}, out, err);
	out << x3.disp_data() << endl;
	liton_sp::debug::exec_except([&]() {x3.copy_from(x3, RA::IN, RA::P, FL::P);}, out, err);
	out << x3.disp_data() << endl;

	out.close();
	err.close();

	return 0;
}