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

	D1::PointData_H<float> x1;
	out << x1.disp() << endl;
	D1::PointData_H<double, 2> x2;
	out << x2.disp() << endl;
	D1::PointData_H<double, 2> x3(10, 2, 3);
	out << x3.disp() << endl;
	//D1::PointData_C<double, 0> x4;
	//D1::PointData_C<double, -1> x5;
	//D1::PointData_C<float, 1> x6(x1);
	//x3 = x2;
	out << endl;

	liton_sp::debug::exec_except([&]() {x1.alloc(0, 0, 0); }, out, err);
	liton_sp::debug::exec_except([&]() {x1.alloc(0, 1, 0); }, out, err);
	out << 3u * 1024u * 1024u * 1024u / sizeof(float) / x1.N << endl;
	liton_sp::debug::exec_except([&]() {x1.alloc(3u * 1024u * 1024u * 1024u / sizeof(float) / x1.N, 2, 3); }, out, err);
	liton_sp::debug::exec_except([&]() {x1.alloc(10, 5, 5); }, out, err);
	liton_sp::debug::exec_except([&]() {x1.alloc(10, 2, 3); }, out, err);
	out << endl;

	liton_sp::debug::exec_except([&]() {x1.realloc(10, 2, 3); }, out, err);
	liton_sp::debug::exec_except([&]() {x2.realloc(10, 2, 3); }, out, err);
	out << endl;

	liton_sp::debug::exec_except([&]() {x2.clear(); }, out, err);
	liton_sp::debug::exec_except([&]() {x3.clear(); x3.alloc(10, 2, 3); }, out, err);
	liton_sp::debug::exec_except([&]() {x2.alloc(10, 2, 3); x3.alloc(10, 2, 3); }, out, err);
	out << endl;

	out << x3.size().disp() << endl;
	out << endl;

	D1::For_PD_1D(x1.size().range(RA::IN), [&x1]PD_F_i(i) { x1(i, 0, FL::N) = static_cast<float>(i); });
	D1::For_PD_1D(x1.size().range(RA::N), [&x1]PD_F_i(i) { x1(i, 0, FL::N) = static_cast<float>(i - 5); });
	D1::For_PD_1D(x1.size().range(RA::P), [&x1]PD_F_i(i) { x1(i, 0, FL::N) = static_cast<float>(i + 6); });
	x1(x1.size().end(0, RA::ALL) - 1, 0, FL::P) = -10;
	out << x1.disp_data() << endl;

	liton_sp::debug::exec_except([&]() {out << x1(-1, 0, FL::N) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x1(-5, 0, FL::P) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x1(120, 0, FL::N) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x1(-1, 1, FL::P) << endl; }, out, err);
	out << endl;

	D1::For_PD_1D_N(x2.size().range(RA::ALL), 0, x2.N, [&x1, &x2, &x3]PD_F_i_n(i, n)
	{
		x2(i, n, FL::N) = (x1(i, 0, FL::N) + x1(i, 0, FL::P)) / 2;
		x3(i, n, FL::N) = x2(i, n, FL::N) + x1(i, 0, FL::P);
	});
	For_N(0, x2.N, [&x2]PD_F_n(n) { x2(x2.size().end(0, RA::ALL) - 1, n, FL::P) = -5; });
	For_N(0, x3.N, [&x3]PD_F_n(n) { x3(x3.size().end(0, RA::ALL) - 1, n, FL::P) = 100; });
	out << x2.disp_data() << endl;
	out << x3.disp_data() << endl;
	out << endl;

	double sum = 0;
	double max = -100;
	auto lam_sum = []PD_RF(double, x, xx) { xx += x; };
	auto lam_max = []PD_RF(double, x, xx) { xx = x > xx ? x : xx; };
	D1::Reduce_PD_1D(x3.size().range(RA::IN),
	                 lam_sum,
	                 [&x3]PD_F_i(i)->double { return x3(i, 0, FL::N); },
	                 sum);
	lam_sum(x3(x3.size().end(0, RA::IN) - 1, 0, FL::P), sum);
	D1::Reduce_PD_1D_N(x3.size().range(RA::ALL),
	                   0, x3.N,
	                   lam_max,
	                   [&x3]PD_F_i_n(i, n)->double { return x3(i, n, FL::N); },
	                   max);
	Reduce_N(0, x3.N, lam_max, [&x3]PD_F_n(n)->double { return x3(x3.size().end(0, RA::ALL) - 1, n, FL::P); }, max);
	out << sum << endl;
	out << max << endl;

	out.close();
	err.close();

	return 0;
}