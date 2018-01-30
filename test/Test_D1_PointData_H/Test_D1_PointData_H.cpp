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

	D1::PointData<float, 1, LO::half> x1;
	out << x1.disp() << endl;
	D1::PointData<double, 2, LO::half> x2;
	out << x2.disp() << endl;
	D1::PointData<double, 2, LO::half> x3(2, 10, 3);
	out << x3.disp() << endl;
	//D1::PointData<double, 0, LO::half> x4;
	//D1::PointData<double, -1, LO::half> x5;
	//D1::PointData<float, 1, LO::half> x6(x1);
	//x3 = x2;
	out << endl;

	liton_sp::debug::exec_except([&]() {x1.alloc(0, 0, 0); }, out, err);
	liton_sp::debug::exec_except([&]() {x1.alloc(1, 0, 0); }, out, err);
	out << 3u * 1024u * 1024u * 1024u / sizeof(float) / x1.N << endl;
	liton_sp::debug::exec_except([&]() {x1.alloc(2, 3u * 1024u * 1024u * 1024u / sizeof(float) / x1.N, 3); }, out, err);
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

	D1::PD_For_1D(x1.size().range(RA::IN), [&x1]PD_F_i(i) { x1(0, i, FL::N) = static_cast<float>(i); });
	D1::PD_For_1D(x1.size().range(RA::N), [&x1]PD_F_i(i) { x1(0, i, FL::N) = static_cast<float>(i - 5); });
	D1::PD_For_1D(x1.size().range(RA::P), [&x1]PD_F_i(i) { x1(0, i, FL::N) = static_cast<float>(i + 6); });
	x1(0, x1.size().last(0,RA::ALL), FL::P) = -10;
	out << x1.disp_data() << endl;

	liton_sp::debug::exec_except([&]() {out << x1(-1, 0, FL::N) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x1(-5, 0, FL::P) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x1(120, 0, FL::N) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x1(-1, 1, FL::P) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x1(-1, 0, FL::C) << endl; }, out, err);
	out << endl;

	D1::PD_For_N_1D(0, x2.N, x2.size().range(RA::ALL), [&x1, &x2, &x3]PD_F_n_i(n, i)
	{
		x2(n, i, FL::N) = (x1(0, i, FL::N) + x1(0, i, FL::P)) / 2;
		x3(n, i, FL::N) = x2(n, i, FL::N) + x1(0, i, FL::P);
	});
	PD_For_N(0, x2.N, [&x2]PD_F_n(n) { x2(n, x2.size().last(0, RA::ALL), FL::P) = -5; });
	PD_For_N(0, x3.N, [&x3]PD_F_n(n) { x3(n, x3.size().last(0, RA::ALL), FL::P) = 100; });
	out << x2.disp_data() << endl;
	out << x3.disp_data() << endl;
	out << endl;

	double sum = 0;
	double max = -100;
	auto lam_sum = []PD_RF(double, x, xx) { xx += x; };
	auto lam_max = []PD_RF(double, x, xx) { xx = x > xx ? x : xx; };
	D1::PD_Reduce_1D(x3.size().range(RA::IN),
	                 sum,
	                 lam_sum,
	                 [&x3]PD_F_i(i)->double { return x3(0, i, FL::N); });
	lam_sum(x3(0, x3.size().last(0, RA::IN), FL::P), sum);
	D1::PD_Reduce_N_1D(0, x3.N, x3.size().range(RA::ALL),
	                   max,
	                   lam_max,
	                   [&x3]PD_F_n_i(n, i)->double { return x3(n, i, FL::N); });
	PD_Reduce_N(0, x3.N, max, lam_max, [&x3]PD_F_n(n)->double { return x3(n, x3.size().last(0, RA::ALL), FL::P); });
	out << sum << endl;
	out << max << endl;

	out.close();
	err.close();

	return 0;
}