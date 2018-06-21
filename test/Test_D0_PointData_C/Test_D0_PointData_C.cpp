#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
using namespace std;
#include "../../scr/liton_cpp_snippets/lion_snippets.hpp"

#ifdef _DEBUG
	#define _CHECK_POINTDATA_RANGE
#endif
#include "../../scr/liton_point_data/PointData.hpp"

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

	D0::PointData<float, 1> x1;
	out << x1.disp() << endl;
	D0::PointData<double, 10> x2;
	out << x2.disp() << endl;
	//D0::PointData<double, 0> x3;
	//D0::PointData<double, -1> x4;
	//D0::PointData<float, 1> x5(x1);
	//x1 = x2;
	out << endl;

	PD_For_N(0, x1.N, [&x1]PD_F_n(n) { x1(n) = 10; });
	out << x1.disp_data() << endl;
	PD_For_N(0, x2.N, [&x2]PD_F_n(n) { x2(n) = n; });
	out << x2.disp_data() << endl;
	out << endl;

	liton_sp::debug::exec_except([&]() {out << x1(0) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x1(1) << endl; }, out, err);

	const D0::PointData<double, 1> xc;
	liton_sp::debug::exec_except([&]() {out << xc(0) << endl; }, out, err);
	out << endl;

	double sum = 0;
	PD_Reduce_N(0, x2.N, sum,
	[]PD_RF(double, x, xx) { xx += x; },
	[&x2]PD_F_n(n)->double { return x2(n) / 2; });
	out << sum << endl;

	out.close();
	err.close();

	return 0;
}