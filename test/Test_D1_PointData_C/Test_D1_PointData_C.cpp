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

using namespace liton;

template<typename Fun>
void exec(Fun f, ostream &out_file, ostream &except_file)
{
	try
	{
		f();
		out_file << "no exception" << endl;
	}
	catch (const std::exception &err)
	{
		out_file << "with exception" << endl;
		except_file << err.what() << endl;
	}
}

int main(int argc, char** argv)
{
	string name(__FILE__);
	name.erase(name.find_last_of('.'));
	cout << name << endl;
	ofstream out((name + "_out.txt").c_str());
	ofstream err((name + "_err.txt").c_str());

	liton_sp::env::disp_env(out);
	out << endl;

	D1::PointData_C<float> x1;
	out << x1.disp() << endl;
	D1::PointData_C<double, 2> x2;
	out << x2.disp() << endl;
	D1::PointData_C<double, 2> x3(10, 2, 3);
	out << x3.disp() << endl;
	//D1::PointData_C<double, 0> x4;
	//D1::PointData_C<double, -1> x5;
	//D1::PointData_C<float, 1> x6(x1);
	//x3 = x2;
	out << endl;

	exec([&]() {x1.alloc(0, 0, 0); }, out, err);
	exec([&]() {x1.alloc(0, 1, 0); }, out, err);
	out << 3 * 1024 / sizeof(float) / x1.N * 1024 * 1024  << endl;
	exec([&]() {x1.alloc(3 * 1024 / sizeof(float) / x1.N * 1024 * 1024, 2, 3); }, out, err);
	exec([&]() {x1.alloc(10, 5, 5); }, out, err);
	exec([&]() {x1.alloc(10, 2, 3); }, out, err);
	out << endl;

	exec([&]() {x1.realloc(10, 2, 3); }, out, err);
	exec([&]() {x2.realloc(10, 2, 3); }, out, err);
	out << endl;

	exec([&]() {x2.clear(); }, out, err);
	exec([&]() {x3.clear(); x3.alloc(10, 2, 3); }, out, err);
	exec([&]() {x2.alloc(10, 2, 3); x3.alloc(10, 2, 3); }, out, err);
	out << endl;

	out << x3.size().disp() << endl;
	out << endl;

	D1::For_PD_1D(x1.size().range(RA::IN), [&x1]PD_F_i(i) { x1(i, 0) = static_cast<float>(i); });
	D1::For_PD_1D(x1.size().range(RA::N), [&x1]PD_F_i(i) { x1(i, 0) = static_cast<float>(i - 5); });
	D1::For_PD_1D(x1.size().range(RA::P), [&x1]PD_F_i(i) { x1(i, 0) = static_cast<float>(i + 6); });
	out << x1.disp_data() << endl;

	exec([&]() {out << x1(-1, 0) << endl; }, out, err);
	exec([&]() {out << x1(-5, 0) << endl; }, out, err);
	exec([&]() {out << x1(120, 0) << endl; }, out, err);
	exec([&]() {out << x1(-1, 1) << endl; }, out, err);
	out << endl;

	D1::For_PD_1D_N(x2.size().range(RA::ALL), 0, x2.N, [&x1, &x2, &x3]PD_F_i_n(i, n)
	{
		x2(i, n) = static_cast<double>(i) * static_cast<double>(n + 1);
		x3(i, n) = x2(i, n) + x1(i, 0);
	});
	out << x2.disp_data() << endl;
	out << x3.disp_data() << endl;
	out << endl;

	double sum = 0;
	double max = -100;
	D1::Reduce_PD_1D(x3.size().range(RA::IN),
	[]PD_RF(double, x, xx) { xx += x; },
	[&x3]PD_F_i(i)->double { return x3(i, 0) / 2; }, sum);
	D1::Reduce_PD_1D_N(x3.size().range(RA::ALL), 0, x3.N,
	[]PD_RF(double, x, xx) { xx = x > xx ? x : xx; },
	[&x3]PD_F_i_n(i, n)->double { return x3(i, n); },
	max);
	out << sum << endl;
	out << max << endl;

	out.close();
	err.close();

	return 0;
}