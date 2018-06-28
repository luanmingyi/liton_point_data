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

	D1::RangeT a;
	out << a.disp() << endl;
	D1::RangeT b(2, 10);
	out << b.disp() << endl;
	out << endl;

	liton_sp::debug::exec_except([&]() {out << a.begin(0) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << a.end(0) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << a.size(0) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << a.last(0) << endl; }, out, err);
	out << endl;
	liton_sp::debug::exec_except([&]() {out << b.begin(0) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << b.end(0) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << b.size(0) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << b.last(0) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << b.bound(0, FL::N) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << b.bound(0, FL::P) << endl; }, out, err);
	out << endl;
	liton_sp::debug::exec_except([&]() {out << b.size(1) << endl; }, out, err);
	out << endl;

	out << b.disp() << endl;
	b.cut_head(1);
	out << b.disp() << endl;
	b.cut_tail(1);
	out << b.disp() << endl;
	D1::RangeT c = b;
	c.tran(2);
	out << c.disp() << endl;
	out << b.overlap(c).disp() << endl;
	liton_sp::debug::exec_except([&]() {out << b.overlap(c.tran(6)).disp() << endl;}, out, err);
	liton_sp::debug::exec_except([&]() {out << b.overlap(c.tran(-16)).disp() << endl;}, out, err);
	liton_sp::debug::exec_except([&]() {out << b.overlap(c.tran(1)).disp() << endl;}, out, err);
	out << endl;

	liton_sp::debug::exec_except([&]() {b.check_range(0, 2); }, out, err);
	liton_sp::debug::exec_except([&]() {b.check_range(0, 11); }, out, err);
	liton_sp::debug::exec_except([&]() {b.check_range(0, 10); }, out, err);

	liton_sp::debug::exec_except([&]() {b.check_range(0, LO::half, 2, FL::N.offset); }, out, err);
	liton_sp::debug::exec_except([&]() {b.check_range(0, LO::half, 2, FL::P.offset); }, out, err);
	liton_sp::debug::exec_except([&]() {b.check_range(0, LO::half, 11, FL::N.offset); }, out, err);
	liton_sp::debug::exec_except([&]() {b.check_range(0, LO::half, 11, FL::P.offset); }, out, err);
	out << endl;

	out.close();
	err.close();

	return 0;
}