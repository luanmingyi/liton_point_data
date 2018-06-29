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

	D3::RangeT a;
	out << a.disp() << endl;
	D3::RangeT b(2, 6, 10, 1, 0, 20);
	out << b.disp() << endl;
	out << endl;

	liton_sp::debug::exec_except([&]() {out << a.begin(0) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << a.end(0) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << a.size(0) << endl; }, out, err);
	out << endl;
	liton_sp::debug::exec_except([&]() {out << b.begin(0) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << b.end(0) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << b.size(0) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << b.begin(1) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << b.end(1) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << b.size(1) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << b.begin(2) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << b.end(2) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << b.size(2) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << b.last(2) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << b.bound(2, FL::N) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << b.bound(2, FL::P) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << b.bound(3, FL::P) << endl; }, out, err);
	out << endl;
	liton_sp::debug::exec_except([&]() {out << b.size(3) << endl; }, out, err);

	out << b.disp() << endl;
	out << b.cut_head(-1, -3, 1).disp() << endl;
	out << b.cut_tail(-1, -3, -2).disp() << endl;
	D3::RangeT c = b;
	out << c.tran(2, 1, 0).disp() << endl;
	out << b.overlap(c).disp() << endl;
	liton_sp::debug::exec_except([&]() {out << b.overlap(c.tran(0, 0, 21)).disp() << endl;}, out, err);
	liton_sp::debug::exec_except([&]() {out << b.overlap(c.tran(0, 0, -42)).disp() << endl;}, out, err);
	liton_sp::debug::exec_except([&]() {out << b.overlap(c.tran(0, 0, 1)).disp() << endl;}, out, err);
	out << endl;

	liton_sp::debug::exec_except([&]() {b.check_range(2, LO::half, 0, FL::_N::offset); }, out, err);
	liton_sp::debug::exec_except([&]() {b.check_range(2, LO::half, 0, FL::_P::offset); }, out, err);
	liton_sp::debug::exec_except([&]() {b.check_range(2, LO::half, 22, FL::_N::offset); }, out, err);
	liton_sp::debug::exec_except([&]() {b.check_range(2, LO::half, 22, FL::_P::offset); }, out, err);

	liton_sp::debug::exec_except([&]() {b.check_range(2, 0); }, out, err);
	liton_sp::debug::exec_except([&]() {b.check_range(2, 22); }, out, err);
	liton_sp::debug::exec_except([&]() {b.check_range(2, 21); }, out, err);

	out.close();
	err.close();

	return 0;
}