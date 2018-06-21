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
	string name("test");
	cout << name << endl;
	ofstream out((name + "_out.txt").c_str());
	ofstream err((name + "_err.txt").c_str());

	liton_sp::env::disp_env(out);
	out << endl;

	D1::SizeT s1;
	out << s1.disp() << endl;
	D1::SizeT s2(2, 10, 3);
	out << s2.disp() << endl;
	out << endl;

	liton_sp::debug::exec_except([&]() {out << s1.in(0) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << s1.n(0) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << s1.p(0) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << s2.in(0) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << s2.n(0) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << s2.p(0) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << s2.in(1) << endl; }, out, err);
	out << endl;

	liton_sp::debug::exec_except([&]() {out << s2.begin(0, RA::ALL) << " " << s2.end(0, RA::ALL) << " " << s2.size(0, RA::ALL) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << s2.begin(0, RA::IN) << " " << s2.end(0, RA::IN) << " " << s2.size(0, RA::IN) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << s2.begin(0, RA::N) << " " << s2.end(0, RA::N) << " " << s2.size(0, RA::N) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << s2.begin(0, RA::P) << " " << s2.end(0, RA::P) << " " << s2.size(0, RA::P) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << s2.begin(1, RA::P) << " " << s2.end(1, RA::P) << " " << s2.size(1, RA::P) << endl; }, out, err);
	out << endl;

	out << s2.range(RA::ALL).disp() << endl;
	out << s2.range(RA::IN).disp() << endl;
	out << s2.range(RA::N).disp() << endl;
	out << s2.range(RA::P).disp() << endl;
	out << endl;

	out << s2.range(RA::ALL, true, false).disp() << endl;
	out << s2.range(RA::IN, true, false).disp() << endl;
	out << s2.range(RA::N, true, false).disp() << endl;
	out << s2.range(RA::P, true, false).disp() << endl;
	out << endl;

	out << s2.range(RA::ALL, false, true).disp() << endl;
	out << s2.range(RA::IN, false, true).disp() << endl;
	out << s2.range(RA::N, false, true).disp() << endl;
	out << s2.range(RA::P, false, true).disp() << endl;
	out << endl;

	liton_sp::debug::exec_except([&]() {s1.check(); }, out, err);
	liton_sp::debug::exec_except([&]() {s2.check(); }, out, err);
	liton_sp::debug::exec_except([&]() {D1::SizeT(1, 0, 0).check(); }, out, err);
	out << endl;

	liton_sp::debug::exec_except([&]() {s2.check_range(0, 0); }, out, err);
	liton_sp::debug::exec_except([&]() {s2.check_range(0, -3); }, out, err);
	liton_sp::debug::exec_except([&]() {s2.check_range(0, 13); }, out, err);

	liton_sp::debug::exec_except([&]() {s2.check_range(0, LO::half, -3, FL::N.offset); }, out, err);
	liton_sp::debug::exec_except([&]() {s2.check_range(0, LO::half, -3, FL::P.offset); }, out, err);
	liton_sp::debug::exec_except([&]() {s2.check_range(0, LO::half, 13, FL::N.offset); }, out, err);
	liton_sp::debug::exec_except([&]() {s2.check_range(0, LO::half, 13, FL::P.offset); }, out, err);
	out << endl;

	out << s2.disp() << endl;
	out << s2.mirror(0, FL::N, -2) << endl;
	out << s2.mirror(0, FL::P, 11) << endl;
	out << s2.periodic(0, FL::N, -2) << endl;
	out << s2.periodic(0, FL::P, 11) << endl;

	out.close();
	err.close();

	return 0;
}