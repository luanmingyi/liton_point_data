#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
using namespace std;
#include "../dep/liton_cpp_snippets/lion_snippets.hpp"

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

	liton_sp::debug::exec_except([&]() {s1.check(); }, out, err);
	liton_sp::debug::exec_except([&]() {s2.check(); }, out, err);
	liton_sp::debug::exec_except([&]() {D1::SizeT(1, 0, 0).check(); }, out, err);
	out << endl;

	liton_sp::debug::exec_except([&]() {s2.check_range(0); }, out, err);
	liton_sp::debug::exec_except([&]() {s2.check_range(-3); }, out, err);
	liton_sp::debug::exec_except([&]() {s2.check_range(13); }, out, err);

	out.close();
	err.close();

	return 0;
}