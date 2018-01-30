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

	D2::RangeT a;
	out << a.disp() << endl;
	D2::RangeT b(2, 10, 0, 6);
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
	out << endl;
	liton_sp::debug::exec_except([&]() {out << b.size(2) << endl; }, out, err);

	out.close();
	err.close();

	return 0;
}