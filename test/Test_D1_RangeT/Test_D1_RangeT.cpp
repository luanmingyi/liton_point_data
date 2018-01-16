#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
using namespace std;
#include "../dep/lion_snippets.hpp"

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

	D1::RangeT a;
	out << a.disp() << endl;
	D1::RangeT b(2, 10);
	out << b.disp() << endl;
	out << endl;

	exec([&]() {out << a.begin(0) << endl; }, out, err);
	exec([&]() {out << a.end(0) << endl; }, out, err);
	exec([&]() {out << a.size(0) << endl; }, out, err);
	out << endl;
	exec([&]() {out << b.begin(0) << endl; }, out, err);
	exec([&]() {out << b.end(0) << endl; }, out, err);
	exec([&]() {out << b.size(0) << endl; }, out, err);
	out << endl;
	exec([&]() {out << b.size(1) << endl; }, out, err);

	out.close();
	err.close();

	return 0;
}