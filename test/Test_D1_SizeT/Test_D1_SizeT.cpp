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

	D1::SizeT s1;
	out << s1.disp() << endl;
	D1::SizeT s2(10, 2, 3);
	out << s2.disp() << endl;
	out << endl;

	exec([&]() {out << s1.in(0) << endl; }, out, err);
	exec([&]() {out << s1.n(0) << endl; }, out, err);
	exec([&]() {out << s1.p(0) << endl; }, out, err);
	exec([&]() {out << s2.in(0) << endl; }, out, err);
	exec([&]() {out << s2.n(0) << endl; }, out, err);
	exec([&]() {out << s2.p(0) << endl; }, out, err);
	exec([&]() {out << s2.in(1) << endl; }, out, err);
	out << endl;

	exec([&]() {out << s2.sum(0) << endl; }, out, err);
	exec([&]() {out << s2.sum(1) << endl; }, out, err);
	out << endl;

	exec([&]() {out << s2.begin(0, RA::ALL) << " " << s2.end(0, RA::ALL) << " " << s2.size(0, RA::ALL) << endl; }, out, err);
	exec([&]() {out << s2.begin(0, RA::IN) << " " << s2.end(0, RA::IN) << " " << s2.size(0, RA::IN) << endl; }, out, err);
	exec([&]() {out << s2.begin(0, RA::N) << " " << s2.end(0, RA::N) << " " << s2.size(0, RA::N) << endl; }, out, err);
	exec([&]() {out << s2.begin(0, RA::P) << " " << s2.end(0, RA::P) << " " << s2.size(0, RA::P) << endl; }, out, err);
	exec([&]() {out << s2.begin(1, RA::P) << " " << s2.end(1, RA::P) << " " << s2.size(1, RA::P) << endl; }, out, err);
	out << endl;

	out << s2.range(RA::ALL).disp() << endl;
	out << s2.range(RA::IN).disp() << endl;
	out << s2.range(RA::N).disp() << endl;
	out << s2.range(RA::P).disp() << endl;
	out << endl;

	exec([&]() {s1.check(); }, out, err);
	exec([&]() {s2.check(); }, out, err);
	exec([&]() {D1::SizeT(0, 1, 0).check(); }, out, err);
	out << endl;

	exec([&]() {s2.check_range(0); }, out, err);
	exec([&]() {s2.check_range(-3); }, out, err);
	exec([&]() {s2.check_range(13); }, out, err);

	out.close();
	err.close();

	return 0;
}