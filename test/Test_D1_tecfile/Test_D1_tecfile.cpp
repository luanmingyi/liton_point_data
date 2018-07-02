#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
using namespace std;
#include "../../scr/liton_cpp_snippets/lion_snippets.hpp"
#include "../../scr/liton_ordered_tec/ordered_tec.h"
#define PD_OT
#include "../../scr/liton_point_data/PointData.hpp"

using namespace liton_pd;
using namespace liton_ot;

int main(int argc, char** argv)
{
	try
	{
		string name("test");
		cout << name << endl;
		ofstream out((name + "_out.txt").c_str());

		liton_sp::env::disp_env(out);
		out << endl;

		D1::PointData<double, 1, LO::center> x(3, 100, 3);
		D1::PointData<double, 2, LO::center> y(3, 100, 3);

		TEC_FILE tecfile;
		tecfile.Variables.push_back("x");
		tecfile.Variables.push_back("y1");
		tecfile.Variables.push_back("y2");
		tecfile.Zones.push_back(TEC_ZONE());
		tecfile.Zones[0].Max[0] = x.size().size(0, RA::ALL);
		tecfile.Zones[0].Begin[0] = x.size().size(0, RA::N);
		tecfile.Zones[0].End[0] = x.size().size(0, RA::P);
		tecfile.Zones[0].Data.push_back(TEC_DATA(x.data_pt(0)));
		tecfile.Zones[0].Data.push_back(TEC_DATA(y.data_pt(0)));
		tecfile.Zones[0].Data.push_back(TEC_DATA(y.data_pt(1)));

		D1::PD_For_1D(x.size().range(RA::ALL), [&]PD_F_i(i)
		{
			x(0, i, FL::C) = static_cast<double>(i) / (x.size().size(0, RA::IN) - 1);
			y(0, i, FL::C) = x(0, i, FL::C) / 2;
			y(1, i, FL::C) = x(0, i, FL::C) * x(0, i, FL::C);
		});

		tecfile.write_plt();
		tecfile.last_log.write_echo(out);
		tecfile.last_log.write_xml();

		D1::PointData<double, 3, LO::center> z(0, 100, 0);
		z.read_plt(".", tecfile.last_log, 0, "x", 0);
		z.read_plt(".", tecfile.last_log, 0, "y1", 1);
		z.read_plt(".", tecfile.last_log, 0, "y2", 2);
		out << z.disp_data() << endl;

		liton_sp::debug::exec_except([&]() {z.read_plt(".", tecfile.last_log, 0, "y2", 3);}, out, out);
		liton_sp::debug::exec_except([&]() {z.read_plt(".", tecfile.last_log, 0, "xxx", 0);}, out, out);
		liton_sp::debug::exec_except([&]() {z.read_plt("aaa", tecfile.last_log, 0, "y2", 0);}, out, out);
		z.realloc(0, 101, 0);
		liton_sp::debug::exec_except([&]() {z.read_plt(".", tecfile.last_log, 0, "y2", 2);}, out, out);
		z.realloc(0, 99, 0);
		liton_sp::debug::exec_except([&]() {z.read_plt(".", tecfile.last_log, 0, "y2", 2);}, out, out);
		z.realloc(1, 100, 1);
		liton_sp::debug::exec_except([&]() {z.read_plt(".", tecfile.last_log, 0, "y2", 2);}, out, out);

		out.close();
	}
	catch (const std::exception &err)
	{
		cout << err.what() << endl;
	}
	return 0;
}