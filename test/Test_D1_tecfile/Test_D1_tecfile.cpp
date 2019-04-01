#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
using namespace std;
#include "lion_snippets.hpp"
#include "ordered_tec.h"
#define PD_OT
#include "PointData.hpp"

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
		z.read_plt(".", tecfile.last_log, 0, "x", RA::IN, 0, FL::N);
		z.read_plt(".", tecfile.last_log, 0, "y1", RA::IN, 1, FL::N);
		z.read_plt(".", tecfile.last_log, 0, "y2", RA::IN, 2, FL::N);
		out << z.disp_data() << endl;

		liton_sp::debug::exec_except([&]() {z.read_plt(".", tecfile.last_log, 0, "y2", RA::IN, 3, FL::N);}, out, out);
		liton_sp::debug::exec_except([&]() {z.read_plt(".", tecfile.last_log, 0, "xxx", RA::IN, 0, FL::N);}, out, out);
		liton_sp::debug::exec_except([&]() {z.read_plt("aaa", tecfile.last_log, 0, "y2", RA::IN, 0, FL::N);}, out, out);
		z.realloc(0, 101, 0);
		liton_sp::debug::exec_except([&]() {z.read_plt(".", tecfile.last_log, 0, "y2", RA::IN, 2, FL::N);}, out, out);
		z.realloc(0, 99, 0);
		liton_sp::debug::exec_except([&]() {z.read_plt(".", tecfile.last_log, 0, "y2", RA::IN, 2, FL::N);}, out, out);
		z.realloc(10, 100, 10);
		D1::PD_For_N_1D(0, 3, z.size().range(RA::ALL), [&]PD_F_n_i(n, i) {z(n, i) = 0;});
		liton_sp::debug::exec_except([&]() {z.read_plt(".", tecfile.last_log, 0, "y2", RA::N, 2, FL::P);}, out, out);
		out << z.disp_data() << endl;

		out.close();
	}
	catch (const std::exception &err)
	{
		cout << err.what() << endl;
	}
	return 0;
}