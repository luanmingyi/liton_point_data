#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
using namespace std;
#include "../../scr/liton_cpp_snippets/lion_snippets.hpp"
#define _CRT_SECURE_NO_WARNINGS
#include "../../scr/liton_ordered_tec/ordered_tec.h"

#ifdef _DEBUG
#define _CHECK_POINTDATA_RANGE
#endif
#include "../../scr/liton_point_data/PointData.hpp"

using namespace liton_pd;
using namespace liton_ot;

int main(int argc, char** argv)
{
	try
	{
		string name(__FILE__);
		name.erase(name.find_last_of('.'));
		cout << name << endl;
		ofstream out((name + "_out.txt").c_str());

		liton_sp::env::disp_env(out);
		out << endl;

		D2::PointData<double, 2, LO::center, LO::center> x(3, 100, 3, 2, 50, 2);
		D2::PointData<double, 2, LO::center, LO::center> u(3, 100, 3, 2, 50, 2);

		TEC_FILE tecfile;
		tecfile.Variables.push_back("x");
		tecfile.Variables.push_back("y");
		tecfile.Variables.push_back("u1");
		tecfile.Variables.push_back("u2");
		tecfile.Zones.push_back(TEC_ZONE());
		for (unsigned d = 0; d != D2::DIM::D; ++d)
		{
			tecfile.Zones[0].Max_C(D2::DIM::D, d) = x.size().size(d, RA::ALL);
			tecfile.Zones[0].Begin_C(D2::DIM::D, d) = x.size().size(d, RA::N);
			tecfile.Zones[0].End_C(D2::DIM::D, d) = x.size().size(d, RA::P);
		}
		tecfile.Zones[0].Data.push_back(TEC_DATA(x.data_pt(0)));
		tecfile.Zones[0].Data.push_back(TEC_DATA(x.data_pt(1)));
		tecfile.Zones[0].Data.push_back(TEC_DATA(u.data_pt(0)));
		tecfile.Zones[0].Data.push_back(TEC_DATA(u.data_pt(1)));

		D2::PD_For_2D(x.size().range(RA::ALL, RA::ALL), [&]PD_F_ij(i, j) {
			x(0, i, j, FL::C, FL::C) = static_cast<double>(i) / (x.size().size(0, RA::IN) - 1);
			x(1, i, j, FL::C, FL::C) = static_cast<double>(j) / (x.size().size(1, RA::IN) - 1);
			u(0, i, j, FL::C, FL::C) = x(0, i, j, FL::C, FL::C) + x(1, i, j, FL::C, FL::C);
			u(1, i, j, FL::C, FL::C) = x(0, i, j, FL::C, FL::C)*x(1, i, j, FL::C, FL::C);
		});

		tecfile.set_echo_mode("full", "full");
		tecfile.write_plt();
		tecfile.last_log.write_echo(out);

		out.close();
	}
	catch (const std::exception& err)
	{
		cout << err.what() << endl;
	}
	return 0;
}