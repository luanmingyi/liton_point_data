#include <iostream>
#include <fstream>
#include <stdexcept>
#include <cmath>
#include <ctime>
using namespace std;

#include "../../scr/liton_cpp_snippets/lion_snippets.hpp"
#include "../../scr/liton_ordered_tec/ordered_tec.h"
using namespace liton_ot;

#ifdef _DEBUG
#define _CHECK_POINTDATA_RANGE
#endif
#include "../../scr/liton_point_data/PointData.hpp"
using namespace liton_pd;

int main(int argc, char** argv)
{
	try
	{
		string name("test");
		cout << name << endl;
		ofstream out((name + "_out.txt").c_str(), ios::app);
		liton_sp::env::disp_env(out);
		out << endl;

		const double pi = 3.1415926;
		const unsigned N = 500;

		double ax = 4;
		unsigned Nx = N + 1;
		double x0 = -pi;
		double x1 = pi;
		double ay = 1;
		unsigned Ny = 2 * N + 1;
		double y0 = -2 * pi;
		double y1 = 2 * pi;

		double hx = (x1 - x0) / (Nx - 1);
		double a_hh_x = ax / hx / hx;
		double hy = (y1 - y0) / (Ny - 1);
		double a_hh_y = ay / hy / hy;

		D2::PointData<double, 2, LO::center, LO::center> x(1, Nx, 1, 1, Ny, 1);
		D2::PointData<double, 2, LO::center, LO::center> u(1, Nx, 1, 1, Ny, 1);
		D2::PointData<double, 1, LO::center, LO::center> dudt(0, Nx, 0, 0, Ny, 0);

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
			tecfile.Zones[0].Data.push_back(TEC_DATA(x.data_pt(d)));
		}
		tecfile.Zones[0].Data.push_back(TEC_DATA(u.data_pt(0)));
		tecfile.Zones[0].Data.push_back(TEC_DATA(u.data_pt(1)));

		D2::PD_For_2D(x.size().range(RA::IN, RA::IN), [&x, hx, x0, hy, y0]PD_F_ij(i, j) {
			x(0, i, j) = x0 + static_cast<double>(i)*hx;
			x(1, i, j) = y0 + static_cast<double>(j)*hy;
		});

		double dt = 0.1 / (a_hh_x > a_hh_y ? a_hh_x : a_hh_y);
		double time_stop = 0.1;
		double time_print = time_stop / 20;

		unsigned step;
		double flowtime;
		double time_next_print;

		double usingtime = liton_sp::debug::exec_time(1, [&]() {
			step = 0;
			flowtime = 0;
			time_next_print = 0;
			D2::PD_For_2D(u.size().range(RA::IN,RA::IN), [&x, &u]PD_F_ij(i, j) {
				u(0, i, j) = sin(x(0, i, j))*cos(x(1, i, j) / 2.0);
			});

			for (;;)
			{
				D2::PD_For_2D(u.size().range(RA::IN, RA::N), [&u]PD_F_ij(i, j) { u(0, i, j) = u(0, i, u.size().mirror(1, FL::N, j)); });
				D2::PD_For_2D(u.size().range(RA::IN, RA::P), [&u]PD_F_ij(i, j) { u(0, i, j) = u(0, i, u.size().mirror(1, FL::P, j)); });
				D2::PD_For_2D(u.size().range(RA::N, RA::IN), [&u]PD_F_ij(i, j) { u(0, i, j) = -u(0, u.size().mirror(0, FL::N, i), j); });
				D2::PD_For_2D(u.size().range(RA::P, RA::IN), [&u]PD_F_ij(i, j) { u(0, i, j) = -u(0, u.size().mirror(0, FL::P, i), j); });

				D2::PD_For_2D(dudt.size().range(RA::IN, RA::IN), [&dudt, &u, a_hh_x]PD_F_ij(i, j) {
					dudt(0, i, j) = a_hh_x*(u(0, i - 1, j) - 2 * u(0, i, j) + u(0, i + 1, j));
				});
				D2::PD_For_2D(dudt.size().range(RA::IN, RA::IN), [&dudt, &u, a_hh_y]PD_F_ij(i, j) {
					dudt(0, i, j) += a_hh_y*(u(0, i, j - 1) - 2 * u(0, i, j) + u(0, i, j + 1));
				});
				D2::PD_For_2D(u.size().range(RA::IN, RA::IN), [&dudt, &u, dt]PD_F_ij(i, j) {
					u(0, i, j) += dt*dudt(0, i, j);
				});

				++step;
				flowtime += dt;
				if (flowtime >= time_stop)
				{
					break;
				}
				if (flowtime >= time_next_print)
				{
					cout << "step: " << step << " time: " << flowtime << endl;
					time_next_print += time_print;
				}
			}
		});

		D2::PD_For_2D(u.size().range(RA::IN, RA::IN), [&x, &u, ax, ay, time_stop]PD_F_ij(i, j) {
			u(1, i, j) = exp(-(ax + ay / 4)*time_stop)*sin(x(0, i, j))*cos(x(1, i, j) / 2.0);
		});

		tecfile.set_echo_mode("full", "full");
		tecfile.write_plt();
		tecfile.last_log.write_xml();

		out << "time: " << flowtime << endl;
		out << "step: " << step << endl;
		out << "using time: " << usingtime << endl;
		out << endl;
		out.close();
	}
	catch (const exception &err)
	{
		cout << err.what() << endl;
	}
	return 0;
}