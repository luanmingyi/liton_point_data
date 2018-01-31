#include <iostream>
#include <fstream>
#include <stdexcept>
#include <cmath>
#include <ctime>
using namespace std;
#include "../dep/liton_cpp_snippets/lion_snippets.hpp"

#ifdef _DEBUG
#define _CHECK_POINTDATA_RANGE
#endif
#include "../../scr/PointData.hpp"

using namespace liton_pd;

int main(int argc, char** argv)
{
	try
	{
		string name(__FILE__);
		name.erase(name.find_last_of('.'));
		cout << name << endl;
		ofstream out((name + "_out.txt").c_str(), ios::app);
		liton_sp::env::disp_env(out);
		out << endl;

		double a = 2;

		double pi = 3.1415926;
		unsigned N = 5001;
		double x0 = -pi;
		double x1 = pi;
		double h = (x1 - x0) / (N - 1);
		double a_hh = a / h / h;

		D1::PointData<double, 1, LO::center> x(0, N, 0);
		D1::PointData<double, 2, LO::center> u(1, N, 1);
		D1::PointData<double, 1, LO::center> dudt(0, N, 0);

		D1::PD_For_1D(x.size().range(RA::IN), [&x, &h, &x0]PD_F_i(i) { x(0, i, FL::C) = x0 + static_cast<double>(i)*h; });

		double dt = 0.1 / a_hh;
		double tiem_stop = 0.5;
		double time_print = tiem_stop / 20;

		unsigned step;
		double flowtime;
		double time_next_print;		

		double usingtime = liton_sp::debug::exec_time(2, [&]() {
			step = 0;
			flowtime = 0;
			time_next_print = 0;
			D1::PD_For_1D(u.size().range(RA::IN), [&x, &u]PD_F_i(i) { u(0, i, FL::C) = cos(x(0, i, FL::C)); });

			for (;;)
			{
				D1::PD_For_1D(u.size().range(RA::N), [&u]PD_F_i(i) { u(0, i, FL::C) = u(0, u.size().mirror(0, FL::N, i), FL::C); });
				D1::PD_For_1D(u.size().range(RA::P), [&u]PD_F_i(i) { u(0, i, FL::C) = u(0, u.size().mirror(0, FL::P, i), FL::C); });

				D1::PD_For_1D(dudt.size().range(RA::IN), [&dudt, &u, &a_hh]PD_F_i(i) {
					dudt(0, i, FL::C) = a_hh*(u(0, i - 1, FL::C) - 2 * u(0, i, FL::C) + u(0, i + 1, FL::C));
				});
				D1::PD_For_1D(u.size().range(RA::IN), [&dudt, &u, &dt]PD_F_i(i) {
					u(0, i, FL::C) = u(0, i, FL::C) + dt*dudt(0, i, FL::C);
				});

				++step;
				flowtime += dt;
				if (flowtime >= tiem_stop)
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

		D1::PD_For_1D(u.size().range(RA::IN), [&x, &u, &a, &tiem_stop]PD_F_i(i) { u(1, i, FL::C) = exp(-a*tiem_stop)*cos(x(0, i, FL::C)); });

		ofstream out_data((name + "_data.dat").c_str());
		out_data << u.disp_data() << endl;
		out_data.close();

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