#include <iostream>
#include <fstream>
#include <stdexcept>
#include <cmath>
#include <ctime>
using namespace std;
#include "../../scr/liton_cpp_snippets/lion_snippets.hpp"
#include "../../scr/liton_point_data/PointData.hpp"

using namespace liton_pd;

int main(int argc, char** argv)
{
	try
	{
		string name("test");
		cout << name << endl;
		ofstream out((name + "_out.txt").c_str());

		liton_sp::env::disp_env(out);
		out << endl;

		unsigned N = 10000001;
		double x0 = 0;
		double x1 = 1;
		double h = (x1 - x0) / (N - 1);

		D1::PointData<double, 1, LO::center> x(0, N, 0);
		D1::PointData<double, 2, LO::center> f(0, N, 0);
		D1::PD_For_1D(x.size().range(RA::ALL), [&x, &f, &h]PD_F_i(i)
		{
			x(0, i) = static_cast<double>(i) * h;
			f(0, i) = x(0, i);
			f(1, i) = pow(x(0, i), 0.5);
		});

		out << "N = " << N << endl;

		clock_t clock_begin = clock();
		double time = liton_sp::debug::exec_time(100, [&]()
		{
			double sum[2] = { 0, 0 };
			PD_For_N(0, f.N, [&]PD_F_n(n)
			{
				D1::PD_Reduce_1D(f.size().range(RA::ALL),
				sum[n],
				PD_RF_SUM(double),
				[&]PD_F_i(i)->double { return f(n, i); });
				sum[n] /= (N - 1);
				out << "ans = " << sum[n] << endl;
			});
		});
		out << "time used averaged: " << time << endl;

		out.close();

	}
	catch(const exception &err)
	{
		cout << err.what() << endl;
	}
	return 0;
}