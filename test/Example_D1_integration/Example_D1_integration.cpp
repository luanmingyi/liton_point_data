#include <iostream>
#include <fstream>
#include <stdexcept>
#include <cmath>
#include<ctime>
using namespace std;
#include "../dep/lion_snippets.hpp"

#ifdef _DEBUG
	#define _CHECK_POINTDATA_RANGE
#endif
#include "../../scr/PointData.hpp"

using namespace liton;

int main(int argc, char** argv)
{
	try
	{
		ofstream out("example.txt");

		liton_sp::env::disp_env(out);
		out << endl;

		unsigned N = 10000001;
		double x0 = 0;
		double x1 = 1;
		double h = (x1 - x0) / (N - 1);

		D1::PointData_C<double> x(N, 0, 0);
		D1::PointData_C<double, 2> f(N, 0, 0);
		D1::For_PD_1D(x.size().range(RA::ALL), [&x, &f, &h]PD_F_i(i)
		{
			x(i, 0) = static_cast<double>(i) * h;
			f(i, 0) = x(i, 0);
			f(i, 1) = pow(x(i, 0), 0.5);
		});

		out << "N = " << N << endl;

		clock_t clock_begin = clock();
		double sum[2] = { 0, 0 };
		For_N(0, f.N, [&]PD_F_n(n)
		{
			D1::Reduce_PD_1D(f.size().range(RA::ALL),
			[]PD_RF(double, x, xx) { xx += x; },
			[&]PD_F_i(i)->double { return f(i, n); },
			sum[n]);
			sum[n] /= (N - 1);
			out << "ans = " << sum[n] << endl;
		});
		double time = double(clock() - clock_begin) / CLOCKS_PER_SEC;
		out << "time used: " << time << endl;

		out.close();

	}
	catch(const exception &err)
	{
		cout << err.what() << endl;
	}
	return 0;
}