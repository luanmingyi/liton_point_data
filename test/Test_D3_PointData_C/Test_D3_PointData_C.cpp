#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
using namespace std;
#include "../dep/liton_cpp_snippets/lion_snippets.hpp"

#ifdef _DEBUG
#define _CHECK_POINTDATA_RANGE
#endif
#include "../../scr/liton_point_data/PointData.hpp"

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

	out << "line: " << __LINE__ << endl;
	D3::PointData<float, 1, LO::center, LO::center, LO::center> x1;
	out << x1.disp() << endl;
	out << x1.disp_data() << endl;
	D3::PointData<double, 2, LO::center, LO::center, LO::center> x2;
	out << x2.disp() << endl;
	D3::PointData<double, 2, LO::center, LO::center, LO::center> x3(1, 2, 1, 1, 5, 2, 1, 7, 3);
	out << x3.disp() << endl;
	//D2::PointData_C<double, 0, LO::center, LO::center, LO::center> x4;
	//D2::PointData_C<double, -1, LO::center, LO::center, LO::center> x5;
	//D2::PointData_C<float, LO::center, LO::center, LO::center> x6(x1);
	//x3 = x2;
	out << endl;

	out << "line: " << __LINE__ << endl;
	liton_sp::debug::exec_except([&]() {out << x1(0, 0, -1, 0, FL::C, FL::C, FL::C) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {x1.alloc(0, 0, 0, 0, 0, 0, 0, 0, 0); }, out, err);
	liton_sp::debug::exec_except([&]() {x1.alloc(0, 0, 0, 0, 0, 0, 1, 0, 0); }, out, err);
#ifndef __linux__
	const unsigned big = 0.1 * 1024 / sizeof(float) / x1.N * 1024 * 1024;
	out << big* big << endl;
	liton_sp::debug::exec_except([&]() {x1.alloc(0, big, 0, 0, big, 0, 1, big, 0); }, out, err);
#endif
	liton_sp::debug::exec_except([&]() {x1.alloc(1, 2, 1, 0, 7, 0, 1, 5, 0); }, out, err);
	liton_sp::debug::exec_except([&]() {x1.alloc(1, 2, 1, 0, 7, 0, 1, 5, 0); }, out, err);
	out << endl;

	liton_sp::debug::exec_except([&]() {x1.realloc(1, 2, 1, 1, 5, 2, 1, 7, 3); }, out, err);
	liton_sp::debug::exec_except([&]() {x2.realloc(1, 2, 1, 1, 5, 2, 1, 7, 3); }, out, err);
	out << endl;

	liton_sp::debug::exec_except([&]() {x2.clear(); }, out, err);
	liton_sp::debug::exec_except([&]() {x3.clear(); x3.alloc(1, 2, 1, 1, 5, 2, 1, 7, 3); }, out, err);
	liton_sp::debug::exec_except([&]() {x2.alloc(1, 2, 1, 1, 5, 2, 1, 7, 3); x3.alloc(1, 2, 1, 1, 5, 2, 1, 7, 3); }, out, err);
	out << endl;

	out << "line: " << __LINE__ << endl;
	out << x3.size().disp() << endl;
	out << endl;
	
	out << "line: " << __LINE__ << endl;
	
	D3::PD_For_3D(x1.size().range(RA::N, RA::IN, RA::IN), [&x1]PD_F_ijk(i, j, k) { x1(0, i, j, k, FL::C, FL::C, FL::C) = 122; });
	D3::PD_For_3D(x1.size().range(RA::N, RA::N, RA::N), [&x1]PD_F_ijk(i, j, k) { x1(0, i, j, k, FL::C, FL::C, FL::C) = 111; });
	D3::PD_For_3D(x1.size().range(RA::N, RA::N, RA::P), [&x1]PD_F_ijk(i, j, k) { x1(0, i, j, k, FL::C, FL::C, FL::C) = 113; });
	D3::PD_For_3D(x1.size().range(RA::N, RA::P, RA::N), [&x1]PD_F_ijk(i, j, k) { x1(0, i, j, k, FL::C, FL::C, FL::C) = 131; });
	D3::PD_For_3D(x1.size().range(RA::N, RA::P, RA::P), [&x1]PD_F_ijk(i, j, k) { x1(0, i, j, k, FL::C, FL::C, FL::C) = 133; });
	D3::PD_For_3D(x1.size().range(RA::N, RA::N, RA::IN), [&x1]PD_F_ijk(i, j, k) { x1(0, i, j, k, FL::C, FL::C, FL::C) = 112; });
	D3::PD_For_3D(x1.size().range(RA::N, RA::P, RA::IN), [&x1]PD_F_ijk(i, j, k) { x1(0, i, j, k, FL::C, FL::C, FL::C) = 132; });
	D3::PD_For_3D(x1.size().range(RA::N, RA::IN, RA::N), [&x1]PD_F_ijk(i, j, k) { x1(0, i, j, k, FL::C, FL::C, FL::C) = 121; });
	D3::PD_For_3D(x1.size().range(RA::N, RA::IN, RA::P), [&x1]PD_F_ijk(i, j, k) { x1(0, i, j, k, FL::C, FL::C, FL::C) = 123; });

	D3::PD_For_3D(x1.size().range(RA::IN, RA::IN, RA::IN), [&x1]PD_F_ijk(i, j, k) { x1(0, i, j, k, FL::C, FL::C, FL::C) = i * 100 + j * 10 + k; });
	D3::PD_For_3D(x1.size().range(RA::IN, RA::N, RA::N), [&x1]PD_F_ijk(i, j, k) { x1(0, i, j, k, FL::C, FL::C, FL::C) = 211; });
	D3::PD_For_3D(x1.size().range(RA::IN, RA::N, RA::P), [&x1]PD_F_ijk(i, j, k) { x1(0, i, j, k, FL::C, FL::C, FL::C) = 213; });
	D3::PD_For_3D(x1.size().range(RA::IN, RA::P, RA::N), [&x1]PD_F_ijk(i, j, k) { x1(0, i, j, k, FL::C, FL::C, FL::C) = 231; });
	D3::PD_For_3D(x1.size().range(RA::IN, RA::P, RA::P), [&x1]PD_F_ijk(i, j, k) { x1(0, i, j, k, FL::C, FL::C, FL::C) = 233; });
	D3::PD_For_3D(x1.size().range(RA::IN, RA::N, RA::IN), [&x1]PD_F_ijk(i, j, k) { x1(0, i, j, k, FL::C, FL::C, FL::C) = 212; });
	D3::PD_For_3D(x1.size().range(RA::IN, RA::P, RA::IN), [&x1]PD_F_ijk(i, j, k) { x1(0, i, j, k, FL::C, FL::C, FL::C) = 232; });
	D3::PD_For_3D(x1.size().range(RA::IN, RA::IN, RA::N), [&x1]PD_F_ijk(i, j, k) { x1(0, i, j, k, FL::C, FL::C, FL::C) = 221; });
	D3::PD_For_3D(x1.size().range(RA::IN, RA::IN, RA::P), [&x1]PD_F_ijk(i, j, k) { x1(0, i, j, k, FL::C, FL::C, FL::C) = 223; });

	D3::PD_For_3D(x1.size().range(RA::P, RA::IN, RA::IN), [&x1]PD_F_ijk(i, j, k) { x1(0, i, j, k, FL::C, FL::C, FL::C) = 322; });
	D3::PD_For_3D(x1.size().range(RA::P, RA::N, RA::N), [&x1]PD_F_ijk(i, j, k) { x1(0, i, j, k, FL::C, FL::C, FL::C) = 311; });
	D3::PD_For_3D(x1.size().range(RA::P, RA::N, RA::P), [&x1]PD_F_ijk(i, j, k) { x1(0, i, j, k, FL::C, FL::C, FL::C) = 313; });
	D3::PD_For_3D(x1.size().range(RA::P, RA::P, RA::N), [&x1]PD_F_ijk(i, j, k) { x1(0, i, j, k, FL::C, FL::C, FL::C) = 331; });
	D3::PD_For_3D(x1.size().range(RA::P, RA::P, RA::P), [&x1]PD_F_ijk(i, j, k) { x1(0, i, j, k, FL::C, FL::C, FL::C) = 333; });
	D3::PD_For_3D(x1.size().range(RA::P, RA::N, RA::IN), [&x1]PD_F_ijk(i, j, k) { x1(0, i, j, k, FL::C, FL::C, FL::C) = 312; });
	D3::PD_For_3D(x1.size().range(RA::P, RA::P, RA::IN), [&x1]PD_F_ijk(i, j, k) { x1(0, i, j, k, FL::C, FL::C, FL::C) = 332; });
	D3::PD_For_3D(x1.size().range(RA::P, RA::IN, RA::N), [&x1]PD_F_ijk(i, j, k) { x1(0, i, j, k, FL::C, FL::C, FL::C) = 321; });
	D3::PD_For_3D(x1.size().range(RA::P, RA::IN, RA::P), [&x1]PD_F_ijk(i, j, k) { x1(0, i, j, k, FL::C, FL::C, FL::C) = 323; });
	out << x1.disp_data() << endl;
	out << endl;
	
	out << "line: " << __LINE__ << endl;
	D3::PD_For_3D(x1.size().range(RA::IN, RA::IN, RA::N), [&x1]PD_F_ijk(i, j, k) { x1(0, i, j, k, FL::C, FL::C, FL::C) = -x1(0, i, j, x1.size().mirror(2, FL::N, k), FL::C, FL::C, FL::C); });
	D3::PD_For_3D(x1.size().range(RA::IN, RA::IN, RA::P), [&x1]PD_F_ijk(i, j, k) { x1(0, i, j, k, FL::C, FL::C, FL::C) = -x1(0, i, j, x1.size().mirror(2, FL::P, k), FL::C, FL::C, FL::C); });
	D3::PD_For_3D(x1.size().range(RA::IN, RA::N, RA::IN), [&x1]PD_F_ijk(i, j, k) { x1(0, i, j, k, FL::C, FL::C, FL::C) = -x1(0, i, x1.size().mirror(1, FL::N, j), k, FL::C, FL::C, FL::C); });
	D3::PD_For_3D(x1.size().range(RA::IN, RA::P, RA::IN), [&x1]PD_F_ijk(i, j, k) { x1(0, i, j, k, FL::C, FL::C, FL::C) = -x1(0, i, x1.size().mirror(1, FL::P, j), k, FL::C, FL::C, FL::C); });
	D3::PD_For_3D(x1.size().range(RA::N, RA::IN, RA::IN), [&x1]PD_F_ijk(i, j, k) { x1(0, i, j, k, FL::C, FL::C, FL::C) = -x1(0, x1.size().mirror(0, FL::N, i), j, k, FL::C, FL::C, FL::C); });
	D3::PD_For_3D(x1.size().range(RA::P, RA::IN, RA::IN), [&x1]PD_F_ijk(i, j, k) { x1(0, i, j, k, FL::C, FL::C, FL::C) = -x1(0, x1.size().mirror(0, FL::P, i), j, k, FL::C, FL::C, FL::C); });
	out << x1.disp_data() << endl;
	out << endl;

	out << "line: " << __LINE__ << endl;
	liton_sp::debug::exec_except([&]() {out << x1(0, 0, -1, 0, FL::C, FL::C, FL::C) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x1(0, 0, -2, 0, FL::C, FL::C, FL::C) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x1(0, 0, 7, 0, FL::C, FL::C, FL::C) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x1(0, 0, 0, -2, FL::C, FL::C, FL::C) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x1(0, 0, 0, 13, FL::C, FL::C, FL::C) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x1(0, 3, 7, 13, FL::C, FL::C, FL::C) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x1(1, 0, 0, 0, FL::C, FL::C, FL::C) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x1(0, 0, -1, 0, FL::C, FL::N, FL::C) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x1(0, 0, -1, 0, FL::C, FL::C, FL::P) << endl; }, out, err);
	out << endl;

	out << "line: " << __LINE__ << endl;
	D3::PD_For_N_3D(0, x2.N, x2.size().range(RA::ALL, RA::ALL, RA::ALL), [&x2]PD_F_n_ijk(n, i, j, k) {
		x2(n, i, j, k, FL::C, FL::C, FL::C) = static_cast<int>(n) * 1000 + i * 100 + j * 10 + k;
	});
	out << x2.disp_data() << endl;
	out << endl;

	out << "line: " << __LINE__ << endl;
	double max = 0;
	D3::PD_Reduce_3D(x1.size().range(RA::IN, RA::IN, RA::IN), max,
		[]PD_RF(double, x, xx) { xx = x > xx ? x : xx; },
		[&x1]PD_F_ijk(i, j, k)->double { return x1(0, i, j, k, FL::C, FL::C, FL::C); }
	);
	out << "max = " << max << endl;
	max = 0;
	D3::PD_Reduce_N_3D(0, x2.N, x2.size().range(RA::IN, RA::IN, RA::IN), max,
		[]PD_RF(double, x, xx) { xx = x > xx ? x : xx; },
		[&x2]PD_F_n_ijk(n, i, j, k)->double { return x2(n, i, j, k, FL::C, FL::C, FL::C); }
	);
	out << "max = " << max << endl;
	out << endl;

	out << "line: " << __LINE__ << endl;
	D3::PD_For_2D(1, 2, x1.size().range(RA::ALL, RA::ALL, RA::ALL), [&x1]PD_F_ij(ii, jj) { x1(0, 0, ii, jj, FL::C, FL::C, FL::C) = x1(0, 1, ii, jj, FL::C, FL::C, FL::C); });
	out << x1.disp_data() << endl;
	max = 0;
	D3::PD_Reduce_2D(1, 2, x1.size().range(RA::ALL, RA::ALL, RA::ALL), max,
		[]PD_RF(double, x, xx) { xx = x > xx ? x : xx; },
		[&x1]PD_F_ij(ii, jj)->double { return x1(0, 0, ii, jj, FL::C, FL::C, FL::C); }
	);
	out << "max = " << max << endl;
	out << endl;

	D3::PD_For_N_2D(0, x2.N, 1, 2, x2.size().range(RA::ALL, RA::ALL, RA::ALL), [&x2]PD_F_n_ij(n, ii, jj) {
		x2(n, 0, ii, jj, FL::C, FL::C, FL::C) = x2(n, 1, ii, jj, FL::C, FL::C, FL::C);
	});
	out << x2.disp_data() << endl;
	max = 0;
	D3::PD_Reduce_N_2D(0, x2.N, 1, 2, x2.size().range(RA::ALL, RA::ALL, RA::ALL), max,
		[]PD_RF(double, x, xx) { xx = x > xx ? x : xx; },
		[&x2]PD_F_n_ij(n, ii, jj)->double { return x2(n, 0, ii, jj, FL::C, FL::C, FL::C); }
	);
	out << "max = " << max << endl;
	out << endl;

	out << "line: " << __LINE__ << endl;
	D3::PD_For_1D(1, x1.size().range(RA::ALL, RA::ALL, RA::ALL), [&x1]PD_F_i(ii) { x1(0, 0, ii, 0, FL::C, FL::C, FL::C) = x1(0, 0, ii, 1, FL::C, FL::C, FL::C); });
	out << x1.disp_data() << endl;
	max = 0;
	D3::PD_Reduce_1D(1, x1.size().range(RA::ALL, RA::ALL, RA::ALL), max,
		[]PD_RF(double, x, xx) { xx = x > xx ? x : xx; },
		[&x1]PD_F_i(ii)->double { return x1(0, 0, ii, 0, FL::C, FL::C, FL::C); }
	);
	out << "max = " << max << endl;
	out << endl;

	D3::PD_For_N_1D(0, x2.N, 1, x2.size().range(RA::ALL, RA::ALL, RA::ALL), [&x2]PD_F_n_i(n, ii) {
		x2(n, 0, ii, 0, FL::C, FL::C, FL::C) = x2(n, 0, ii, 1, FL::C, FL::C, FL::C);
	});
	out << x2.disp_data() << endl;
	max = 0;
	D3::PD_Reduce_N_1D(0, x2.N, 1, x2.size().range(RA::ALL, RA::ALL, RA::ALL), max,
		[]PD_RF(double, x, xx) { xx = x > xx ? x : xx; },
		[&x2]PD_F_n_i(n, ii)->double { return x2(n, 0, ii, 0, FL::C, FL::C, FL::C); }
	);
	out << "max = " << max << endl;
	out << endl;

	out.close();
	err.close();

	return 0;
}