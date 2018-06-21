#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
using namespace std;
#include "../../scr/liton_cpp_snippets/lion_snippets.hpp"

#ifdef _DEBUG
#define _CHECK_POINTDATA_RANGE
#endif
#include "../../scr/liton_point_data/PointData.hpp"

using namespace liton_pd;

int main(int argc, char** argv)
{
	string name("test");
	cout << name << endl;
	ofstream out((name + "_out.txt").c_str());
	ofstream err((name + "_err.txt").c_str());

	liton_sp::env::disp_env(out);
	out << endl;

	D3::PointData<double, 2, LO::center, LO::center, LO::center> x0;
	D3::PointData<double, 2, LO::center, LO::center, LO::half> x1;
	D3::PointData<double, 2, LO::center, LO::half, LO::center> x2;
	D3::PointData<double, 2, LO::center, LO::half, LO::half> x3;
	D3::PointData<double, 2, LO::half, LO::center, LO::center> x4;
	D3::PointData<double, 2, LO::half, LO::center, LO::half> x5;
	D3::PointData<double, 2, LO::half, LO::half, LO::center> x6;
	D3::PointData<double, 2, LO::half, LO::half, LO::half> x7;
	//D2::PointData<double, 0, LO::center, LO::half, LO::half> x4;
	//D2::PointData<double, -1, LO::center, LO::half, LO::half> x5;
	//D2::PointData<double, 1, LO::center, LO::half, LO::half> x6(x3);
	//x3 = x3;
	out << endl;

	liton_sp::debug::exec_except([&]() {x1.alloc(0, 0, 0, 0, 0, 0, 0, 0, 0); }, out, err);
	liton_sp::debug::exec_except([&]() {x1.alloc(0, 0, 0, 0, 0, 0, 1, 0, 0); }, out, err);
#ifndef __linux__
	const unsigned big = 0.1 * 1024 / sizeof(float) / x1.N * 1024 * 1024;
	out << big* big << endl;
	liton_sp::debug::exec_except([&]() {x1.alloc(0, big, 0, 0, big, 0, 1, big, 0); }, out, err);
#endif
	liton_sp::debug::exec_except([&]() {x1.alloc(1, 2, 1, 0, 4, 0, 1, 5, 0); }, out, err);
	liton_sp::debug::exec_except([&]() {x1.alloc(1, 2, 1, 0, 4, 0, 1, 5, 0); }, out, err);
	out << endl;
	liton_sp::debug::exec_except([&]() {x1.realloc(1, 2, 1, 1, 5, 2, 1, 4, 2); }, out, err);
	liton_sp::debug::exec_except([&]() {x2.realloc(1, 2, 1, 1, 5, 2, 1, 4, 2); }, out, err);
	out << endl;
	liton_sp::debug::exec_except([&]() {x2.clear(); }, out, err);
	out << endl;

	x0.alloc(1, 2, 1, 1, 5, 2, 1, 4, 2);
	x2.alloc(1, 2, 1, 1, 5, 2, 1, 4, 2);
	x3.alloc(1, 2, 1, 1, 5, 2, 1, 4, 2);
	x4.alloc(1, 2, 1, 1, 5, 2, 1, 4, 2);
	x5.alloc(1, 2, 1, 1, 5, 2, 1, 4, 2);
	x6.alloc(1, 2, 1, 1, 5, 2, 1, 4, 2);
	x7.alloc(1, 2, 1, 1, 5, 2, 1, 4, 2);
	out << x0.disp() << endl;
	out << x1.disp() << endl;
	out << x2.disp() << endl;
	out << x3.disp() << endl;
	out << x4.disp() << endl;
	out << x5.disp() << endl;
	out << x6.disp() << endl;
	out << x7.disp() << endl;
	out << endl;

	D3::PD_For_N_3D(0, x0.N, x0.size().range(RA::ALL, RA::ALL, RA::ALL), [&]PD_F_n_ijk(n, i, j, k) {
		x0(n, i, j, k, FL::C, FL::C, FL::C) = static_cast<int>(n) * 1000 + i * 100 + j * 10 + k;
	});
	out << "x0:" << endl;
	out << x0.disp_data() << endl;

	D3::PD_For_N_3D(0, x1.N, x1.size().range(RA::ALL, RA::ALL, RA::ALL), [&]PD_F_n_ijk(n, i, j, k) {
		x1(n, i, j, k, FL::C, FL::C, FL::N) = static_cast<int>(n) * 1000 + i * 100 + j * 10 + k;
	});
	D3::PD_For_N_3D(0, x1.N, x1.size().range(RA::ALL, RA::ALL, RA::ALL), [&]PD_F_n_ijk(n, i, j, k) {
		x1(n, i, j, x1.size().last(2, RA::ALL), FL::C, FL::C, FL::P) = 0;
	});
	out << "x1:" << endl;
	out << x1.disp_data() << endl;

	D3::PD_For_N_3D(0, x2.N, x2.size().range(RA::ALL, RA::ALL, RA::ALL), [&]PD_F_n_ijk(n, i, j, k) {
		x2(n, i, j, k, FL::C, FL::N, FL::C) = static_cast<int>(n) * 1000 + i * 100 + j * 10 + k;
	});
	D3::PD_For_N_3D(0, x2.N, x2.size().range(RA::ALL, RA::ALL, RA::ALL), [&]PD_F_n_ijk(n, i, j, k) {
		x2(n, i, x2.size().last(1, RA::ALL), k, FL::C, FL::P, FL::C) = 0;
	});
	out << "x2:" << endl;
	out << x2.disp_data() << endl;

	D3::PD_For_N_3D(0, x3.N, x3.size().range(RA::ALL, RA::ALL, RA::ALL), [&]PD_F_n_ijk(n, i, j, k) {
		x3(n, i, j, k, FL::C, FL::N, FL::N) = static_cast<int>(n) * 1000 + i * 100 + j * 10 + k;
	});
	D3::PD_For_N_3D(0, x3.N, x3.size().range(RA::ALL, RA::ALL, RA::ALL), [&]PD_F_n_ijk(n, i, j, k) {
		x3(n, i, x3.size().last(1, RA::ALL), k, FL::C, FL::P, FL::N) = 1;
		x3(n, i, x3.size().last(1, RA::ALL), k, FL::C, FL::P, FL::P) = 1;
		x3(n, i, j, x3.size().last(2, RA::ALL), FL::C, FL::N, FL::P) = 2;
		x3(n, i, j, x3.size().last(2, RA::ALL), FL::C, FL::P, FL::P) = 2;
		x3(n, i, x3.size().last(1, RA::ALL), x3.size().last(2, RA::ALL), FL::C, FL::P, FL::P) = 3;
	});
	out << "x3:" << endl;
	out << x3.disp_data() << endl;

	D3::PD_For_N_3D(0, x4.N, x4.size().range(RA::ALL, RA::ALL, RA::ALL), [&]PD_F_n_ijk(n, i, j, k) {
		x4(n, i, j, k, FL::N, FL::C, FL::C) = static_cast<int>(n) * 1000 + i * 100 + j * 10 + k;
	});
	D3::PD_For_N_3D(0, x4.N, x4.size().range(RA::ALL, RA::ALL, RA::ALL), [&]PD_F_n_ijk(n, i, j, k) {
		x4(n, x4.size().last(0, RA::ALL), j, k, FL::P, FL::C, FL::C) = 0;
	});
	out << "x4:" << endl;
	out << x4.disp_data() << endl;

	D3::PD_For_N_3D(0, x5.N, x5.size().range(RA::ALL, RA::ALL, RA::ALL), [&]PD_F_n_ijk(n, i, j, k) {
		x5(n, i, j, k, FL::N, FL::C, FL::N) = static_cast<int>(n) * 1000 + i * 100 + j * 10 + k;
	});
	D3::PD_For_N_3D(0, x5.N, x5.size().range(RA::ALL, RA::ALL, RA::ALL), [&]PD_F_n_ijk(n, i, j, k) {
		x5(n, x5.size().last(0, RA::ALL), j, k, FL::P, FL::C, FL::N) = 1;
		x5(n, x5.size().last(0, RA::ALL), j, k, FL::P, FL::C, FL::P) = 1;
		x5(n, i, j, x5.size().last(2, RA::ALL), FL::N, FL::C, FL::P) = 2;
		x5(n, i, j, x5.size().last(2, RA::ALL), FL::P, FL::C, FL::P) = 2;
		x5(n, x5.size().last(0, RA::ALL), j, x5.size().last(2, RA::ALL), FL::P, FL::C, FL::P) = 3;
	});
	out << "x5:" << endl;
	out << x5.disp_data() << endl;

	D3::PD_For_N_3D(0, x6.N, x6.size().range(RA::ALL, RA::ALL, RA::ALL), [&]PD_F_n_ijk(n, i, j, k) {
		x6(n, i, j, k, FL::N, FL::N, FL::C) = static_cast<int>(n) * 1000 + i * 100 + j * 10 + k;
	});
	D3::PD_For_N_3D(0, x6.N, x6.size().range(RA::ALL, RA::ALL, RA::ALL), [&]PD_F_n_ijk(n, i, j, k) {
		x6(n, i, x6.size().last(1, RA::ALL), k, FL::N, FL::P, FL::C) = 1;
		x6(n, i, x6.size().last(1, RA::ALL), k, FL::P, FL::P, FL::C) = 1;
		x6(n, x6.size().last(0, RA::ALL), j, k, FL::P, FL::N, FL::C) = 2;
		x6(n, x6.size().last(0, RA::ALL), j, k, FL::P, FL::P, FL::C) = 2;
		x6(n, x6.size().last(0, RA::ALL), x6.size().last(1, RA::ALL), k, FL::P, FL::P, FL::C) = 3;
	});
	out << "x6:" << endl;
	out << x6.disp_data() << endl;

	D3::PD_For_N_3D(0, x7.N, x7.size().range(RA::ALL, RA::ALL, RA::ALL), [&]PD_F_n_ijk(n, i, j, k) {
		x7(n, i, j, k, FL::N, FL::N, FL::N) = static_cast<int>(n) * 1000 + i * 100 + j * 10 + k;
	});
	D3::PD_For_N_3D(0, x7.N, x7.size().range(RA::ALL, RA::ALL, RA::ALL), [&]PD_F_n_ijk(n, i, j, k) {
		x7(n, x7.size().last(0, RA::ALL), j, k, FL::P, FL::N, FL::N) = 1;
		x7(n, x7.size().last(0, RA::ALL), j, k, FL::P, FL::N, FL::P) = 1;
		x7(n, x7.size().last(0, RA::ALL), j, k, FL::P, FL::P, FL::N) = 1;
		x7(n, x7.size().last(0, RA::ALL), j, k, FL::P, FL::P, FL::P) = 1;
		x7(n, i, x7.size().last(1, RA::ALL), k, FL::N, FL::P, FL::N) = 2;
		x7(n, i, x7.size().last(1, RA::ALL), k, FL::N, FL::P, FL::P) = 2;
		x7(n, i, x7.size().last(1, RA::ALL), k, FL::P, FL::P, FL::N) = 2;
		x7(n, i, x7.size().last(1, RA::ALL), k, FL::P, FL::P, FL::P) = 2;
		x7(n, i, j, x7.size().last(2, RA::ALL), FL::N, FL::N, FL::P) = 3;
		x7(n, i, j, x7.size().last(2, RA::ALL), FL::N, FL::P, FL::P) = 3;
		x7(n, i, j, x7.size().last(2, RA::ALL), FL::P, FL::N, FL::P) = 3;
		x7(n, i, j, x7.size().last(2, RA::ALL), FL::P, FL::P, FL::P) = 3;
		x7(n, x7.size().last(0, RA::ALL), x7.size().last(1, RA::ALL), x7.size().last(2, RA::ALL), FL::P, FL::P, FL::P) = 4;
	});
	out << "x7:" << endl;
	out << x7.disp_data() << endl;

	liton_sp::debug::exec_except([&]() {out << x7(0, 0, 0, 0, FL::P, FL::P, FL::P) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x7(0, -2, 0, 0, FL::P, FL::P, FL::P) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x7(0, 3, 0, 0, FL::P, FL::P, FL::P) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x7(0, 0, -2, 0, FL::P, FL::P, FL::P) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x7(0, 0, 7, 0, FL::P, FL::P, FL::P) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x7(0, 0, 0, -2, FL::P, FL::P, FL::P) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x7(0, 0, 0, 6, FL::P, FL::P, FL::P) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x7(2, 0, 0, 0, FL::P, FL::P, FL::P) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x7(0, 0, 0, 0, FL::C, FL::P, FL::P) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x7(0, 0, 0, 0, FL::P, FL::C, FL::P) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x7(0, 0, 0, 0, FL::P, FL::P, FL::C) << endl; }, out, err);

	const D3::PointData<double, 1, LO::half, LO::half, LO::center> xc(0, 1, 0, 0, 5, 0, 0, 10, 0);
	liton_sp::debug::exec_except([&]() {out << xc(0, 0, 0, 0, FL::N, FL::N, FL::C) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << xc(0, 0, 0, 0) << endl; }, out, err);

	liton_sp::debug::exec_except([&]() {out << x7(0, -2, 0, 0, FL::N, FL::N, FL::N) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x7(0, -2, 0, 0, FL::P, FL::N, FL::N) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x7(0, 3, 0, 0, FL::N, FL::N, FL::N) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x7(0, 3, 0, 0, FL::P, FL::N, FL::N) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x7(0, 0, -2, 0, FL::N, FL::N, FL::N) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x7(0, 0, -2, 0, FL::N, FL::P, FL::N) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x7(0, 0, 7, 0, FL::N, FL::N, FL::N) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x7(0, 0, 7, 0, FL::N, FL::P, FL::N) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x7(0, 0, 0, -2, FL::N, FL::N, FL::N) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x7(0, 0, 0, -2, FL::N, FL::N, FL::P) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x7(0, 0, 0, 6, FL::N, FL::N, FL::N) << endl; }, out, err);
	liton_sp::debug::exec_except([&]() {out << x7(0, 0, 0, 6, FL::N, FL::N, FL::P) << endl; }, out, err);

	out.close();
	err.close();

	return 0;
}