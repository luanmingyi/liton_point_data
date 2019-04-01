#include<iostream>
using namespace std;

#ifdef _DEBUG
	#define _CHECK_POINTDATA_RANGE
#endif
#include "PointData.hpp"

#include "Test_HPP_sub.h"

int main(int argc, char** argv)
{
	liton_pd::D1::DIM::check_d(0);
	cout << liton_pd::D1::DIM::D << endl;
	test_hpp a;
	a.init_and_disp();
	return 0;
}