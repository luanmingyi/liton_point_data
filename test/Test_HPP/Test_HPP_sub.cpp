#include "Test_HPP_sub.h"
#include <iostream>

#ifdef _DEBUG
#define _CHECK_POINTDATA_RANGE
#endif
#include "../../scr/liton_point_data/PointData.hpp"
using namespace liton_pd;

void test_hpp::init_and_disp()
{
	x.alloc(2, 100, 2);
	liton_pd::D1::DIM::check_d(0);
	std::cout << D1::DIM::D << std::endl;
	std::cout << x.disp() << std::endl;
}
