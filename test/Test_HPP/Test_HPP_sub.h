#ifndef TEST_HPP_SUB_H
#define TEST_HPP_SUB_H

#include "../../scr/liton_point_data/PointData.hpp"
using namespace liton_pd;

class test_hpp
{
public:
	D1::PointData<double, 1, LO::center> x;
public:
	void init_and_disp();
};

#endif
