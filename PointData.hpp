#ifndef POINTDATA_H
#define POINTDATA_H

namespace liton_pd
{
	namespace LO
	{
		enum LOCATION
		{
			center = 0,
			half = 1
		};
	}
	namespace FL
	{
		enum FLAG
		{
			N = 0,
			P = 1
		};
	}
	namespace RA
	{
		class _N {} N;
		class _P {} P;
		class _ALL {} ALL;
		class _IN {} IN;
	}

	template <typename Function>
	inline void For_N(const unsigned &N_b, const unsigned &N_e, const Function &fun)
	{
		for (unsigned n = N_b; n != N_e; ++n)
		{
			fun(n);
		}
	}

	template <typename T, typename Reducer, typename Function>
	inline void Reduce_N(const unsigned &N_b, const unsigned &N_e, T &ans, const Reducer &reduce, const Function &fun)
	{
		T temp = ans;
		for (unsigned n = N_b; n != N_e; ++n)
		{
			reduce(fun(n), temp);
		}
		reduce(temp, ans);
	}
}

#define PD_F_n(n_name) (const unsigned &n_name)
#define PD_F_i(i_name) (const int &i_name)
#define PD_F_ij(i_name, j_name) (const int &i_name, const int &j_name)
#define PD_F_ijk(i_name, j_name, k_name) (const int &i_name, const int &j_name, const int &k_name)
#define PD_F_i_n(i_name, n_name) (const int &i_name, const unsigned &n_name)
#define PD_F_ij_n(i_name, j_name, n_name) (const int &i_name, const int &j_name, const unsigned &n_name)
#define PD_F_ijk_n(i_name, j_name, k_name, n_name) (const int &i_name, const int &j_name, const int &k_name, const unsigned &n_name)

#define PD_RF(type, x_name, xx_name) (const type &x_name, type &xx_name)

#include "PointData_0D.hpp"
#include "PointData_1D.hpp"
#include "PointData_2D.hpp"

#endif
