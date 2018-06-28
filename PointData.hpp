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
		class _C { public: static const int offset = 0; static const int sign = 0; };
		class _N { public: static const int offset = 0; static const int sign = -1; };
		class _P { public: static const int offset = 1; static const int sign = 1; };
		const _C C;
		const _N N;
		const _P P;
	}
	namespace RA
	{
		class _N {public: _N() {} _N(FL::_N fl) {}};
		class _P {public: _P() {} _P(FL::_P fl) {}};
		class _IN {public: _IN() {} _IN(FL::_C fl) {}};
		class _ALL {};
		const _N N;
		const _P P;
		const _IN IN;
		const _ALL ALL;
	}

	template <typename Function>
	inline void PD_For_N(const unsigned N_b, const unsigned N_e, const Function fun)
	{
		for (unsigned n = N_b; n != N_e; ++n)
		{
			fun(n);
		}
	}

	template <typename T, typename Reducer, typename Function>
	inline void PD_Reduce_N(const unsigned N_b, const unsigned N_e, T &ans, const Reducer reduce, const Function fun)
	{
		T temp = ans;
		for (unsigned n = N_b; n != N_e; ++n)
		{
			reduce(fun(n), temp);
		}
		reduce(temp, ans);
	}
}

#define PD_F_n(n_name) (const unsigned n_name)
#define PD_F_i(i_name) (const int i_name)
#define PD_F_ij(i_name, j_name) (const int i_name, const int j_name)
#define PD_F_ijk(i_name, j_name, k_name) (const int i_name, const int j_name, const int k_name)
#define PD_F_n_i(n_name, i_name) (const unsigned n_name, const int i_name)
#define PD_F_n_ij(n_name, i_name, j_name) (const unsigned n_name, const int i_name, const int j_name)
#define PD_F_n_ijk(n_name, i_name, j_name, k_name) (const unsigned n_name, const int i_name, const int j_name, const int k_name)

#define PD_RF(type, x_name, xx_name) (const type x_name, type &xx_name)
#define PD_RF_SUM(type) [](const type x, type &xx){ xx += x; }
#define PD_RF_SUM_ABS(type) [](const type x, type &xx){ xx += (x >= 0 ? x : -x); }
#define PD_RF_PRO(type) [](const type x, type &xx){ xx *= x; }
#define PD_RF_PRO_ABS(type) [](const type x, type &xx){ xx *= (x >= 0 ? x : -x); }
#define PD_RF_MAX(type) [](const type x, type &xx){ xx = xx < x ? x : xx; }
#define PD_RF_MAX_ABS(type) [](const type x, type &xx){ type x_abs = x >= 0 ? x : -x; xx = xx < x_abs ? x_abs : xx; }
#define PD_RF_MIN(type) [](const type x, type &xx){ xx = xx > x ? x : xx; }
#define PD_RF_MIN_ABS(type) [](const type x, type &xx){ type x_abs = x >= 0 ? x : -x; xx = xx > x_abs ? x_abs : xx; }

#include "PointData_0D.hpp"
#include "PointData_1D.hpp"
#include "PointData_2D.hpp"
#include "PointData_3D.hpp"

#endif
