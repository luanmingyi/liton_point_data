#ifndef POINTDATA_0D_HPP
#define POINTDATA_0D_HPP

#include <stdexcept>
#include <typeinfo>
#include <string>
#include <sstream>

namespace liton_pd
{
	namespace D0
	{
		static const int DIM = 0;

		inline void check_d(const unsigned &d)
		{
#ifdef _CHECK_POINTDATA_RANGE
			if(d >= DIM)
			{
				std::ostringstream errlog;
				errlog << "out of DIM range: d:[" << d << "] range:[" << 0 << "," << static_cast<int>
				       (DIM) - 1 << "]";
				throw(std::runtime_error(errlog.str()));
			}
#endif
		}

		template <typename _NUMT, LO::LOCATION _LOC = LO::center, unsigned _N = 1>
		class PointData
		{
		  public:
			typedef _NUMT num_type;
			const LO::LOCATION LOC = _LOC;
			const unsigned N = _N;
		  protected:
			_NUMT pt0[_N];

		  public:
			PointData();
			PointData(const PointData<_NUMT, _LOC, _N> &) = delete;
			const PointData<_NUMT, _LOC, _N> &operator=(const PointData<_NUMT, _LOC, _N> &) = delete;

			std::string disp() const;
			std::string disp_data() const;
		  protected:
			inline void check_n(const unsigned &n) const
			{
#ifdef _CHECK_POINTDATA_RANGE
				if(n >= _N)
				{
					std::ostringstream errlog;
					errlog << "out of N range: n:[" << n << "] range:[" << 0 << "," << static_cast<int>(_N) - 1 << "]";
					throw(std::runtime_error(errlog.str()));
				}
#endif
			}
		};

		template <typename _NUMT, unsigned _N = 1>
		class PointData_C : public PointData<_NUMT, LO::center, _N>
		{
		  public:
			using PointData<_NUMT, LO::center, _N>::num_type;
			using PointData<_NUMT, LO::center, _N>::LOC;
			using PointData<_NUMT, LO::center, _N>::N;
		  protected:
			using PointData<_NUMT, LO::center, _N>::pt0;
		  public:
			PointData_C(): PointData<_NUMT, LO::center, _N>() {}
			PointData_C(const PointData_C<_NUMT, _N> &) = delete;
			const PointData_C<_NUMT, _N> &operator=(const PointData_C<_NUMT, _N> &) = delete;

			inline _NUMT &operator()(const unsigned &n)
			{
				this->check_n(n);
				return pt0[n];
			}
			inline const _NUMT &operator()(const int &i, const unsigned &n) const { return *this(n); }
		};

		template <typename Function>
		inline void For_PD_0D_N(const unsigned &N_b, const unsigned &N_e, const Function &fun)
		{For_N(N_b, N_e, fun);}
		template <typename T, typename Reducer, typename Function>
		inline void Reduce_PD_0D_N(const unsigned &N_b, const unsigned &N_e, T &ans, const Reducer &reduce, const Function &fun)
		{Reduce_N(N_b, N_e, ans, reduce, fun);}

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		template <typename _NUMT, LO::LOCATION _LOC, unsigned _N>
		PointData<_NUMT, _LOC, _N>::PointData()
		{
			if(_N == 0)
			{
				throw(std::runtime_error("N is not allowed to be zero"));
			}
		}

		template<typename _NUMT, LO::LOCATION _LOC, unsigned _N>
		inline std::string PointData<_NUMT, _LOC, _N>::disp() const
		{
			char loc_str[2][50] = { "center\0", "half\0" };
			std::ostringstream displog;
			displog << "dimension:" << DIM
			        << "  location:[" << loc_str[_LOC] << "]"
			        << "  type:[" << typeid(_NUMT).name() << "]"
			        << "  N = " << _N;
			return displog.str();
		}

		template<typename _NUMT, LO::LOCATION _LOC, unsigned _N>
		inline std::string PointData<_NUMT, _LOC, _N>::disp_data() const
		{
			std::ostringstream displog;
			for(unsigned n = 0; n != _N; ++n)
			{
				displog << pt0[n] << ", ";
			}
			displog << std::endl;
			return displog.str();
		}
	}
}

#endif
