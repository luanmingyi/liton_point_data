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
		class DIM
		{
		  public:
			static const int D = 0;

			inline static void check_d(const unsigned d)
			{
#ifdef _CHECK_POINTDATA_RANGE
				if (d >= D)
				{
					std::ostringstream errlog;
					errlog << "out of DIM range: " << d << " range:[" << 0 << "," << static_cast<int>
					       (D) - 1 << "]";
					throw(std::runtime_error(errlog.str()));
				}
#endif
			}
		};

		template <typename _NUMT, unsigned _N>
		class PointData
		{
		  public:
			typedef _NUMT num_type;
			const unsigned N = _N;
		  protected:
			_NUMT pt0[_N];

		  public:
			PointData();
			PointData(const PointData<_NUMT, _N> &) = delete;
			const PointData<_NUMT, _N> &operator=(const PointData<_NUMT, _N> &) = delete;

			std::string disp() const;
			std::string disp_data() const;

			inline _NUMT &operator()(const unsigned n)
			{
				this->check_n(n);
				return pt0[n];
			}
			inline const _NUMT &operator()(const unsigned n) const { return *this(n); }

		  protected:
			inline void check_n(const unsigned n) const
			{
#ifdef _CHECK_POINTDATA_RANGE
				if(n >= _N)
				{
					std::ostringstream errlog;
					errlog << "out of N range: " << n << " range:[" << 0 << "," << static_cast<int>(_N) - 1 << "]";
					throw(std::runtime_error(errlog.str()));
				}
#endif
			}
		};

		template <typename _NUMT, unsigned _N>
		PointData<_NUMT, _N>::PointData()
		{
			if(_N == 0)
			{
				throw(std::runtime_error("N is not allowed to be zero"));
			}
		}

		template <typename _NUMT, unsigned _N>
		inline std::string PointData<_NUMT, _N>::disp() const
		{
			std::ostringstream displog;
			displog << "dimension:" << DIM::D
			        << "  type:[" << typeid(_NUMT).name() << "]"
			        << "  N = " << _N;
			return displog.str();
		}

		template <typename _NUMT, unsigned _N>
		inline std::string PointData<_NUMT, _N>::disp_data() const
		{
			std::ostringstream displog;
			for(unsigned n = 0; n != _N; ++n)
			{
				displog << pt0[n] << ", ";
			}
			return displog.str();
		}
	}
}

#endif
