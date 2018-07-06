#ifndef POINTDATA_3D_HPP
#define POINTDATA_3D_HPP

#include <stdexcept>
#include <typeinfo>
#include <string>
#include <sstream>

#ifdef PD_OT
	#include "../liton_ordered_tec/ordered_tec.h"
	#include <algorithm>
#endif

namespace liton_pd
{
	namespace D3
	{
		class DIM
		{
		  public:
			static const int D = 3;

			inline static void check_d(const unsigned d)
			{
#ifdef _CHECK_POINTDATA_RANGE
				if (d >= D)
				{
					std::ostringstream errlog;
					errlog << "DIM is out of range: " << d << " range:[" << 0 << "," << static_cast<int>
					       (D) - 1 << "]";
					throw(std::runtime_error(errlog.str()));
				}
#endif
			}
		};
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		class RangeT
		{
		  public:
			int _begin[DIM::D] = { 0, 0, 0 };
			int _size[DIM::D] = { 0, 0, 0 };

		  public:
			RangeT() = default;
			RangeT(const int b0, const unsigned s0,
			       const int b1, const unsigned s1,
			       const int b2, const unsigned s2):
				_begin{ b0, b1, b2 },
				_size{ static_cast<int>(s0), static_cast<int>(s1), static_cast<int>(s2) } {}

			inline const int* begin_pt() const { return _begin; }
			inline const int* size_pt() const { return _size; }

			inline int begin(const unsigned d) const { DIM::check_d(d); return _begin[d]; }
			inline int size(const unsigned d) const { DIM::check_d(d); return _size[d]; }
			inline int end(const unsigned d) const { return begin(d) + size(d); }
			inline int last(const unsigned d) const { return end(d) - 1; }
			inline int bound(const unsigned d, FL::_N fl) const { return begin(d); }
			inline int bound(const unsigned d, FL::_P fl) const { return last(d); }

			inline RangeT &cut_head(int i, int j, int k)
			{
				_begin[0] += i; _size[0] -= i;
				_begin[1] += j; _size[1] -= j;
				_begin[2] += k; _size[2] -= k;
				return *this;
			}
			inline RangeT &cut_tail(int i, int j, int k) { _size[0] -= i; _size[1] -= j; _size[2] -= k; return *this; }
			inline RangeT &tran(int i, int j, int k) { _begin[0] += i; _begin[1] += j; _begin[2] += k; return *this; }

			inline bool is_overlap(const RangeT r) const
			{
				if(begin(0) >= r.end(0) || end(0) <= r.begin(0)
				        || begin(1) >= r.end(1) || end(1) <= r.begin(1)
				        || begin(2) >= r.end(2) || end(2) <= r.begin(2))
				{
					return false;
				}
				else
				{
					return true;
				}
			}
			inline RangeT overlap(const RangeT r) const
			{
#ifdef _CHECK_POINTDATA_RANGE
				if(!is_overlap(r))
				{
					throw(std::runtime_error("overlap not found"));
				}
#endif
				int __begin0 = begin(0) > r.begin(0) ? begin(0) : r.begin(0);
				int __end0 = end(0) < r.end(0) ? end(0) : r.end(0);
				int __begin1 = begin(1) > r.begin(1) ? begin(1) : r.begin(1);
				int __end1 = end(1) < r.end(1) ? end(1) : r.end(1);
				int __begin2 = begin(2) > r.begin(2) ? begin(2) : r.begin(2);
				int __end2 = end(2) < r.end(2) ? end(2) : r.end(2);
				return RangeT(__begin0, __end0 - __begin0,
				              __begin1, __end1 - __begin1,
				              __begin2, __end2 - __begin2);
			}

			inline void check_range(const unsigned d, const LO::LOCATION loc, const int ii, const int offset) const
			{
#ifdef _CHECK_POINTDATA_RANGE
				DIM::check_d(d);
				if (ii < begin(d) - offset || ii >= end(d) + loc - offset)
				{
					char ijk_str[3][50] = { "i\0", "j\0", "k\0" };
					std::ostringstream errlog;
					if (loc == LO::center)
					{
						errlog << ijk_str[d] << " is out of range: " << ii << " range:[" << begin(d) << "," << last(d) << "]";
					}
					else
					{
						char flag_str[2][50] = { "N\0", "P\0" };
						errlog << ijk_str[d] << " is out of range: " << ii << "(" << flag_str[offset] << ")"
						       << " range:[" << begin(d) - 1 << "(P)," << begin(d) << "," << last(d) << "," << end(d) << "(N)]";
					}
					throw(std::runtime_error(errlog.str()));
				}
#endif
			}
			inline void check_range(const unsigned d, const int ii) const { check_range(d, LO::center, ii, 0); }

			inline std::string disp() const
			{
				std::ostringstream displog;
				displog << "range[0] = [" << begin(0) << " , " << last(0) << "]  "
				        << "size[0] = " << size(0) << "    ";
				displog << "range[1] = [" << begin(1) << " , " << last(1) << "]  "
				        << "size[1] = " << size(1) << "    ";
				displog << "range[2] = [" << begin(2) << " , " << last(2) << "]  "
				        << "size[2] = " << size(2);
				return displog.str();
			}
		};
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		template <typename Function>
		inline void PD_For_3D(const RangeT r, const Function fun)
		{
			const int begin0 = r.begin(0);
			const int end0 = r.end(0);
			const int begin1 = r.begin(1);
			const int end1 = r.end(1);
			const int begin2 = r.begin(2);
			const int end2 = r.end(2);
			for (int i = begin0; i != end0; ++i)
			{
				for (int j = begin1; j != end1; ++j)
				{
					for (int k = begin2; k != end2; ++k)
					{
						fun(i, j, k);
					}
				}
			}
		}

		template <typename Function>
		inline void PD_For_N_3D(const unsigned N_b, const unsigned N_e, const RangeT &r, const Function fun)
		{
			const int begin0 = r.begin(0);
			const int end0 = r.end(0);
			const int begin1 = r.begin(1);
			const int end1 = r.end(1);
			const int begin2 = r.begin(2);
			const int end2 = r.end(2);
			for (unsigned n = N_b; n != N_e; ++n)
			{
				for (int i = begin0; i != end0; ++i)
				{
					for (int j = begin1; j != end1; ++j)
					{
						for (int k = begin2; k != end2; ++k)
						{
							fun(n, i, j, k);
						}
					}
				}
			}
		}

		template <typename T, typename Reducer, typename Function>
		inline void PD_Reduce_3D(const RangeT r, T &ans, const Reducer reduce, const Function fun)
		{
			T temp = ans;
			const int begin0 = r.begin(0);
			const int end0 = r.end(0);
			const int begin1 = r.begin(1);
			const int end1 = r.end(1);
			const int begin2 = r.begin(2);
			const int end2 = r.end(2);
			for (int i = begin0; i != end0; ++i)
			{
				for (int j = begin1; j != end1; ++j)
				{
					for (int k = begin2; k != end2; ++k)
					{
						reduce(fun(i, j, k), temp);
					}
				}
			}
			reduce(temp, ans);
		}

		template <typename T, typename Reducer, typename Function>
		inline void PD_Reduce_N_3D(const unsigned N_b, const unsigned N_e, const RangeT r, T &ans, const Reducer reduce, const Function fun)
		{
			T temp = ans;
			const int begin0 = r.begin(0);
			const int end0 = r.end(0);
			const int begin1 = r.begin(1);
			const int end1 = r.end(1);
			const int begin2 = r.begin(2);
			const int end2 = r.end(2);
			for (unsigned n = N_b; n != N_e; ++n)
			{
				for (int i = begin0; i != end0; ++i)
				{
					for (int j = begin1; j != end1; ++j)
					{
						for (int k = begin2; k != end2; ++k)
						{
							reduce(fun(n, i, j, k), temp);
						}
					}
				}
			}
			reduce(temp, ans);
		}

		template <typename Function>
		inline void PD_For_2D(const unsigned d0, const unsigned d1, const RangeT r, const Function fun)
		{
			const int begin0 = r.begin(d0);
			const int end0 = r.end(d0);
			const int begin1 = r.begin(d1);
			const int end1 = r.end(d1);
			for (int ii = begin0; ii != end0; ++ii)
			{
				for (int jj = begin1; jj != end1; ++jj)
				{
					fun(ii, jj);
				}
			}
		}

		template <typename Function>
		inline void PD_For_N_2D(const unsigned N_b, const unsigned N_e, const unsigned d0, const unsigned d1, const RangeT &r, const Function fun)
		{
			const int begin0 = r.begin(d0);
			const int end0 = r.end(d0);
			const int begin1 = r.begin(d1);
			const int end1 = r.end(d1);
			for (unsigned n = N_b; n != N_e; ++n)
			{
				for (int ii = begin0; ii != end0; ++ii)
				{
					for (int jj = begin1; jj != end1; ++jj)
					{
						fun(n, ii, jj);
					}
				}
			}
		}

		template <typename T, typename Reducer, typename Function>
		inline void PD_Reduce_2D(const unsigned d0, const unsigned d1, const RangeT r, T &ans, const Reducer reduce, const Function fun)
		{
			T temp = ans;
			const int begin0 = r.begin(d0);
			const int end0 = r.end(d0);
			const int begin1 = r.begin(d1);
			const int end1 = r.end(d1);
			for (int ii = begin0; ii != end0; ++ii)
			{
				for (int jj = begin1; jj != end1; ++jj)
				{
					reduce(fun(ii, jj), temp);
				}
			}
			reduce(temp, ans);
		}

		template <typename T, typename Reducer, typename Function>
		inline void PD_Reduce_N_2D(const unsigned N_b, const unsigned N_e, const unsigned d0, const unsigned d1, const RangeT r, T &ans, const Reducer reduce, const Function fun)
		{
			T temp = ans;
			const int begin0 = r.begin(d0);
			const int end0 = r.end(d0);
			const int begin1 = r.begin(d1);
			const int end1 = r.end(d1);
			for (unsigned n = N_b; n != N_e; ++n)
			{
				for (int ii = begin0; ii != end0; ++ii)
				{
					for (int jj = begin1; jj != end1; ++jj)
					{
						reduce(fun(n, ii, jj), temp);
					}
				}
			}
			reduce(temp, ans);
		}

		template <typename Function>
		inline void PD_For_1D(const unsigned d0, const RangeT r, const Function fun)
		{
			const int begin0 = r.begin(d0);
			const int end0 = r.end(d0);
			for (int ii = begin0; ii != end0; ++ii)
			{
				fun(ii);
			}
		}

		template <typename Function>
		inline void PD_For_N_1D(const unsigned N_b, const unsigned N_e, const unsigned d0, const RangeT r, const Function fun)
		{
			const int begin0 = r.begin(d0);
			const int end0 = r.end(d0);
			for (unsigned n = N_b; n != N_e; ++n)
			{
				for (int ii = begin0; ii != end0; ++ii)
				{
					fun(n, ii);
				}
			}
		}

		template <typename T, typename Reducer, typename Function>
		inline void PD_Reduce_1D(const unsigned d0, const RangeT r, T &ans, const Reducer reduce, const Function fun)
		{
			T temp = ans;
			const int begin0 = r.begin(d0);
			const int end0 = r.end(d0);
			for (int ii = begin0; ii != end0; ++ii)
			{
				reduce(fun(ii), temp);
			}
			reduce(temp, ans);
		}

		template <typename T, typename Reducer, typename Function>
		inline void PD_Reduce_N_1D(const unsigned N_b, const unsigned N_e, const unsigned d0, const RangeT &r, T &ans, const Reducer reduce, const Function fun)
		{
			T temp = ans;
			const int begin0 = r.begin(d0);
			const int end0 = r.end(d0);
			for (unsigned n = N_b; n != N_e; ++n)
			{
				for (int ii = begin0; ii != end0; ++ii)
				{
					reduce(fun(n, ii), temp);
				}
			}
			reduce(temp, ans);
		}
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		class SizeT
		{
		  public:
			int _in[DIM::D] = { 0, 0, 0 };
			int _n[DIM::D] = { 0, 0, 0 };
			int _p[DIM::D] = { 0, 0, 0 };

		  public:
			SizeT() = default;
			SizeT(const unsigned in0, const unsigned iin0, const unsigned ip0,
			      const unsigned in1, const unsigned iin1, const unsigned ip1,
			      const unsigned in2, const unsigned iin2, const unsigned ip2):
				_in{ static_cast<int>(iin0), static_cast<int>(iin1), static_cast<int>(iin2) },
				_n{ static_cast<int>(in0), static_cast<int>(in1), static_cast<int>(in2) },
				_p{ static_cast<int>(ip0), static_cast<int>(ip1), static_cast<int>(ip2) } {}

			inline int in(const unsigned d) const { DIM::check_d(d); return _in[d]; }
			inline int n(const unsigned d) const { DIM::check_d(d); return _n[d]; }
			inline int p(const unsigned d) const { DIM::check_d(d); return _p[d]; }

			inline int begin(const unsigned d, RA::_ALL r) const { DIM::check_d(d); return -_n[d]; }
			inline int end(const unsigned d, RA::_ALL r) const { DIM::check_d(d); return _in[d] + _p[d]; }
			inline int begin(const unsigned d, RA::_IN r) const { DIM::check_d(d); return 0; }
			inline int end(const unsigned d, RA::_IN r) const { DIM::check_d(d); return _in[d]; }
			inline int begin(const unsigned d, RA::_N r) const { DIM::check_d(d); return -_n[d]; }
			inline int end(const unsigned d, RA::_N r) const { DIM::check_d(d); return 0; }
			inline int begin(const unsigned d, RA::_P r) const { DIM::check_d(d); return _in[d]; }
			inline int end(const unsigned d, RA::_P r) const { DIM::check_d(d); return _in[d] + _p[d]; }
			template<typename T0>
			inline unsigned int size(const unsigned d, T0 r) const { return end(d, r) - begin(d, r); }
			template<typename T0>
			inline int last(const unsigned d, T0 r) const { return end(d, r) - 1; }
			template<typename T0>
			inline int bound(const unsigned d, T0 r, FL::_N fl) const { return begin(d, r); }
			template<typename T0>
			inline int bound(const unsigned d, T0 r, FL::_P fl) const { return last(d, r); }

			template<typename _FL>
			inline int mirror(const unsigned d, _FL fl, int ii) const { return 2 * bound(d, RA::IN, fl) - ii; }
			template<typename _FL>
			inline int periodic(const unsigned d, _FL fl, int ii) const { return -_FL::sign * last(d, RA::IN) + ii; }

			template<typename T0, typename T1, typename T2>
			inline RangeT range(T0 r0, T1 r1, T2 r2) const
			{
				return RangeT(begin(0, r0), size(0, r0),
				              begin(1, r1), size(1, r1),
				              begin(2, r2), size(2, r2));
			}

			std::string disp() const
			{
				std::ostringstream displog;
				displog << "size[0] = [" << _n[0] << " " << _in[0] << " " << _p[0] << "]" << "  ";
				displog << "size[1] = [" << _n[1] << " " << _in[1] << " " << _p[1] << "]" << "  ";
				displog << "size[2] = [" << _n[2] << " " << _in[2] << " " << _p[2] << "]";
				return displog.str();
			}

			inline void check() const
			{
				if (_in[0] == 0 && _n[0] + _p[0] != 0)
				{
					throw(std::runtime_error("dim[0]: size_in can not be zero when size_n or size_p is non-zero"));
				}
				if (_in[1] == 0 && _n[1] + _p[1] != 0)
				{
					throw(std::runtime_error("dim[1]: size_in can not be zero when size_n or size_p is non-zero"));
				}
				if (_in[2] == 0 && _n[2] + _p[2] != 0)
				{
					throw(std::runtime_error("dim[2]: size_in can not be zero when size_n or size_p is non-zero"));
				}
			}
		};
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		template <typename _NUMT, unsigned _N, LO::LOCATION _LOC0, LO::LOCATION _LOC1, LO::LOCATION _LOC2>
		class PointData
		{
		  public:
			typedef _NUMT num_type;
			const LO::LOCATION LOC0 = _LOC0;
			const LO::LOCATION LOC1 = _LOC1;
			const LO::LOCATION LOC2 = _LOC2;
			const unsigned N = _N;
		  protected:
			SizeT _size;
			_NUMT* data = nullptr;
			_NUMT*** pt0[_N];
			_NUMT*** pt1 = nullptr;
			_NUMT** pt2 = nullptr;

		  public:
			PointData();
			PointData(const unsigned in0, const unsigned iin0, const unsigned ip0,
			          const unsigned in1, const unsigned iin1, const unsigned ip1,
			          const unsigned in2, const unsigned iin2, const unsigned ip2);
			PointData(const PointData<_NUMT, _N, _LOC0, _LOC1, _LOC2> &) = delete;
			const PointData<_NUMT, _N, _LOC0, _LOC1, _LOC2> &operator=(const PointData<_NUMT, _N, _LOC0, _LOC1, _LOC2> &) = delete;
			~PointData();

			void alloc(const unsigned in0, const unsigned iin0, const unsigned ip0,
			           const unsigned in1, const unsigned iin1, const unsigned ip1,
			           const unsigned in2, const unsigned iin2, const unsigned ip2);
			void realloc(const unsigned in0, const unsigned iin0, const unsigned ip0,
			             const unsigned in1, const unsigned iin1, const unsigned ip1,
			             const unsigned in2, const unsigned iin2, const unsigned ip2);
			void clear();

			std::string disp() const;
			std::string disp_data() const;

			inline SizeT size() const { return _size; }

			template<typename F0 = FL::_C, typename F1 = FL::_C, typename F2 = FL::_C>
			inline _NUMT & operator()(const unsigned n, const int i, const int j, const int k, const F0 flag0 = FL::C, const F1 flag1 = FL::C, const F2 flag2 = FL::C)
			{
				check_data();
				check_n(n);
				check_flag(flag0, flag1, flag2);
				_size.range(RA::ALL, RA::ALL, RA::ALL).check_range(0, LOC0, i, F0::offset);
				_size.range(RA::ALL, RA::ALL, RA::ALL).check_range(1, LOC1, j, F1::offset);
				_size.range(RA::ALL, RA::ALL, RA::ALL).check_range(2, LOC2, k, F2::offset);
				return pt0[n][i + F0::offset][j + F1::offset][k + F2::offset];
			}
			template<typename F0 = FL::_C, typename F1 = FL::_C, typename F2 = FL::_C>
			inline const _NUMT & operator()(const unsigned n, const int i, const int j, const int k, const F0 flag0 = FL::C, const F1 flag1 = FL::C, const F2 flag2 = FL::C) const
			{
				check_data();
				check_n(n);
				check_flag(flag0, flag1, flag2);
				_size.range(RA::ALL, RA::ALL, RA::ALL).check_range(0, LOC0, i, F0::offset);
				_size.range(RA::ALL, RA::ALL, RA::ALL).check_range(1, LOC1, j, F1::offset);
				_size.range(RA::ALL, RA::ALL, RA::ALL).check_range(2, LOC2, k, F2::offset);
				return pt0[n][i + F0::offset][j + F1::offset][k + F2::offset];
			}

			inline _NUMT* data_pt(const unsigned n) { check_n(n); return pt0[n][-_size.n(0)][-_size.n(1)] - _size.n(2); }
			inline const _NUMT* data_pt(const unsigned n) const { check_n(n); return pt0[n][-_size.n(0)][-_size.n(1)] - _size.n(2); }
			inline _NUMT*** _pt0(const unsigned n) { check_n(n); return pt0[n]; }
			inline const _NUMT*** _pt0(const unsigned n) const { check_n(n); return const_cast<const _NUMT*** >(pt0[n]); }

			void copy_from(const PointData<_NUMT, _N, _LOC0, _LOC1, _LOC2> &pd);
			template<typename __NUMT, unsigned __N>
			void copy_from(const PointData<__NUMT, __N, _LOC0, _LOC1, _LOC2> &pd,
			               const RangeT R, const RangeT r_,
			               unsigned N, unsigned n, int I, int J, int K, int i, int j, int k);
			template<typename __NUMT, unsigned __N,
			         typename _R0, typename _R1, typename _R2, typename _r0, typename _r1, typename _r2,
			         typename _FL0, typename _FL1, typename _FL2>
			inline void copy_from(const PointData<__NUMT, __N, _LOC0, _LOC1, _LOC2> &pd,
			                      _R0 R0, _R1 R1, _R2 R2, _r0 r0, _r1 r1, _r2 r2,
			                      unsigned N, unsigned n, _FL0 fl0, _FL1 fl1, _FL2 fl2)
			{
				copy_from(pd,
				          _size.range(R0, R1, R2), pd.size().range(r0, r1, r2),
				          N, n,
				          _size.bound(0, R0, fl0), _size.bound(1, R1, fl1), _size.bound(2, R2, fl2),
				          pd.size().bound(0, r0, fl0), pd.size().bound(1, r1, fl1), pd.size().bound(2, r2, fl2));
			}

#ifdef PD_OT
			template<typename _R0, typename _R1, typename _R2, typename _FL0, typename _FL1, typename _FL2>
			void read_plt(const std::string &root, const liton_ot::TEC_FILE_LOG &teclog,
			              unsigned zone, const std::string &name,
			              _R0 R0, _R1 R1, _R2 R2, unsigned nn, _FL0 fl0, _FL1 fl1, _FL2 fl2);
#endif

			inline void check_n(const unsigned n) const
			{
#ifdef _CHECK_POINTDATA_RANGE
				if(n >= _N)
				{
					std::ostringstream errlog;
					errlog << "N is out of range: " << n << " range:[" << 0 << "," << static_cast<int>(_N) - 1 << "]";
					throw(std::runtime_error(errlog.str()));
				}
#endif
			}

			template<typename F0, typename F1, typename F2>
			inline void check_flag(const F0 flag0, const F1 flag1, const F2 flag2) const
			{
#ifdef _CHECK_POINTDATA_RANGE
				if (_LOC0 == LO::center && typeid(F0) != typeid(FL::_C))
				{
					throw(std::runtime_error("dim[0]: flag must be [C] when location is [center]"));
				}
				if (_LOC0 == LO::half && typeid(F0) == typeid(FL::_C))
				{
					throw(std::runtime_error("dim[0]: flag must be [N] or [P] when location is [half]"));
				}
				if (_LOC1 == LO::center && typeid(F1) != typeid(FL::_C))
				{
					throw(std::runtime_error("dim[1]: flag must be [C] when location is [center]"));
				}
				if (_LOC1 == LO::half && typeid(F1) == typeid(FL::_C))
				{
					throw(std::runtime_error("dim[1]: flag must be [N] or [P] when location is [half]"));
				}
				if (_LOC2 == LO::center && typeid(F2) != typeid(FL::_C))
				{
					throw(std::runtime_error("dim[2]: flag must be [C] when location is [center]"));
				}
				if (_LOC2 == LO::half && typeid(F2) == typeid(FL::_C))
				{
					throw(std::runtime_error("dim[2]: flag must be [N] or [P] when location is [half]"));
				}
#endif
			}

			inline void check_data() const
			{
#ifdef _CHECK_POINTDATA_RANGE
				if (data == nullptr)
				{
					throw(std::runtime_error("data do not exist"));
				}
#endif
			}
		};

		template <typename _NUMT, unsigned _N, LO::LOCATION _LOC0, LO::LOCATION _LOC1, LO::LOCATION _LOC2>
		PointData<_NUMT, _N, _LOC0, _LOC1, _LOC2>::PointData()
		{
			if(_N == 0)
			{
				throw(std::runtime_error("N is not allowed to be zero"));
			}
			for(unsigned n = 0; n != _N; ++n)
			{
				pt0[n] = nullptr;
			}
		}

		template <typename _NUMT, unsigned _N, LO::LOCATION _LOC0, LO::LOCATION _LOC1, LO::LOCATION _LOC2>
		PointData<_NUMT, _N, _LOC0, _LOC1, _LOC2>::PointData(const unsigned in0, const unsigned iin0, const unsigned ip0,
		        const unsigned in1, const unsigned iin1, const unsigned ip1,
		        const unsigned in2, const unsigned iin2, const unsigned ip2)
		{
			if(_N == 0)
			{
				throw(std::runtime_error("N is not allowed to be zero"));
			}
			for(unsigned n = 0; n != _N; ++n)
			{
				pt0[n] = nullptr;
			}
			alloc(in0, iin0, ip0, in1, iin1, ip1, in2, iin2, ip2);
		}

		template <typename _NUMT, unsigned _N, LO::LOCATION _LOC0, LO::LOCATION _LOC1, LO::LOCATION _LOC2>
		PointData<_NUMT, _N, _LOC0, _LOC1, _LOC2>::~PointData()
		{
			if(data != nullptr)
			{
				clear();
			}
		}

		template <typename _NUMT, unsigned _N, LO::LOCATION _LOC0, LO::LOCATION _LOC1, LO::LOCATION _LOC2>
		void PointData<_NUMT, _N, _LOC0, _LOC1, _LOC2>::alloc(const unsigned in0, const unsigned iin0, const unsigned ip0,
		        const unsigned in1, const unsigned iin1, const unsigned ip1,
		        const unsigned in2, const unsigned iin2, const unsigned ip2)
		{
			if(data != nullptr)
			{
				throw(std::runtime_error("allocating memory do not allowed when data is not empty"));
			}
			else
			{
				SizeT s(in0, iin0, ip0, in1, iin1, ip1, in2, iin2, ip2);
				s.check();
				unsigned size_point = (s.size(0, RA::ALL) + _LOC0) * (s.size(1, RA::ALL) + _LOC1) * (s.size(2, RA::ALL) + _LOC2);
				if(size_point != 0)
				{
					try
					{
						data = new _NUMT[size_point * _N];
						_size = s;

						pt1 = new _NUMT** [N * (s.size(0, RA::ALL) + _LOC0)];
						pt2 = new _NUMT*[N * (s.size(0, RA::ALL) + _LOC0) * (s.size(1, RA::ALL) + _LOC1)];
						for (unsigned n = 0; n != _N; ++n)
						{
							pt0[n] = pt1 + n * (s.size(0, RA::ALL) + _LOC0);
							for (unsigned ii = 0; ii != s.size(0, RA::ALL) + _LOC0; ++ii)
							{
								pt0[n][ii] = pt2 + n * (s.size(0, RA::ALL) + _LOC0) * (s.size(1, RA::ALL) + _LOC1)
								             + ii * (s.size(1, RA::ALL) + _LOC1);
								for (unsigned jj = 0; jj != s.size(1, RA::ALL) + _LOC1; ++jj)
								{
									pt0[n][ii][jj] = data + n * (s.size(0, RA::ALL) + _LOC0) * (s.size(1, RA::ALL) + _LOC1) * (s.size(2, RA::ALL) + _LOC2)
									                 + ii * (s.size(1, RA::ALL) + _LOC1) * (s.size(2, RA::ALL) + _LOC2)
									                 + jj * (s.size(2, RA::ALL) + _LOC2);
									pt0[n][ii][jj] += s.n(2);
								}
								pt0[n][ii] += s.n(1);
							}
							pt0[n] += s.n(0);
						}
					}
					catch(const std::bad_alloc &)
					{
						throw(std::runtime_error("no more memory for allocating"));
					}
				}
				else
				{
					return;
				}
			}
		}

		template <typename _NUMT, unsigned _N, LO::LOCATION _LOC0, LO::LOCATION _LOC1, LO::LOCATION _LOC2>
		void PointData<_NUMT, _N, _LOC0, _LOC1, _LOC2>::realloc(const unsigned in0, const unsigned iin0, const unsigned ip0,
		        const unsigned in1, const unsigned iin1, const unsigned ip1,
		        const unsigned in2, const unsigned iin2, const unsigned ip2)
		{
			SizeT(in0, iin0, ip0, in1, iin1, ip1, in2, iin2, ip2).check();
			clear();
			alloc(in0, iin0, ip0, in1, iin1, ip1, in2, iin2, ip2);
		}

		template <typename _NUMT, unsigned _N, LO::LOCATION _LOC0, LO::LOCATION _LOC1, LO::LOCATION _LOC2>
		void PointData<_NUMT, _N, _LOC0, _LOC1, _LOC2>::clear()
		{
			if(data == nullptr)
			{
				throw(std::runtime_error("clearing data do not allowed when data is empty"));
			}
			else
			{
				delete[] data;
				data = nullptr;
				delete[] pt1;
				pt1 = nullptr;
				delete[] pt2;
				pt2 = nullptr;
				_size = SizeT();
			}
		}

		template <typename _NUMT, unsigned _N, LO::LOCATION _LOC0, LO::LOCATION _LOC1, LO::LOCATION _LOC2>
		inline std::string PointData<_NUMT, _N, _LOC0, _LOC1, _LOC2>::disp() const
		{
			char loc_str[2][50] = { "center\0", "half\0" };
			std::ostringstream displog;
			displog << "dimension:" << DIM::D
			        << "  location:[" << loc_str[_LOC0] << ", " << loc_str[_LOC1] << ", " << loc_str[_LOC2] << "]"
			        << "  type:[" << typeid(_NUMT).name() << "]"
			        << "  N = " << _N
			        << "  " << _size.disp();
			return displog.str();
		}

		template <typename _NUMT, unsigned _N, LO::LOCATION _LOC0, LO::LOCATION _LOC1, LO::LOCATION _LOC2>
		inline std::string PointData<_NUMT, _N, _LOC0, _LOC1, _LOC2>::disp_data() const
		{
			std::ostringstream displog;
			const int dln = 10;
			const int begin0 = _size.begin(0, RA::ALL);
			const int end0 = _size.end(0, RA::ALL) + _LOC0;
			const int begin1 = _size.begin(1, RA::ALL);
			const int end1 = _size.end(1, RA::ALL) + _LOC1;
			const int begin2 = _size.begin(2, RA::ALL);
			const int end2 = _size.end(2, RA::ALL) + _LOC2;
			for (unsigned n = 0; n != _N; ++n)
			{
				displog << "N = " << n << std::endl;
				for (int i = begin0; i != end0; ++i)
				{
					displog << "i = " << i;
					if(i > _size.last(0, RA::IN))
					{displog << " (+" << i - _size.last(0, RA::IN) << ")";}
					displog << std::endl;
					for(int ii = 0; ii != dln * _size.size(2, RA::ALL); ++ii) { displog << "="; } displog << std::endl;
					for (int j = begin1; j != end1; ++j)
					{
						for (int k = begin2; k != end2; ++k)
						{
							displog << pt0[n][i][j][k];
							if (k != _size.last(2, RA::N)
							        && k != _size.last(2, RA::IN))
							{displog << ", ";}
							else {displog << " | ";}
						}
						displog << std::endl;
						if (j == _size.last(1, RA::N)
						        || j == _size.last(1, RA::IN))
						{for(int ii = 0; ii != dln * _size.size(2, RA::ALL); ++ii) { displog << "-"; } displog << std::endl;}
					}
					for(int ii = 0; ii != dln * _size.size(2, RA::ALL); ++ii) { displog << "="; } displog << std::endl;
				}
				displog << std::endl;
			}
			return displog.str();
		}

		template <typename _NUMT, unsigned _N, LO::LOCATION _LOC0, LO::LOCATION _LOC1, LO::LOCATION _LOC2>
		void PointData<_NUMT, _N, _LOC0, _LOC1, _LOC2>::copy_from(const PointData<_NUMT, _N, _LOC0, _LOC1, _LOC2> &pd)
		{
#ifdef _CHECK_POINTDATA_RANGE
			if(&pd == this)
			{
				throw(std::runtime_error("copy from self"));
			}
#endif
			RangeT overlap = _size.range(RA::ALL, RA::ALL, RA::ALL).overlap(pd.size().range(RA::ALL, RA::ALL, RA::ALL));
			overlap.cut_tail(-static_cast<int>(_LOC0), -static_cast<int>(_LOC1), -static_cast<int>(_LOC2));
			PD_For_N_3D(0, N, overlap, [&]PD_F_n_ijk(n, i, j, k)
			{
				pt0[n][i][j][k] = pd.pt0[n][i][j][k];
			});
		}

		template <typename _NUMT, unsigned _N, LO::LOCATION _LOC0, LO::LOCATION _LOC1, LO::LOCATION _LOC2>
		template<typename __NUMT, unsigned __N>
		void PointData<_NUMT, _N, _LOC0, _LOC1, _LOC2>::copy_from(const PointData<__NUMT, __N, _LOC0, _LOC1, _LOC2> &pd,
		        const RangeT R, const RangeT r_,
		        unsigned N, unsigned n, int I, int J, int K, int i, int j, int k)
		{
			check_n(N);
			pd.check_n(n);
			RangeT r = r_;
			RangeT overlap = R.overlap(r.tran(I - i, J - j, K - k));
#ifdef _CHECK_POINTDATA_RANGE
			if(reinterpret_cast<const void*>(&pd) == reinterpret_cast<const void*>(this))
			{
				RangeT RR = overlap;
				RR.tran(i - I, j - J, k - K);
				RR.cut_tail(-static_cast<int>(_LOC0), -static_cast<int>(_LOC1), -static_cast<int>(_LOC2));
				RR.cut_head(-static_cast<int>(_LOC0), -static_cast<int>(_LOC1), -static_cast<int>(_LOC2));
				if(RR.is_overlap(overlap))
				{
					throw(std::runtime_error("copy and write zone overlap"));
				}
			}
#endif
			overlap.cut_tail(-static_cast<int>(_LOC0), -static_cast<int>(_LOC1), -static_cast<int>(_LOC2));
			int ddi = i - I, ddj = j - J, ddk = k - K;
			PD_For_3D(overlap, [&]PD_F_ijk(II, JJ, KK) {pt0[N][II][JJ][KK] = static_cast<_NUMT>(pd._pt0(n)[II + ddi][JJ + ddj][KK + ddk]);});
		}

#ifdef PD_OT
		template <typename _NUMT, unsigned _N, LO::LOCATION _LOC0, LO::LOCATION _LOC1, LO::LOCATION _LOC2>
		template<typename _R0, typename _R1, typename _R2, typename _FL0, typename _FL1, typename _FL2>
		void PointData<_NUMT, _N, _LOC0, _LOC1, _LOC2>::read_plt(const std::string &root, const liton_ot::TEC_FILE_LOG &teclog,
		        unsigned zone, const std::string &name,
		        _R0 R0, _R1 R1, _R2 R2, unsigned nn, _FL0 fl0, _FL1 fl1, _FL2 fl2)
		{
			check_n(nn);
			if(zone >= teclog.Zones.size())
			{
				throw(std::runtime_error("zone number is too big"));
			}

			unsigned N0_file = teclog.Zones[zone].Real_Max_C(DIM::D, 0);
			unsigned N1_file = teclog.Zones[zone].Real_Max_C(DIM::D, 1);
			unsigned N2_file = teclog.Zones[zone].Real_Max_C(DIM::D, 2);
			if(N0_file == 0 || N1_file == 0 || N2_file == 0)
			{
				throw(std::runtime_error("N_file is equal to 0"));
			}
			int vv = std::find(teclog.Variables.begin(), teclog.Variables.end(), name) - teclog.Variables.begin();
			if (vv == teclog.Variables.size())
			{
				throw(std::runtime_error("error when finding variables " + name + " in: " + teclog.FileName));
			}

			std::string plt_file_name = root + "/" + teclog.FileName + ".plt";
			std::ifstream data_in;
			data_in.open(plt_file_name.c_str(), std::ios::binary);
			data_in.seekg(teclog.Zones[zone].Data[vv].file_pt);
			if (!data_in)
			{
				throw(std::runtime_error("error when opening file: " + plt_file_name));
			}

			if (teclog.Zones[zone].Data[vv].type == 1)
			{
				PointData<float, 1, LO::center, LO::center, LO::center> buf(0, N0_file, 0, 0, N1_file, 0, 0, N2_file, 0);
				data_in.read(reinterpret_cast<char*>(buf.data_pt(0)), N0_file * N1_file * N2_file * teclog.Zones[zone].Data[vv].size);
				if (!data_in)
				{
					data_in.close();
					throw(std::runtime_error("error when reading file " + plt_file_name));
				}
				data_in.close();
				copy_from(buf, R0, R1, R2, RA::IN, RA::IN, RA::IN, nn, 0, fl0, fl1, fl2);
			}
			else if (teclog.Zones[zone].Data[vv].type == 2)
			{
				PointData<double, 1, LO::center, LO::center, LO::center> buf(0, N0_file, 0, 0, N1_file, 0, 0, N2_file, 0);
				data_in.read(reinterpret_cast<char*>(buf.data_pt(0)), N0_file * N1_file * N2_file * teclog.Zones[zone].Data[vv].size);
				if (!data_in)
				{
					data_in.close();
					throw(std::runtime_error("error when reading file " + plt_file_name));
				}
				data_in.close();
				copy_from(buf, R0, R1, R2, RA::IN, RA::IN, RA::IN, nn, 0, fl0, fl1, fl2);
			}
			else
			{
				throw(std::runtime_error("tec_file data type error"));
			}
		}
#endif
	}
}

#endif
