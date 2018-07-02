#ifndef POINTDATA_2D_HPP
#define POINTDATA_2D_HPP

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
	namespace D2
	{
		class DIM
		{
		  public:
			static const int D = 2;

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
			int _begin[DIM::D] = { 0, 0 };
			int _size[DIM::D] = { 0, 0 };

		  public:
			RangeT() = default;
			RangeT(const int b0, const unsigned s0,
			       const int b1, const unsigned s1):
				_begin{ b0, b1 },
				_size{ static_cast<int>(s0), static_cast<int>(s1) } {}

			inline const int* begin_pt() const { return _begin; }
			inline const int* size_pt() const { return _size; }

			inline int begin(const unsigned d) const { DIM::check_d(d); return _begin[d]; }
			inline int size(const unsigned d) const { DIM::check_d(d); return _size[d]; }
			inline int end(const unsigned d) const { DIM::check_d(d); return begin(d) + size(d); }
			inline int last(const unsigned d) const { DIM::check_d(d); return end(d) - 1; }
			inline int bound(const unsigned d, FL::_N fl) const { return begin(d); }
			inline int bound(const unsigned d, FL::_P fl) const { return last(d); }

			inline RangeT &cut_head(int i, int j)
			{
				_begin[0] += i; _size[0] -= i;
				_begin[1] += j; _size[1] -= j;
				return *this;
			}
			inline RangeT &cut_tail(int i, int j) { _size[0] -= i; _size[1] -= j; return *this; }
			inline RangeT &tran(int i, int j) { _begin[0] += i; _begin[1] += j; return *this; }

			inline bool is_overlap(const RangeT r) const
			{
				if(begin(0) >= r.end(0) || end(0) <= r.begin(0) ||
				        begin(1) >= r.end(1) || end(1) <= r.begin(1) )
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
				return RangeT(__begin0, __end0 - __begin0,
				              __begin1, __end1 - __begin1);
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
				        << "size[1] = " << size(1) ;
				return displog.str();
			}
		};
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		template <typename Function>
		inline void PD_For_2D(const RangeT r, const Function fun)
		{
			const int begin0 = r.begin(0);
			const int end0 = r.end(0);
			const int begin1 = r.begin(1);
			const int end1 = r.end(1);
			for (int i = begin0; i != end0; ++i)
			{
				for (int j = begin1; j != end1; ++j)
				{
					fun(i, j);
				}
			}
		}

		template <typename Function>
		inline void PD_For_N_2D(const unsigned N_b, const unsigned N_e, const RangeT &r, const Function fun)
		{
			const int begin0 = r.begin(0);
			const int end0 = r.end(0);
			const int begin1 = r.begin(1);
			const int end1 = r.end(1);
			for (unsigned n = N_b; n != N_e; ++n)
			{
				for (int i = begin0; i != end0; ++i)
				{
					for (int j = begin1; j != end1; ++j)
					{
						fun(n, i, j);
					}
				}
			}
		}

		template <typename T, typename Reducer, typename Function>
		inline void PD_Reduce_2D(const RangeT r, T &ans, const Reducer reduce, const Function fun)
		{
			T temp = ans;
			const int begin0 = r.begin(0);
			const int end0 = r.end(0);
			const int begin1 = r.begin(1);
			const int end1 = r.end(1);
			for (int i = begin0; i != end0; ++i)
			{
				for (int j = begin1; j != end1; ++j)
				{
					reduce(fun(i, j), temp);
				}
			}
			reduce(temp, ans);
		}

		template <typename T, typename Reducer, typename Function>
		inline void PD_Reduce_N_2D(const unsigned N_b, const unsigned N_e, const RangeT r, T &ans, const Reducer reduce, const Function fun)
		{
			T temp = ans;
			const int begin0 = r.begin(0);
			const int end0 = r.end(0);
			const int begin1 = r.begin(1);
			const int end1 = r.end(1);
			for (unsigned n = N_b; n != N_e; ++n)
			{
				for (int i = begin0; i != end0; ++i)
				{
					for (int j = begin1; j != end1; ++j)
					{
						reduce(fun(n, i, j), temp);
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
			int _in[DIM::D] = { 0, 0 };
			int _n[DIM::D] = { 0, 0 };
			int _p[DIM::D] = { 0, 0 };

		  public:
			SizeT() = default;
			SizeT(const unsigned in0, const unsigned iin0, const unsigned ip0,
			      const unsigned in1, const unsigned iin1, const unsigned ip1):
				_in{ static_cast<int>(iin0), static_cast<int>(iin1) },
				_n{ static_cast<int>(in0), static_cast<int>(in1) },
				_p{ static_cast<int>(ip0), static_cast<int>(ip1) } {}

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

			template<typename T0, typename T1>
			inline RangeT range(T0 r0, T1 r1) const
			{
				return RangeT(begin(0, r0), size(0, r0), begin(1, r1), size(1, r1));
			}

			std::string disp() const
			{
				std::ostringstream displog;
				displog << "size[0] = [" << _n[0] << " " << _in[0] << " " << _p[0] << "]" << "  ";
				displog << "size[1] = [" << _n[1] << " " << _in[1] << " " << _p[1] << "]";
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
			}
		};
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		template <typename _NUMT, unsigned _N, LO::LOCATION _LOC0, LO::LOCATION _LOC1>
		class PointData
		{
		  public:
			typedef _NUMT num_type;
			const LO::LOCATION LOC0 = _LOC0;
			const LO::LOCATION LOC1 = _LOC1;
			const unsigned N = _N;
		  protected:
			SizeT _size;
			_NUMT* data = nullptr;
			_NUMT** pt0[_N];
			_NUMT** pt1 = nullptr;

		  public:
			PointData();
			PointData(const unsigned in0, const unsigned iin0, const unsigned ip0,
			          const unsigned in1, const unsigned iin1, const unsigned ip1);
			PointData(const PointData<_NUMT, _N, _LOC0, _LOC1> &) = delete;
			const PointData<_NUMT, _N, _LOC0, _LOC1> &operator=(const PointData<_NUMT, _N, _LOC0, _LOC1> &) = delete;
			~PointData();

			void alloc(const unsigned in0, const unsigned iin0, const unsigned ip0,
			           const unsigned in1, const unsigned iin1, const unsigned ip1);
			void realloc(const unsigned in0, const unsigned iin0, const unsigned ip0,
			             const unsigned in1, const unsigned iin1, const unsigned ip1);
			void clear();

			std::string disp() const;
			std::string disp_data() const;

			inline SizeT size() const { return _size; }

			template<typename F0 = FL::_C, typename F1 = FL::_C>
			inline _NUMT & operator()(const unsigned n, const int i, const int j, const F0 flag0 = FL::C, const F1 flag1 = FL::C)
			{
				check_data();
				check_n(n);
				check_flag(flag0, flag1);
				_size.range(RA::ALL, RA::ALL).check_range(0, LOC0, i, F0::offset);
				_size.range(RA::ALL, RA::ALL).check_range(1, LOC1, j, F1::offset);
				return pt0[n][i + F0::offset][j + F1::offset];
			}
			template<typename F0 = FL::_C, typename F1 = FL::_C>
			inline const _NUMT & operator()(const unsigned n, const int i, const int j, const F0 flag0 = FL::C, const F1 flag1 = FL::C) const
			{
				check_data();
				check_n(n);
				check_flag(flag0, flag1);
				_size.range(RA::ALL, RA::ALL).check_range(0, LOC0, i, F0::offset);
				_size.range(RA::ALL, RA::ALL).check_range(1, LOC1, j, F1::offset);
				return pt0[n][i + F0::offset][j + F1::offset];
			}

			inline _NUMT* data_pt(const unsigned n) { check_n(n); return pt0[n][-_size.n(0)] - _size.n(1); }
			inline const _NUMT* data_pt(const unsigned n) const { check_n(n); return pt0[n][-_size.n(0)] - _size.n(1); }

			void copy_from(const PointData<_NUMT, _N, _LOC0, _LOC1> &pd,
			               const RangeT R, const RangeT r_,
			               int I, int J, int i, int j);
			template<typename _R0, typename _R1, typename _r0, typename _r1>
			inline void copy_from(const PointData<_NUMT, _N, _LOC0, _LOC1> &pd,
			                      _R0 R0, _R1 R1, _r0 r0, _r1 r1,
			                      int I, int J, int i, int j)
			{copy_from(pd, _size.range(R0, R1), pd.size().range(r0, r1), I, J, i, j);}
			template<typename _FL0, typename _FL1>
			inline void copy_from(const PointData<_NUMT, _N, _LOC0, _LOC1> &pd,
			                      const RangeT R, const RangeT r_,
			                      _FL0 fl0, _FL1 fl1)
			{copy_from(pd, R, r_, R.bound(0, fl0), R.bound(1, fl1), r_.bound(0, fl0), r_.bound(1, fl1));}
			template<typename _R0, typename _R1, typename _r0, typename _r1, typename _FL0, typename _FL1>
			inline void copy_from(const PointData<_NUMT, _N, _LOC0, _LOC1> &pd,
			                      _R0 R0, _R1 R1, _r0 r0, _r1 r1,
			                      _FL0 fl0, _FL1 fl1)
			{
				copy_from(pd,
				          _size.range(R0, R1), pd.size().range(r0, r1),
				          _size.bound(0, R0, fl0), _size.bound(1, R1, fl1),
				          pd.size().bound(0, r0, fl0), pd.size().bound(1, r1, fl1));
			}

#ifdef PD_OT
			void read_plt(const std::string &root, const liton_ot::TEC_FILE_LOG &teclog,
			              unsigned zone, const std::string &name, unsigned nn);
#endif

		  protected:
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

			template<typename F0, typename F1>
			inline void check_flag(const F0 flag0, const F1 flag1) const
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

		template <typename _NUMT, unsigned _N, LO::LOCATION _LOC0, LO::LOCATION _LOC1>
		PointData<_NUMT, _N, _LOC0, _LOC1>::PointData()
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

		template <typename _NUMT, unsigned _N, LO::LOCATION _LOC0, LO::LOCATION _LOC1>
		PointData<_NUMT, _N, _LOC0, _LOC1>::PointData(const unsigned in0, const unsigned iin0, const unsigned ip0,
		        const unsigned in1, const unsigned iin1, const unsigned ip1)
		{
			if(_N == 0)
			{
				throw(std::runtime_error("N is not allowed to be zero"));
			}
			for(unsigned n = 0; n != _N; ++n)
			{
				pt0[n] = nullptr;
			}
			alloc(in0, iin0, ip0, in1, iin1, ip1);
		}

		template <typename _NUMT, unsigned _N, LO::LOCATION _LOC0, LO::LOCATION _LOC1>
		PointData<_NUMT, _N, _LOC0, _LOC1>::~PointData()
		{
			if(data != nullptr)
			{
				clear();
			}
		}

		template <typename _NUMT, unsigned _N, LO::LOCATION _LOC0, LO::LOCATION _LOC1>
		void PointData<_NUMT, _N, _LOC0, _LOC1>::alloc(const unsigned in0, const unsigned iin0, const unsigned ip0,
		        const unsigned in1, const unsigned iin1, const unsigned ip1)
		{
			if(data != nullptr)
			{
				throw(std::runtime_error("allocating memory do not allowed when data is not empty"));
			}
			else
			{
				SizeT s(in0, iin0, ip0, in1, iin1, ip1);
				s.check();
				unsigned size_point = (s.size(0, RA::ALL) + _LOC0) * (s.size(1, RA::ALL) + _LOC1);
				if(size_point != 0)
				{
					try
					{
						data = new _NUMT[size_point * _N];
						_size = s;

						pt1 = new _NUMT*[N * (s.size(0, RA::ALL) + _LOC0)];
						for (unsigned n = 0; n != _N; ++n)
						{
							pt0[n] = pt1 + n * (s.size(0, RA::ALL) + _LOC0);
							for (unsigned ii = 0; ii != s.size(0, RA::ALL) + _LOC0; ++ii)
							{
								pt0[n][ii] = data + n * (s.size(0, RA::ALL) + _LOC0) * (s.size(1, RA::ALL) + _LOC1)
								             + ii * (s.size(1, RA::ALL) + _LOC1);
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

		template <typename _NUMT, unsigned _N, LO::LOCATION _LOC0, LO::LOCATION _LOC1>
		void PointData<_NUMT, _N, _LOC0, _LOC1>::realloc(const unsigned in0, const unsigned iin0, const unsigned ip0,
		        const unsigned in1, const unsigned iin1, const unsigned ip1)
		{
			SizeT(in0, iin0, ip0, in1, iin1, ip1).check();
			clear();
			alloc(in0, iin0, ip0, in1, iin1, ip1);
		}

		template <typename _NUMT, unsigned _N, LO::LOCATION _LOC0, LO::LOCATION _LOC1>
		void PointData<_NUMT, _N, _LOC0, _LOC1>::clear()
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
				_size = SizeT();
			}
		}

		template <typename _NUMT, unsigned _N, LO::LOCATION _LOC0, LO::LOCATION _LOC1>
		inline std::string PointData<_NUMT, _N, _LOC0, _LOC1>::disp() const
		{
			char loc_str[2][50] = { "center\0", "half\0" };
			std::ostringstream displog;
			displog << "dimension:" << DIM::D
			        << "  location:[" << loc_str[_LOC0] << ", " << loc_str[_LOC1] << "]"
			        << "  type:[" << typeid(_NUMT).name() << "]"
			        << "  N = " << _N
			        << "  " << _size.disp();
			return displog.str();
		}

		template <typename _NUMT, unsigned _N, LO::LOCATION _LOC0, LO::LOCATION _LOC1>
		inline std::string PointData<_NUMT, _N, _LOC0, _LOC1>::disp_data() const
		{
			std::ostringstream displog;
			const int begin0 = _size.begin(0, RA::ALL);
			const int end0 = _size.end(0, RA::ALL) + _LOC0;
			const int begin1 = _size.begin(1, RA::ALL);
			const int end1 = _size.end(1, RA::ALL) + _LOC1;;
			for (unsigned n = 0; n != _N; ++n)
			{
				displog << "N = " << n << std::endl;
				for (int i = begin0; i != end0; ++i)
				{
					for (int j = begin1; j != end1; ++j)
					{
						displog << pt0[n][i][j] << ", ";
					}
					displog << std::endl;
				}
				displog << std::endl;
			}
			return displog.str();
		}

		template <typename _NUMT, unsigned _N, LO::LOCATION _LOC0, LO::LOCATION _LOC1>
		void PointData<_NUMT, _N, _LOC0, _LOC1>::copy_from(const PointData<_NUMT, _N, _LOC0, _LOC1> &pd,
		        const RangeT R, const RangeT r_,
		        int I, int J, int i, int j)
		{
			RangeT r = r_;
			RangeT overlap = R.overlap(r.tran(I - i, J - j));
#ifdef _CHECK_POINTDATA_RANGE
			RangeT RR = R;
			RR.tran(i - I, j - J);
			RR.cut_tail(-static_cast<int>(_LOC0), -static_cast<int>(_LOC1));
			RR.cut_head(-static_cast<int>(_LOC0), -static_cast<int>(_LOC1));
			if(&pd == this)
			{
				if(RR.is_overlap(overlap))
				{
					throw(std::runtime_error("copy and write zone overlap"));
				}
			}
#endif
			overlap.cut_tail(-static_cast<int>(_LOC0), -static_cast<int>(_LOC1));
			int ddi = i - I, ddj = j - J;
			int begin0 = overlap.begin(0);
			int end0 = overlap.end(0);
			int begin1 = overlap.begin(1);
			int end1 = overlap.end(1);

			for (int II = begin0; II != end0; ++II)
			{
				for (int JJ = begin1; JJ != end1; ++JJ)
				{
					for (unsigned n = 0; n != N; ++n)
					{
						pt0[n][II][JJ] = pd.pt0[n][II + ddi][JJ + ddj];
					}
				}
			}
		}

#ifdef PD_OT
		template <typename _NUMT, unsigned _N, LO::LOCATION _LOC0, LO::LOCATION _LOC1>
		void PointData<_NUMT, _N, _LOC0, _LOC1>::read_plt(const std::string &root, const liton_ot::TEC_FILE_LOG &teclog,
		        unsigned zone, const std::string &name, unsigned nn)
		{
			check_n(nn);
			if(zone >= teclog.Zones.size())
			{
				throw(std::runtime_error("zone number is too big"));
			}

			unsigned N0_file = teclog.Zones[zone].Real_Max_C(DIM::D, 0);
			unsigned N1_file = teclog.Zones[zone].Real_Max_C(DIM::D, 1);
			if(N0_file > static_cast<unsigned>(_size._in[0]) || N0_file == 0
			        || N1_file > static_cast<unsigned>(_size._in[1]) || N1_file == 0)
			{
				throw(std::runtime_error("N_file is bigger than N_in or equal to 0"));
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
			char* buf = nullptr;
			try
			{
				buf = new char[N0_file * N1_file * teclog.Zones[zone].Data[vv].size];
			}
			catch (const std::bad_alloc &)
			{
				data_in.close();
				throw(std::runtime_error("no more memory for allocating"));
			}
			data_in.read(buf, N0_file * N1_file * teclog.Zones[zone].Data[vv].size);
			if (!data_in)
			{
				data_in.close();
				delete []buf;
				throw(std::runtime_error("error when reading file " + plt_file_name));
			}
			data_in.close();

			if (teclog.Zones[zone].Data[vv].type == 1)
			{
				float* temp = reinterpret_cast<float*>(buf);
				for(unsigned ii = 0; ii != N0_file; ++ii)
				{
					for(unsigned jj = 0; jj != N1_file; ++jj)
					{
						pt0[nn][ii][jj] = temp[ii * N1_file + jj];
					}
				}
			}
			else if (teclog.Zones[zone].Data[vv].type == 2)
			{
				double* temp = reinterpret_cast<double*>(buf);
				for(unsigned ii = 0; ii != N0_file; ++ii)
				{
					for(unsigned jj = 0; jj != N1_file; ++jj)
					{
						pt0[nn][ii][jj] = temp[ii * N1_file + jj];
					}
				}
			}
			else
			{
				delete []buf;
				throw(std::runtime_error("tec_file data type error"));
			}
			delete []buf;
		}
#endif
	}
}

#endif
