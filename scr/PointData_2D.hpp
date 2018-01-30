#ifndef POINTDATA_2D_HPP
#define POINTDATA_2D_HPP

#include <stdexcept>
#include <typeinfo>
#include <string>
#include <sstream>

namespace liton_pd
{
	namespace D2
	{
		static const int DIM = 2;

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

		class RangeT
		{
		  protected:
			int _begin[DIM] = { 0, 0 };
			int _end[DIM] = { 0, 0 };
			int _size[DIM] = { 0, 0 };

		  public:
			RangeT() = default;
			RangeT(const int &b0, const int &b1,
			       const unsigned &s0, const unsigned &s1):
				_begin{ b0, b1 },
				_size{ static_cast<int>(s0), static_cast<int>(s1) }
			{
				_end[0] = _begin[0] + _size[0];
				_end[1] = _begin[1] + _size[1];
			}

			inline int begin(const unsigned &d) const { check_d(d); return _begin[d]; }
			inline int end(const unsigned &d) const { check_d(d); return _end[d]; }
			inline int size(const unsigned &d) const { check_d(d); return _size[d]; }

			inline const int* begin_pt() const { return _begin; }
			inline const int* end_pt() const { return _end; }
			inline const int* size_pt() const { return _size; }

			inline std::string disp() const
			{
				std::ostringstream displog;
				displog << "range_0 = [" << _begin[0] << " , " << _end[0] << "]  "
				        << "size_0 = " << _size[0] << "    ";
				displog << "range_1 = [" << _begin[1] << " , " << _end[1] << "]  "
				        << "size_1 = " << _size[1] ;
				return displog.str();
			}
		};

		class SizeT
		{
		  protected:
			int _in[DIM] = { 0, 0 };
			int _n[DIM] = { 0, 0 };
			int _p[DIM] = { 0, 0 };

		  public:
			SizeT() = default;
			SizeT(const unsigned &iin0, const unsigned &iin1,
			      const unsigned &in0, const unsigned &in1,
			      const unsigned &ip0, const unsigned &ip1):
				_in{ static_cast<int>(iin0), static_cast<int>(iin1) },
				_n{ static_cast<int>(in0), static_cast<int>(in1) },
				_p{ static_cast<int>(ip0), static_cast<int>(ip1) } {}

			inline int in(const unsigned &d) const { check_d(d); return _in[d]; }
			inline int n(const unsigned &d) const { check_d(d); return _n[d]; }
			inline int p(const unsigned &d) const { check_d(d); return _p[d]; }

			inline int sum(const unsigned &d) const { check_d(d); return _in[d] + _n[d] + _p[d]; }

			inline int begin(const unsigned &d, RA::_ALL r) const { check_d(d); return -static_cast<int>(_n[d]); }
			inline int end(const unsigned &d, RA::_ALL r) const { check_d(d); return static_cast<int>(_in[d] + _p[d]); }
			inline int size(const unsigned &d, RA::_ALL r) const { check_d(d); return static_cast<int>(_in[d] + _n[d] + _p[d]); }
			inline int begin(const unsigned &d, RA::_IN r) const { check_d(d); return 0; }
			inline int end(const unsigned &d, RA::_IN r) const { check_d(d); return static_cast<int>(_in[d]); }
			inline int size(const unsigned &d, RA::_IN r) const { check_d(d); return static_cast<int>(_in[d]); }
			inline int begin(const unsigned &d, RA::_N r) const { check_d(d); return -static_cast<int>(_n[d]); }
			inline int end(const unsigned &d, RA::_N r) const { check_d(d); return 0; }
			inline int size(const unsigned &d, RA::_N r) const { check_d(d); return static_cast<int>(_n[d]); }
			inline int begin(const unsigned &d, RA::_P r) const { check_d(d); return static_cast<int>(_in[d]); }
			inline int end(const unsigned &d, RA::_P r) const { check_d(d); return static_cast<int>(_in[d] + _p[d]); }
			inline int size(const unsigned &d, RA::_P r) const { check_d(d); return static_cast<int>(_p[d]); }

			template<typename T0>
			inline int last(const unsigned &d, T0 r) const { return end(d, r) - 1; }
			template<typename T0, typename T1>
			inline RangeT range(T0 r0, T1 r1) const { return RangeT(begin(0, r0), begin(1, r1), size(0, r0), size(1, r1)); }

			std::string disp() const
			{
				std::ostringstream displog;
				displog << "size_0 = [" << _n[0] << " " << _in[0] << " " << _p[0] << "]" << "  ";
				displog << "size_1 = [" << _n[1] << " " << _in[1] << " " << _p[1] << "]";
				return displog.str();
			}

			inline void check() const
			{
				if (_in[0] == 0 && _n[0] + _p[0] != 0)
				{
					throw(std::runtime_error("dim 0: size_in can not be zero when size_n or size_p is non-zero"));
				}
				if (_in[1] == 0 && _n[1] + _p[1] != 0)
				{
					throw(std::runtime_error("dim 1: size_in can not be zero when size_n or size_p is non-zero"));
				}
			}

			inline void check_range(const int &i, const int &j) const
			{
#ifdef _CHECK_POINTDATA_RANGE
				if (i < -static_cast<int>(_n[0]) || i >= static_cast<int>(_in[0] + _p[0]))
				{
					std::ostringstream errlog;
					errlog << "out of sub range: i:[" << i << "] range:[" << -static_cast<int>(_n[0]) << "," << static_cast<int>(_in[0] + _p[0]) - 1 << "]";
					throw(std::runtime_error(errlog.str()));
				}
				if (j < -static_cast<int>(_n[1]) || j >= static_cast<int>(_in[1] + _p[1]))
				{
					std::ostringstream errlog;
					errlog << "out of sub range: j:[" << j << "] range:[" << -static_cast<int>(_n[1]) << "," << static_cast<int>(_in[1] + _p[1]) - 1 << "]";
					throw(std::runtime_error(errlog.str()));
				}
#endif
			}
		};

		template <typename _NUMT, LO::LOCATION _LOC0, LO::LOCATION _LOC1, unsigned _N>
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

		  public:
			PointData();
			PointData(const unsigned &iin0, const unsigned &iin1,
			          const unsigned &in0, const unsigned &in1,
			          const unsigned &ip0, const unsigned &ip1);
			PointData(const PointData<_NUMT, _LOC0, _LOC1, _N> &) = delete;
			const PointData<_NUMT, _LOC0, _LOC1, _N> &operator=(const PointData<_NUMT, _LOC0, _LOC1, _N> &) = delete;
			~PointData();

			void alloc(const unsigned &iin0, const unsigned &iin1,
			           const unsigned &in0, const unsigned &in1,
			           const unsigned &ip0, const unsigned &ip1);
			void realloc(const unsigned &iin0, const unsigned &iin1,
			             const unsigned &in0, const unsigned &in1,
			             const unsigned &ip0, const unsigned &ip1);
			void clear();

			std::string disp() const;
			std::string disp_data() const;

			inline SizeT size() const { return _size; }

			template<typename F0, typename F1>
			inline _NUMT &operator()(const int &i, const int &j, const unsigned &n, const F0 &flag0, const F1 &flag1)
			{
				check_n(n);
				check_flag(flag0, flag1);
				_size.check_range(i, j);
				return pt0[n][i + _size.n(0) + F0::offset][j + F1::offset];
			}
			template<typename F0, typename F1>
			inline const _NUMT &operator()(const int &i, const int &j, const unsigned &n, const F0 &flag0, const F1 &flag1) const { return *this(i, j, n, flag0, flag1); }

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

			template<typename F0, typename F1>
			inline void check_flag(const F0 &flag0, const F1 &flag1)
			{
#ifdef _CHECK_POINTDATA_RANGE
				if (_LOC0 == LO::center && typeid(F0) != typeid(FL::_C))
				{
					throw(std::runtime_error("dim 0: flag must be [C] when location is [center]"));
				}
				if (_LOC0 == LO::half && typeid(F0) == typeid(FL::_C))
				{
					throw(std::runtime_error("dim 0: flag must be [N] or [P] when location is [half]"));
				}
				if (_LOC1 == LO::center && typeid(F1) != typeid(FL::_C))
				{
					throw(std::runtime_error("dim 1: flag must be [C] when location is [center]"));
				}
				if (_LOC1 == LO::half && typeid(F1) == typeid(FL::_C))
				{
					throw(std::runtime_error("dim 1: flag must be [N] or [P] when location is [half]"));
				}
#endif
			}
		};

		template <typename _NUMT, LO::LOCATION _LOC0, LO::LOCATION _LOC1, unsigned _N>
		PointData<_NUMT, _LOC0, _LOC1, _N>::PointData()
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

		template <typename _NUMT, LO::LOCATION _LOC0, LO::LOCATION _LOC1, unsigned _N>
		PointData<_NUMT, _LOC0, _LOC1, _N>::PointData(const unsigned &iin0, const unsigned &iin1,
		        const unsigned &in0, const unsigned &in1,
		        const unsigned &ip0, const unsigned &ip1)
		{
			if(_N == 0)
			{
				throw(std::runtime_error("N is not allowed to be zero"));
			}
			for(unsigned n = 0; n != _N; ++n)
			{
				pt0[n] = nullptr;
			}
			alloc(iin0, iin1, in0, in1, ip0, ip1);
		}

		template <typename _NUMT, LO::LOCATION _LOC0, LO::LOCATION _LOC1, unsigned _N>
		PointData<_NUMT, _LOC0, _LOC1, _N>::~PointData()
		{
			if(data != nullptr)
			{
				clear();
			}
		}

		template <typename _NUMT, LO::LOCATION _LOC0, LO::LOCATION _LOC1, unsigned _N>
		void PointData<_NUMT, _LOC0, _LOC1, _N>::alloc(const unsigned &iin0, const unsigned &iin1,
		        const unsigned &in0, const unsigned &in1,
		        const unsigned &ip0, const unsigned &ip1)
		{
			if(data != nullptr)
			{
				throw(std::runtime_error("allocating memory do not allowed when data is not empty"));
			}
			else
			{
				SizeT s(iin0, iin1, in0, in1, ip0, ip1);
				s.check();
				unsigned size_point = (s.sum(0) + _LOC0) * (s.sum(1) + _LOC1);
				if(size_point != 0)
				{
					try
					{
						data = new _NUMT[size_point * _N];
						_size = s;
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
				for(unsigned n = 0; n != _N; ++n)
				{
					pt0[n] = new _NUMT*[s.sum(0) + _LOC0];
					for (unsigned ii = 0; ii != s.sum(0) + _LOC0; ++ii)
					{
						pt0[n][ii] = data + n * size_point + ii * (s.sum(1) + _LOC1) + s.n(1);
					}
				}
			}
		}

		template <typename _NUMT, LO::LOCATION _LOC0, LO::LOCATION _LOC1, unsigned _N>
		void PointData<_NUMT, _LOC0, _LOC1, _N>::realloc(const unsigned &iin0, const unsigned &iin1,
		        const unsigned &in0, const unsigned &in1,
		        const unsigned &ip0, const unsigned &ip1)
		{
			SizeT(iin0, iin1, in0, in1, ip0, ip1).check();
			clear();
			alloc(iin0, iin1, in0, in1, ip0, ip1);
		}

		template <typename _NUMT, LO::LOCATION _LOC0, LO::LOCATION _LOC1, unsigned _N>
		void PointData<_NUMT, _LOC0, _LOC1, _N>::clear()
		{
			if(data == nullptr)
			{
				throw(std::runtime_error("clearing data do not allowed when data is empty"));
			}
			else
			{
				delete[] data;
				data = nullptr;
				for(unsigned n = 0; n != _N; ++n)
				{
					delete[] pt0[n];
					pt0[n] = nullptr;
				}
				_size = SizeT();
			}
		}

		template <typename _NUMT, LO::LOCATION _LOC0, LO::LOCATION _LOC1, unsigned _N>
		inline std::string PointData<_NUMT, _LOC0, _LOC1, _N>::disp() const
		{
			char loc_str[2][50] = { "center\0", "half\0" };
			std::ostringstream displog;
			displog << "dimension:" << DIM
			        << "  location:[" << loc_str[_LOC0] << ", " << loc_str[_LOC1] << "]"
			        << "  type:[" << typeid(_NUMT).name() << "]"
			        << "  N = " << _N
			        << "  " << _size.disp();
			return displog.str();
		}

		template <typename _NUMT, LO::LOCATION _LOC0, LO::LOCATION _LOC1, unsigned _N>
		inline std::string PointData<_NUMT, _LOC0, _LOC1, _N>::disp_data() const
		{
			std::ostringstream displog;
			const int begin0 = _size.begin(0, RA::ALL);
			const int end0 = _size.end(0, RA::ALL) + _LOC0;
			const int begin1 = _size.begin(1, RA::ALL);
			const int end1 = _size.end(1, RA::ALL) + _LOC1;;
			for (unsigned n = 0; n != _N; ++n)
			{
				for (int i = begin0; i != end0; ++i)
				{
					for (int j = begin1; j != end1; ++j)
					{
						displog << pt0[n][i + _size.n(0)][j] << ", ";
					}
					displog << std::endl;
				}
				displog << std::endl;
			}
			return displog.str();
		}

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		template <typename Function>
		inline void For_PD_2D(const RangeT &r, const Function &fun)
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
		inline void For_PD_2D_N(const RangeT &r, const unsigned &N_b, const unsigned &N_e, const Function &fun)
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
						fun(i, j, n);
					}
				}
			}
		}

		template <typename T, typename Reducer, typename Function>
		inline void Reduce_PD_2D(const RangeT &r, T &ans, const Reducer &reduce, const Function &fun)
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
		inline void Reduce_PD_2D_N(const RangeT &r, const unsigned &N_b, const unsigned &N_e, T &ans, const Reducer &reduce, const Function &fun)
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
						reduce(fun(i, j, n), temp);
					}
				}
			}
			reduce(temp, ans);
		}

		template <typename Function>
		inline void For_PD_1D(const RangeT &r, const unsigned &d0, const Function &fun)
		{
			const int begin0 = r.begin(d0);
			const int end0 = r.end(d0);
			for (int ii = begin0; ii != end0; ++ii)
			{
				fun(ii);
			}
		}

		template <typename Function>
		inline void For_PD_1D_N(const RangeT &r, const unsigned &d0, const unsigned &N_b, const unsigned &N_e, const Function &fun)
		{
			const int begin0 = r.begin(d0);
			const int end0 = r.end(d0);
			for (unsigned n = N_b; n != N_e; ++n)
			{
				for (int ii = begin0; ii != end0; ++ii)
				{
					fun(ii, n);
				}
			}
		}

		template <typename T, typename Reducer, typename Function>
		inline void Reduce_PD_1D(const RangeT &r, const unsigned &d0, T &ans, const Reducer &reduce, const Function &fun)
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
		inline void Reduce_PD_1D_N(const RangeT &r, const unsigned &d0, const unsigned &N_b, const unsigned &N_e, T &ans, const Reducer &reduce, const Function &fun)
		{
			T temp = ans;
			const int begin0 = r.begin(d0);
			const int end0 = r.end(d0);
			for (unsigned n = N_b; n != N_e; ++n)
			{
				for (int ii = begin0; ii != end0; ++ii)
				{
					reduce(fun(ii, n), temp);
				}
			}
			reduce(temp, ans);
		}

	}
}

#endif
