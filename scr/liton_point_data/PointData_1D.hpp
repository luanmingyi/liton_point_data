#ifndef POINTDATA_1D_HPP
#define POINTDATA_1D_HPP

#include <stdexcept>
#include <typeinfo>
#include <string>
#include <sstream>

namespace liton_pd
{
	namespace D1
	{
		static const int DIM = 1;

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
			int _begin = 0;
			int _end = 0;
			int _size = 0;

		  public:
			RangeT() = default;
			RangeT(const int &b, const unsigned &s): _begin(b), _size(static_cast<int>(s)) { _end = _begin + _size; }

			inline int begin(const unsigned &d) const { check_d(d); return _begin; }
			inline int end(const unsigned &d) const { check_d(d); return _end; }
			inline int size(const unsigned &d) const { check_d(d); return _size; }

			inline const int* begin_pt() const { return &_begin; }
			inline const int* end_pt() const { return &_end; }
			inline const int* size_pt() const { return &_size; }

			inline std::string disp() const
			{
				std::ostringstream displog;
				displog << "range_0 = [" << _begin << " , " << _end << "]  "
				        << "size_0 = " << _size;
				return displog.str();
			}
		};

		class SizeT
		{
		  protected:
			int _in = 0;
			int _n = 0;
			int _p = 0;

		  public:
			SizeT() = default;
			SizeT(const unsigned &in, const unsigned &iin, const unsigned &ip)
				: _in(static_cast<int>(iin)), _n(static_cast<int>(in)), _p(static_cast<int>(ip)) {}

			inline int in(const unsigned &d) const { check_d(d); return _in; }
			inline int n(const unsigned &d) const { check_d(d); return _n; }
			inline int p(const unsigned &d) const { check_d(d); return _p; }

			inline int sum(const unsigned &d) const { check_d(d); return _in + _n + _p; }

			inline int begin(const unsigned &d, RA::_ALL r) const { check_d(d); return -_n; }
			inline int end(const unsigned &d, RA::_ALL r) const { check_d(d); return _in + _p; }
			inline int size(const unsigned &d, RA::_ALL r) const { check_d(d); return _in + _n + _p; }
			inline int begin(const unsigned &d, RA::_IN r) const { check_d(d); return 0; }
			inline int end(const unsigned &d, RA::_IN r) const { check_d(d); return _in; }
			inline int size(const unsigned &d, RA::_IN r) const { check_d(d); return _in; }
			inline int begin(const unsigned &d, RA::_N r) const { check_d(d); return -_n; }
			inline int end(const unsigned &d, RA::_N r) const { check_d(d); return 0; }
			inline int size(const unsigned &d, RA::_N r) const { check_d(d); return _n; }
			inline int begin(const unsigned &d, RA::_P r) const { check_d(d); return _in; }
			inline int end(const unsigned &d, RA::_P r) const { check_d(d); return _in + _p; }
			inline int size(const unsigned &d, RA::_P r) const { check_d(d); return _p; }

			template<typename T0>
			inline int last(const unsigned &d, T0 r) const { return end(d, r) - 1; }
			inline int mirror(const unsigned &d, FL::_N fl, int ii) const { return 2 * begin(d, RA::IN) - ii; }
			inline int mirror(const unsigned &d, FL::_P fl, int ii) const { return 2 * last(d, RA::IN) - ii; }

			template<typename T0>
			inline RangeT range(T0 r) const { return RangeT(begin(0, r), size(0, r)); }

			std::string disp() const
			{
				std::ostringstream displog;
				displog << "size_0 = [" << _n << " " << _in << " " << _p << "]";
				return displog.str();
			}

			inline void check() const
			{
				if(_in == 0 && _n + _p != 0)
				{
					throw(std::runtime_error("dim 0: size_in can not be zero when size_n or size_p is non-zero"));
				}
			}

			inline void check_range(const int &i) const
			{
#ifdef _CHECK_POINTDATA_RANGE
				if(i < -_n || i >= _in + _p)
				{
					std::ostringstream errlog;
					errlog << "out of sub range: i:[" << i << "] range:[" << -_n << "," << _in + _p - 1 << "]";
					throw(std::runtime_error(errlog.str()));
				}
#endif
			}
		};

		template <typename _NUMT, unsigned _N, LO::LOCATION _LOC0>
		class PointData
		{
		  public:
			typedef _NUMT num_type;
			const LO::LOCATION LOC0 = _LOC0;
			const unsigned N = _N;
		  protected:
			SizeT _size;
			_NUMT* data = nullptr;
			_NUMT* pt0[_N];

		  public:
			PointData();
			PointData(const unsigned &in, const unsigned &iin, const unsigned &ip);
			PointData(const PointData<_NUMT, _N, _LOC0> &) = delete;
			const PointData<_NUMT, _N, _LOC0> &operator=(const PointData<_NUMT, _N, _LOC0> &) = delete;
			~PointData();

			void alloc(const unsigned &in, const unsigned &iin, const unsigned &ip);
			void realloc(const unsigned &in, const unsigned &iin, const unsigned &ip);
			void clear();

			std::string disp() const;
			std::string disp_data() const;

			inline const SizeT & size() const { return _size; }

			template<typename F0>
			inline _NUMT &operator()(const unsigned &n, const int &i, const F0 &flag0)
			{
				check_n(n);
				check_flag(flag0);
				_size.check_range(i);
				return pt0[n][i + F0::offset];
			}
			template<typename F0>
			inline const _NUMT &operator()(const unsigned &n, const int &i, const F0 &flag0) const { return *this(i, n, flag0); }

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

			template<typename F0>
			inline void check_flag(const F0 &flag0)
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
#endif
			}
		};

		template <typename _NUMT, unsigned _N, LO::LOCATION _LOC0>
		PointData<_NUMT, _N, _LOC0>::PointData()
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

		template <typename _NUMT, unsigned _N, LO::LOCATION _LOC0>
		PointData<_NUMT, _N, _LOC0>::PointData(const unsigned &in, const unsigned &iin, const unsigned &ip)
		{
			if(_N == 0)
			{
				throw(std::runtime_error("N is not allowed to be zero"));
			}
			for(unsigned n = 0; n != _N; ++n)
			{
				pt0[n] = nullptr;
			}
			alloc(in, iin, ip);
		}

		template <typename _NUMT, unsigned _N, LO::LOCATION _LOC0>
		PointData<_NUMT, _N, _LOC0>::~PointData()
		{
			if(data != nullptr)
			{
				clear();
			}
		}

		template<typename _NUMT, unsigned _N, LO::LOCATION _LOC0>
		void PointData<_NUMT, _N, _LOC0>::alloc(const unsigned &in, const unsigned &iin, const unsigned &ip)
		{
			if(data != nullptr)
			{
				throw(std::runtime_error("allocating memory do not allowed when data is not empty"));
			}
			else
			{
				SizeT s(in, iin, ip);
				s.check();
				unsigned size_point = (s.sum(0) + _LOC0);
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
					pt0[n] = data + n * size_point + s.n(0);
				}
			}
		}

		template<typename _NUMT, unsigned _N, LO::LOCATION _LOC0>
		void PointData<_NUMT, _N, _LOC0>::realloc(const unsigned &in, const unsigned &iin, const unsigned &ip)
		{
			SizeT(in, iin, ip).check();
			clear();
			alloc(in, iin, ip);
		}

		template<typename _NUMT, unsigned _N, LO::LOCATION _LOC0>
		void PointData<_NUMT, _N, _LOC0>::clear()
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
					pt0[n] = nullptr;
				}
				_size = SizeT();
			}
		}

		template<typename _NUMT, unsigned _N, LO::LOCATION _LOC0>
		inline std::string PointData<_NUMT, _N, _LOC0>::disp() const
		{
			char loc_str[2][50] = { "center\0", "half\0" };
			std::ostringstream displog;
			displog << "dimension:" << DIM
			        << "  location:[" << loc_str[_LOC0] << "]"
			        << "  type:[" << typeid(_NUMT).name() << "]"
			        << "  N = " << _N
			        << "  " << _size.disp();
			return displog.str();
		}

		template<typename _NUMT, unsigned _N, LO::LOCATION _LOC0>
		inline std::string PointData<_NUMT, _N, _LOC0>::disp_data() const
		{
			std::ostringstream displog;
			int begin0 = _size.begin(0, RA::ALL);
			int end0 = _size.end(0, RA::ALL) + _LOC0;
			for (int i = begin0; i != end0; ++i)
			{
				for (unsigned n = 0; n != N; ++n)
				{
					displog << pt0[n][i] << ", ";
				}
				displog << std::endl;
			}
			return displog.str();
		}

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		template <typename Function>
		inline void PD_For_1D(const RangeT &r, const Function &fun)
		{
			const int begin = r.begin(0);
			const int end = r.end(0);
			for (int i = begin; i != end; ++i)
			{
				fun(i);
			}
		}

		template <typename Function>
		inline void PD_For_N_1D(const unsigned &N_b, const unsigned &N_e, const RangeT &r, const Function &fun)
		{
			const int begin = r.begin(0);
			const int end = r.end(0);
			for (unsigned n = N_b; n != N_e; ++n)
			{
				for (int i = begin; i != end; ++i)
				{
					fun(n, i);
				}
			}
		}

		template <typename T, typename Reducer, typename Function>
		inline void PD_Reduce_1D(const RangeT &r, T &ans, const Reducer &reduce, const Function &fun)
		{
			T temp = ans;
			int begin = r.begin(0);
			int end = r.end(0);
			for (int i = begin; i != end; ++i)
			{
				reduce(fun(i), temp);
			}
			reduce(temp, ans);
		}

		template <typename T, typename Reducer, typename Function>
		inline void PD_Reduce_N_1D(const unsigned &N_b, const unsigned &N_e, const RangeT &r, T &ans, const Reducer &reduce, const Function &fun)
		{
			T temp = ans;
			int begin = r.begin(0);
			int end = r.end(0);
			for (unsigned n = N_b; n != N_e; ++n)
			{
				for (int i = begin; i != end; ++i)
				{
					reduce(fun(n, i), temp);
				}
			}
			reduce(temp, ans);
		}
	}
}

#endif