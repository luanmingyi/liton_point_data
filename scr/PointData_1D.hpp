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
			SizeT(const unsigned &iin, const unsigned &in, const unsigned &ip)
				: _in(static_cast<int>(iin)), _n(static_cast<int>(in)), _p(static_cast<int>(ip)) {}

			inline int in(const unsigned &d) const { check_d(d); return _in; }
			inline int n(const unsigned &d) const { check_d(d); return _n; }
			inline int p(const unsigned &d) const { check_d(d); return _p; }

			inline int sum(const unsigned &d) const { check_d(d); return _in + _n + _p; }

			inline int begin(const unsigned &d, RA::_ALL r) const { check_d(d); return -static_cast<int>(_n); }
			inline int end(const unsigned &d, RA::_ALL r) const { check_d(d); return static_cast<int>(_in + _p); }
			inline int size(const unsigned &d, RA::_ALL r) const { check_d(d); return static_cast<int>(_in + _n + _p); }
			inline int begin(const unsigned &d, RA::_IN r) const { check_d(d); return 0; }
			inline int end(const unsigned &d, RA::_IN r) const { check_d(d); return static_cast<int>(_in); }
			inline int size(const unsigned &d, RA::_IN r) const { check_d(d); return static_cast<int>(_in); }
			inline int begin(const unsigned &d, RA::_N r) const { check_d(d); return -static_cast<int>(_n); }
			inline int end(const unsigned &d, RA::_N r) const { check_d(d); return 0; }
			inline int size(const unsigned &d, RA::_N r) const { check_d(d); return static_cast<int>(_n); }
			inline int begin(const unsigned &d, RA::_P r) const { check_d(d); return static_cast<int>(_in); }
			inline int end(const unsigned &d, RA::_P r) const { check_d(d); return static_cast<int>(_in + _p); }
			inline int size(const unsigned &d, RA::_P r) const { check_d(d); return static_cast<int>(_p); }

			template<typename T0>
			inline int last(const unsigned &d, T0 r) const { return end(d, r) - 1; }
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
				if(i < -static_cast<int>(_n) || i >= static_cast<int>(_in + _p))
				{
					std::ostringstream errlog;
					errlog << "out of sub range: i:[" << i << "] range:[" << -static_cast<int>(_n) << "," << static_cast<int>(_in + _p) - 1 << "]";
					throw(std::runtime_error(errlog.str()));
				}
#endif
			}
		};

		template <typename _NUMT, LO::LOCATION _LOC0, unsigned _N>
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
			PointData(const unsigned &iin, const unsigned &in, const unsigned &ip);
			PointData(const PointData<_NUMT, _LOC0, _N> &) = delete;
			const PointData<_NUMT, _LOC0, _N> &operator=(const PointData<_NUMT, _LOC0, _N> &) = delete;
			~PointData();

			void alloc(const unsigned &iin, const unsigned &in, const unsigned &ip);
			void realloc(const unsigned &iin, const unsigned &in, const unsigned &ip);
			void clear();

			std::string disp() const;
			std::string disp_data() const;

			inline SizeT size() const { return _size; }

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
			using PointData<_NUMT, LO::center, _N>::LOC0;
			using PointData<_NUMT, LO::center, _N>::N;
		  protected:
			using PointData<_NUMT, LO::center, _N>::_size;
			using PointData<_NUMT, LO::center, _N>::pt0;
		  public:
			PointData_C(): PointData<_NUMT, LO::center, _N>() {}
			PointData_C(const unsigned &iin, const unsigned &in,
			            const unsigned &ip): PointData<_NUMT, LO::center, _N>(iin, in, ip) {}

			inline _NUMT &operator()(const int &i, const unsigned &n)
			{
				this->check_n(n);
				_size.check_range(i);
				return pt0[n][i];
			}
			inline const _NUMT &operator()(const int &i, const unsigned &n) const { return *this(i, n); }
		};

		template <typename _NUMT, unsigned _N = 1>
		class PointData_H : public PointData<_NUMT, LO::half, _N>
		{
		  public:
			using PointData<_NUMT, LO::half, _N>::num_type;
			using PointData<_NUMT, LO::half, _N>::LOC0;
			using PointData<_NUMT, LO::half, _N>::N;
		  protected:
			using PointData<_NUMT, LO::half, _N>::_size;
			using PointData<_NUMT, LO::half, _N>::pt0;
		  public:
			PointData_H(): PointData<_NUMT, LO::half, _N>() {}
			PointData_H(const unsigned &iin, const unsigned &in,
			            const unsigned &ip): PointData<_NUMT, LO::half, _N>(iin, in, ip) {}

			inline _NUMT &operator()(const int &i, const unsigned &n, const FL::FLAG &flag0)
			{
				this->check_n(n);
				_size.check_range(i);
				return pt0[n][i + flag0];
			}
			inline const _NUMT &operator()(const int &i, const unsigned &n, const FL::FLAG &flag0) const { return *this(i, n, flag0); }
		};

		template <typename Function>
		inline void For_PD_1D(const RangeT &r, const Function &fun);
		template <typename Function>
		inline void For_PD_1D_N(const RangeT &r, const unsigned &N_b, const unsigned &N_e, const Function &fun);
		template <typename T, typename Reducer, typename Function>
		inline void Reduce_PD_1D(const RangeT &r, T &ans, const Reducer &reduce, const Function &fun);
		template <typename T, typename Reducer, typename Function>
		inline void Reduce_PD_1D_N(const RangeT &r, const unsigned &N_b, const unsigned &N_e, T &ans, const Reducer &reduce, const Function &fun);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		template <typename _NUMT, LO::LOCATION _LOC0, unsigned _N>
		PointData<_NUMT, _LOC0, _N>::PointData()
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

		template <typename _NUMT, LO::LOCATION _LOC0, unsigned _N>
		PointData<_NUMT, _LOC0, _N>::PointData(const unsigned &iin, const unsigned &in, const unsigned &ip)
		{
			if(_N == 0)
			{
				throw(std::runtime_error("N is not allowed to be zero"));
			}
			for(unsigned n = 0; n != _N; ++n)
			{
				pt0[n] = nullptr;
			}
			alloc(iin, in, ip);
		}

		template <typename _NUMT, LO::LOCATION _LOC0, unsigned _N>
		PointData<_NUMT, _LOC0, _N>::~PointData()
		{
			if(data != nullptr)
			{
				clear();
			}
		}

		template<typename _NUMT, LO::LOCATION _LOC0, unsigned _N>
		void PointData<_NUMT, _LOC0, _N>::alloc(const unsigned &iin, const unsigned &in, const unsigned &ip)
		{
			if(data != nullptr)
			{
				throw(std::runtime_error("allocating memory do not allowed when data is not empty"));
			}
			else
			{
				SizeT s(iin, in, ip);
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

		template<typename _NUMT, LO::LOCATION _LOC0, unsigned _N>
		void PointData<_NUMT, _LOC0, _N>::realloc(const unsigned &iin, const unsigned &in, const unsigned &ip)
		{
			SizeT(iin, in, ip).check();
			clear();
			alloc(iin, in, ip);
		}

		template<typename _NUMT, LO::LOCATION _LOC0, unsigned _N>
		void PointData<_NUMT, _LOC0, _N>::clear()
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

		template<typename _NUMT, LO::LOCATION _LOC0, unsigned _N>
		inline std::string PointData<_NUMT, _LOC0, _N>::disp() const
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

		template<typename _NUMT, LO::LOCATION _LOC0, unsigned _N>
		inline std::string PointData<_NUMT, _LOC0, _N>::disp_data() const
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
		inline void For_PD_1D(const RangeT &r, const Function &fun)
		{
			const int begin = r.begin(0);
			const int end = r.end(0);
			for (int i = begin; i != end; ++i)
			{
				fun(i);
			}
		}

		template <typename Function>
		inline void For_PD_1D_N(const RangeT &r, const unsigned &N_b, const unsigned &N_e, const Function &fun)
		{
			const int begin = r.begin(0);
			const int end = r.end(0);
			for (unsigned n = N_b; n != N_e; ++n)
			{
				for (int i = begin; i != end; ++i)
				{
					fun(i, n);
				}
			}
		}

		template <typename T, typename Reducer, typename Function>
		inline void Reduce_PD_1D(const RangeT &r, T &ans, const Reducer &reduce, const Function &fun)
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
		inline void Reduce_PD_1D_N(const RangeT &r, const unsigned &N_b, const unsigned &N_e, T &ans, const Reducer &reduce, const Function &fun)
		{
			T temp = ans;
			int begin = r.begin(0);
			int end = r.end(0);
			for (unsigned n = N_b; n != N_e; ++n)
			{
				for (int i = begin; i != end; ++i)
				{
					reduce(fun(i, n), temp);
				}
			}
			reduce(temp, ans);
		}
	}
}

#endif
