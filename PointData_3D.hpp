#ifndef POINTDATA_3D_HPP
#define POINTDATA_3D_HPP

#include <stdexcept>
#include <typeinfo>
#include <string>
#include <sstream>

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
					errlog << "out of DIM range: " << d << " range:[" << 0 << "," << static_cast<int>
					       (D) - 1 << "]";
					throw(std::runtime_error(errlog.str()));
				}
#endif
			}
		};

		class RangeT
		{
		  protected:
			int _begin[DIM::D] = { 0, 0, 0 };
			int _end[DIM::D] = { 0, 0, 0 };
			int _size[DIM::D] = { 0, 0, 0 };

		  public:
			RangeT() = default;
			RangeT(const int b0, const unsigned s0,
			       const int b1, const unsigned s1,
			       const int b2, const unsigned s2):
				_begin{ b0, b1, b2 },
				_size{ static_cast<int>(s0), static_cast<int>(s1), static_cast<int>(s2) }
			{
				_end[0] = _begin[0] + _size[0];
				_end[1] = _begin[1] + _size[1];
				_end[2] = _begin[2] + _size[2];
			}

			inline int begin(const unsigned d) const { DIM::check_d(d); return _begin[d]; }
			inline int end(const unsigned d) const { DIM::check_d(d); return _end[d]; }
			inline int size(const unsigned d) const { DIM::check_d(d); return _size[d]; }

			inline const int* begin_pt() const { return _begin; }
			inline const int* end_pt() const { return _end; }
			inline const int* size_pt() const { return _size; }

			inline std::string disp() const
			{
				std::ostringstream displog;
				displog << "range[0] = [" << _begin[0] << " , " << _end[0] << "]  "
				        << "size[0] = " << _size[0] << "    ";
				displog << "range[1] = [" << _begin[1] << " , " << _end[1] << "]  "
				        << "size[1] = " << _size[1] << "    ";
				displog << "range[2] = [" << _begin[2] << " , " << _end[2] << "]  "
				        << "size[2] = " << _size[2];
				return displog.str();
			}
		};

		class SizeT
		{
		  protected:
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
			inline int size(const unsigned d, RA::_ALL r) const { DIM::check_d(d); return _in[d] + _n[d] + _p[d]; }
			inline int begin(const unsigned d, RA::_IN r) const { DIM::check_d(d); return 0; }
			inline int end(const unsigned d, RA::_IN r) const { DIM::check_d(d); return _in[d]; }
			inline int size(const unsigned d, RA::_IN r) const { DIM::check_d(d); return _in[d]; }
			inline int begin(const unsigned d, RA::_N r) const { DIM::check_d(d); return -_n[d]; }
			inline int end(const unsigned d, RA::_N r) const { DIM::check_d(d); return 0; }
			inline int size(const unsigned d, RA::_N r) const { DIM::check_d(d); return _n[d]; }
			inline int begin(const unsigned d, RA::_P r) const { DIM::check_d(d); return _in[d]; }
			inline int end(const unsigned d, RA::_P r) const { DIM::check_d(d); return _in[d] + _p[d]; }
			inline int size(const unsigned d, RA::_P r) const { DIM::check_d(d); return _p[d]; }

			template<typename T0>
			inline int last(const unsigned d, T0 r) const { return end(d, r) - 1; }			inline int mirror(const unsigned d, FL::_N fl, int ii) const { return 2 * begin(d, RA::IN) - ii; }
			inline int mirror(const unsigned d, FL::_P fl, int ii) const { return 2 * last(d, RA::IN) - ii; }

			template<typename T0, typename T1, typename T2>
			inline RangeT range(T0 r0, T1 r1, T2 r2) const { return RangeT(begin(0, r0), size(0, r0), begin(1, r1), size(1, r1), begin(2, r2), size(2, r2)); }

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

			inline void check_range(const int i, const int j, const int k) const
			{
#ifdef _CHECK_POINTDATA_RANGE
				if (i < -_n[0] || i >= _in[0] + _p[0])
				{
					std::ostringstream errlog;
					errlog << "out of i range: " << i << " range:[" << -_n[0] << "," << _in[0] + _p[0] - 1 << "]";
					throw(std::runtime_error(errlog.str()));
				}
				if (j < -_n[1] || j >= _in[1] + _p[1])
				{
					std::ostringstream errlog;
					errlog << "out of j range: " << j << " range:[" << -_n[1] << "," << _in[1] + _p[1] - 1 << "]";
					throw(std::runtime_error(errlog.str()));
				}
				if (k < -_n[2] || k >= _in[2] + _p[2])
				{
					std::ostringstream errlog;
					errlog << "out of k range: " << k << " range:[" << -_n[2] << "," << _in[2] + _p[2] - 1 << "]";
					throw(std::runtime_error(errlog.str()));
				}
#endif
			}
		};

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

			template<typename F0, typename F1, typename F2>
			inline _NUMT &operator()(const unsigned n, const int i, const int j, const int k, const F0 flag0, const F1 flag1, const F2 flag2)
			{
				check_n(n);
				check_flag(flag0, flag1, flag2);
				_size.check_range(i, j, k);
				return pt0[n][i + F0::offset][j + F1::offset][k + F2::offset];
			}
			template<typename F0, typename F1, typename F2>
			inline const _NUMT &operator()(const unsigned n, const int i, const int j, const int k, const F0 flag0, const F1 flag1, const F2 flag2) const { return *this(i, j, k, n, flag0, flag1, flag2); }

			inline const _NUMT* data_pt(const unsigned n) const { check_n(n); return pt0[n][-_size.n(0)][-_size.n(1)] - _size.n(2); }
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

			template<typename F0, typename F1, typename F2>
			inline void check_flag(const F0 flag0, const F1 flag1, const F2 flag2)
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
					displog << "i = " << i << std::endl;
					for (int j = begin1; j != end1; ++j)
					{
						for (int k = begin2; k != end2; ++k)
						{
							displog << pt0[n][i][j][k] << ", ";
						}
						displog << std::endl;
					}
				}
				displog << std::endl;
			}
			return displog.str();
		}

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

	}
}

#endif
