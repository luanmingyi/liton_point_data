#ifndef POINTDATA_1D_HPP
#define POINTDATA_1D_HPP

#include <stdexcept>
#include <typeinfo>
#include <string>
#include <sstream>

#ifdef PD_OT
	#include "ordered_tec.h"
	#include <algorithm>
#endif

namespace liton_pd
{
	namespace D1
	{
		class DIM
		{
		  public:
			static const int D = 1;

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
			int _begin = 0;
			int _size = 0;

		  public:
			RangeT() = default;
			RangeT(const int b, const unsigned s): _begin(b), _size(static_cast<int>(s)) {}

			inline const int* begin_pt() const { return &_begin; }
			inline const int* size_pt() const { return &_size; }

			inline int begin(const unsigned d) const { DIM::check_d(d); return _begin; }
			inline int size(const unsigned d) const { DIM::check_d(d); return _size; }
			inline int end(const unsigned d) const { return begin(d) + size(d); }
			inline int last(const unsigned d) const { return end(d) - 1; }
			inline int bound(const unsigned d, FL::_N fl) const { return begin(d); }
			inline int bound(const unsigned d, FL::_P fl) const { return last(d); }

			inline RangeT &cut_head(int i) { _begin += i; _size -= i; return *this; }
			inline RangeT &cut_tail(int i) { _size -= i; return *this; }
			inline RangeT &tran(int i) { _begin += i; return *this; }

			inline bool is_overlap(const RangeT r) const
			{
				if(begin(0) >= r.end(0) || end(0) <= r.begin(0))
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
				int __begin = begin(0) > r.begin(0) ? begin(0) : r.begin(0);
				int __end = end(0) < r.end(0) ? end(0) : r.end(0);
				return RangeT(__begin, __end - __begin);
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
				        << "size[0] = " << size(0);
				return displog.str();
			}
		};
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		template <typename Function>
		inline void PD_For_1D(const RangeT r, const Function fun)
		{
			const int begin = r.begin(0);
			const int end = r.end(0);
			for (int i = begin; i != end; ++i)
			{
				fun(i);
			}
		}

		template <typename Function>
		inline void PD_For_N_1D(const unsigned N_b, const unsigned N_e, const RangeT r, const Function fun)
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
		inline void PD_Reduce_1D(const RangeT r, T &ans, const Reducer reduce, const Function fun)
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
		inline void PD_Reduce_N_1D(const unsigned N_b, const unsigned N_e, const RangeT r, T &ans, const Reducer reduce, const Function fun)
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
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		class SizeT
		{
		  public:
			int _in = 0;
			int _n = 0;
			int _p = 0;

		  public:
			SizeT() = default;
			SizeT(const unsigned in, const unsigned iin, const unsigned ip)
				: _in(static_cast<int>(iin)), _n(static_cast<int>(in)), _p(static_cast<int>(ip)) {}

			inline int in(const unsigned d) const { DIM::check_d(d); return _in; }
			inline int n(const unsigned d) const { DIM::check_d(d); return _n; }
			inline int p(const unsigned d) const { DIM::check_d(d); return _p; }

			inline int begin(const unsigned d, RA::_ALL r) const { DIM::check_d(d); return -_n; }
			inline int end(const unsigned d, RA::_ALL r) const { DIM::check_d(d); return _in + _p; }
			inline int begin(const unsigned d, RA::_IN r) const { DIM::check_d(d); return 0; }
			inline int end(const unsigned d, RA::_IN r) const { DIM::check_d(d); return _in; }
			inline int begin(const unsigned d, RA::_N r) const { DIM::check_d(d); return -_n; }
			inline int end(const unsigned d, RA::_N r) const { DIM::check_d(d); return 0; }
			inline int begin(const unsigned d, RA::_P r) const { DIM::check_d(d); return _in; }
			inline int end(const unsigned d, RA::_P r) const { DIM::check_d(d); return _in + _p; }
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

			template<typename T0>
			inline RangeT range(T0 r0) const { return RangeT(begin(0, r0), size(0, r0)); }

			std::string disp() const
			{
				std::ostringstream displog;
				displog << "size[0] = [" << _n << " " << _in << " " << _p << "]";
				return displog.str();
			}

			inline void check() const
			{
				if(_in == 0 && _n + _p != 0)
				{
					throw(std::runtime_error("dim[0]: size_in can not be zero when size_n or size_p is non-zero"));
				}
			}
		};
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
			PointData(const unsigned in, const unsigned iin, const unsigned ip);
			PointData(const PointData<_NUMT, _N, _LOC0> &) = delete;
			const PointData<_NUMT, _N, _LOC0> &operator=(const PointData<_NUMT, _N, _LOC0> &) = delete;
			~PointData();

			void alloc(const unsigned in, const unsigned iin, const unsigned ip);
			void realloc(const unsigned in, const unsigned iin, const unsigned ip);
			void clear();

			std::string disp() const;
			std::string disp_data() const;

			inline SizeT size() const { return _size; }

			template<typename F0 = FL::_C>
			inline _NUMT & operator()(const unsigned n, const int i, const F0 flag0 = FL::C)
			{
				check_data();
				check_n(n);
				check_flag(flag0);
				_size.range(RA::ALL).check_range(0, LOC0, i, F0::offset);
				return pt0[n][i + F0::offset];
			}
			template<typename F0 = FL::_C>
			inline const _NUMT & operator()(const unsigned n, const int i, const F0 flag0 = FL::C) const
			{
				check_data();
				check_n(n);
				check_flag(flag0);
				_size.range(RA::ALL).check_range(0, LOC0, i, F0::offset);
				return pt0[n][i + F0::offset];
			}

			inline _NUMT* data_pt(const unsigned n) { check_n(n); return pt0[n] - _size.n(0); }
			inline const _NUMT* data_pt(const unsigned n) const { check_n(n); return pt0[n] - _size.n(0); }
			inline _NUMT* _pt0(const unsigned n) { check_n(n); return pt0[n]; }
			inline const _NUMT* _pt0(const unsigned n) const { check_n(n); return const_cast<const _NUMT* >(pt0[n]); }

			void copy_from(const PointData<_NUMT, _N, _LOC0> &pd);
			template<typename __NUMT, unsigned __N>
			void copy_from(const PointData<__NUMT, __N, _LOC0> &pd,
			               const RangeT R, const RangeT r_,
			               unsigned N, unsigned n, int I, int i);
			template<typename __NUMT, unsigned __N, typename _R0, typename _r0, typename _FL0>
			inline void copy_from(const PointData<__NUMT, __N, _LOC0> &pd,
			                      _R0 R0, _r0 r0,
			                      unsigned N, unsigned n, _FL0 fl0)
			{
				copy_from(pd, _size.range(R0), pd.size().range(r0),
				          N, n, _size.bound(0, R0, fl0), pd.size().bound(0, r0, fl0));
			}

#ifdef PD_OT
			template<typename _R0, typename _FL0>
			void read_plt(const std::string &root, const liton_ot::TEC_FILE_LOG &teclog,
			              unsigned zone, const std::string &name,
			              _R0 R0, unsigned nn, _FL0 fl0);
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

			template<typename F0>
			inline void check_flag(const F0 flag0) const
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
		PointData<_NUMT, _N, _LOC0>::PointData(const unsigned in, const unsigned iin, const unsigned ip)
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
		void PointData<_NUMT, _N, _LOC0>::alloc(const unsigned in, const unsigned iin, const unsigned ip)
		{
			if(data != nullptr)
			{
				throw(std::runtime_error("allocating memory do not allowed when data is not empty"));
			}
			else
			{
				SizeT s(in, iin, ip);
				s.check();
				unsigned size_point = (s.size(0, RA::ALL) + _LOC0);
				if(size_point != 0)
				{
					try
					{
						data = new _NUMT[size_point * _N];
						_size = s;

						for (unsigned n = 0; n != _N; ++n)
						{
							pt0[n] = data + n * size_point + s.n(0);
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

		template<typename _NUMT, unsigned _N, LO::LOCATION _LOC0>
		void PointData<_NUMT, _N, _LOC0>::realloc(const unsigned in, const unsigned iin, const unsigned ip)
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
			displog << "dimension:" << DIM::D
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
			const int dln = 10;
			const int begin0 = _size.begin(0, RA::ALL);
			const int end0 = _size.end(0, RA::ALL) + _LOC0;
			for(int ii = 0; ii != dln * N; ++ii) { displog << "="; } displog << std::endl;
			for (int i = begin0; i != end0; ++i)
			{
				for (unsigned n = 0; n != N; ++n)
				{
					displog << pt0[n][i] << ", ";
				}
				displog << std::endl;
				if(i == _size.last(0, RA::N)
				        || i == _size.last(0, RA::IN))
				{for(int ii = 0; ii != dln * N; ++ii) { displog << "-"; } displog << std::endl;}
			}
			for(int ii = 0; ii != dln * N; ++ii) { displog << "="; } displog << std::endl;
			return displog.str();
		}

		template<typename _NUMT, unsigned _N, LO::LOCATION _LOC0>
		void PointData<_NUMT, _N, _LOC0>::copy_from(const PointData<_NUMT, _N, _LOC0> &pd)
		{
#ifdef _CHECK_POINTDATA_RANGE
			if(&pd == this)
			{
				throw(std::runtime_error("copy from self"));
			}
#endif
			RangeT overlap = _size.range(RA::ALL).overlap(pd.size().range(RA::ALL));
			overlap.cut_tail(-static_cast<int>(_LOC0));
			PD_For_N_1D(0, N, overlap, [&]PD_F_n_i(n, i)
			{
				pt0[n][i] = pd.pt0[n][i];
			});
		}

		template<typename _NUMT, unsigned _N, LO::LOCATION _LOC0>
		template<typename __NUMT, unsigned __N>
		void PointData<_NUMT, _N, _LOC0>::copy_from(const PointData<__NUMT, __N, _LOC0> &pd,
		        const RangeT R, const RangeT r_,
		        unsigned N, unsigned n, int I, int i)
		{
			check_n(N);
			pd.check_n(n);
			RangeT r = r_;
			RangeT overlap = R.overlap(r.tran(I - i));
#ifdef _CHECK_POINTDATA_RANGE
			if(reinterpret_cast<const void*>(&pd) == reinterpret_cast<const void*>(this))
			{
				RangeT RR = overlap;
				RR.tran(i - I);
				RR.cut_tail(-static_cast<int>(_LOC0));
				RR.cut_head(-static_cast<int>(_LOC0));
				if(RR.is_overlap(overlap))
				{
					throw(std::runtime_error("copy and write zone overlap"));
				}
			}
#endif
			overlap.cut_tail(-static_cast<int>(_LOC0));
			int ddi = i - I;
			PD_For_1D(overlap, [&]PD_F_i(II) {pt0[N][II] = static_cast<_NUMT>(pd._pt0(n)[II + ddi]);});
		}

#ifdef PD_OT
		template<typename _NUMT, unsigned _N, LO::LOCATION _LOC0>
		template<typename _R0, typename _FL0>
		void PointData<_NUMT, _N, _LOC0>::read_plt(
		    const std::string &root, const liton_ot::TEC_FILE_LOG &teclog,
		    unsigned zone, const std::string &name,
		    _R0 R0, unsigned nn, _FL0 fl0)
		{
			check_n(nn);
			if(zone >= teclog.Zones.size())
			{
				throw(std::runtime_error("zone number is too big"));
			}

			unsigned N0_file = teclog.Zones[zone].Real_Max_C(DIM::D, 0);
			if(N0_file == 0)
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
				PointData<float, 1, LO::center> buf(0, N0_file, 0);
				data_in.read(reinterpret_cast<char*>(buf.data_pt(0)), N0_file * teclog.Zones[zone].Data[vv].size);
				if (!data_in)
				{
					data_in.close();
					throw(std::runtime_error("error when reading file " + plt_file_name));
				}
				data_in.close();
				copy_from(buf, R0, RA::IN, nn, 0, fl0);
			}
			else if (teclog.Zones[zone].Data[vv].type == 2)
			{
				PointData<double, 1, LO::center> buf(0, N0_file, 0);
				data_in.read(reinterpret_cast<char*>(buf.data_pt(0)), N0_file * teclog.Zones[zone].Data[vv].size);
				if (!data_in)
				{
					data_in.close();
					throw(std::runtime_error("error when reading file " + plt_file_name));
				}
				data_in.close();
				copy_from(buf, R0, RA::IN, nn, 0, fl0);
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
