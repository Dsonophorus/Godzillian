#ifndef GDZ_HPP
#define GDZ_HPP

/*
 * Gdz.hpp  
 * version 2.0
 * 
 * Description:
 * ------------
 * This header file defines the Gdz - 'Godzillian' template class, which offers a highly flexible and efficient 
 * arbitrary-precision integer type. The `Gdz` class seamlessly integrates into C++ code, with syntax 
 * usage nearly indistinguishable from built-in types. It supports a wide range of operations, including 
 * arithmetic, bitwise, and shift operations, all while maintaining performance and ease of use.
 * 
 * Key Features:
 * -------------
 * - **Arbitrary Precision**: The `Gdz` class template allows for the creation of integers with arbitrary bit-length, 
 *   providing greater control and precision in mathematical computations.
 * - **Seamless Integration**: `Gdz` objects can be used just like standard integer types, including being used in 
 *   arithmetic expressions, comparisons, and type conversions. The syntax is intuitive and similar to built-in types.
 * - **Type Traits Specialization**: The `Gdz` class is specialized for use with standard type traits, allowing it to 
 *   be treated as an integral type by the standard library.
 * - **Literal Suffixes**: Supports custom literal suffixes (`_g2`, `_g8`, `_g10`, `_g16`) for binary, octal, decimal, 
 *   and hexadecimal initialization.
 * - **Comprehensive Operations**: Provides a full range of operations including addition, subtraction, multiplication, 
 *   division, modulus, bitwise operations, shifts, and rotations, all implemented with attention to performance.
 * - **Utility Functions**: Includes utility functions for operations such as counting bits, finding the most significant 
 *   bit (MSB) and least significant bit (LSB), reverse ordering, bit rotation.
 * - **Conversion Operators**: Easily convert `Gdz` objects to various integral and floating-point types, ensuring 
 *   compatibility with existing codebases.
 * - **Formatted String Conversion**: Use the `GdzTxt` class to format `Gdz` objects as strings with custom radix, 
 *   block size, and separators, enhancing the readability of large numbers.

 * Usage Examples:
 * ---------------
 * #define DO_CREATE
 * #include "Gdz.hpp"
 * 
 * int main() {
 *     Gdz<8> a = 123456789.987;			// Large floating-point initialization
 *     Gdz<8> b = "1ABC:1122:FFEE"_g16;     // Hexadecimal initialization using UDT
 *     Gdz<8> c = a + b + 10;               // Arithmetic operations
 * 
 *     double d = (double) c;				// Floating-point conversion
 * 
 *	   // Setting GdzT for formatting - 16 radix, 4 digit blocks, 0 : MSB,MSD left, no zero pad, : seperator, {} wrap
 *     GdzT.Set(16, 4, 0, ":", "{\a}");  
 *     std::cout << c << std::endl;			// output: {1ABC:187E:CD0D}
 
 * }
 * 
 * `DO_CREATE` should be defined once in the primary translation unit before including "Gdz.hpp" to instantiate
 * necessary static and global variables. Afterward, you can include "Gdz.hpp" anywhere else as needed.
 * 
 * License:
 * --------
 * This code is licensed under the Creative Commons Attribution-NonCommercial 4.0 International License (CC BY-NC 4.0).
 * You are free to use, share, and adapt the code for personal and educational purposes, provided that appropriate credit 
 * is given, and it is not used for commercial purposes without explicit permission.
 * 
 * For commercial use, please contact Mark Riphenburg at vonriphenburger@yahoo.com.
 * 
 * This code is copyright Mark Riphenburg, all rights reserved. [Date: August 2024]
 */

#ifdef DO_CREATE
#define mCreate( ... ) __VA_ARGS__
#define mNotCreate( ... )
#define mCreateElse( A, B ) A
#else
#define mCreate( ... )
#define mNotCreate( ... ) __VA_ARGS__
#define mCreateElse( A, B ) B
#endif


#include <type_traits>
#define _CX constexpr
using std::min;
using std::max;
typedef char    CHR;     /* ASCII 8 bit character */
typedef char*   STR;     /* null terminated 8 bit char pointer */
typedef const char* STR_;    /* null terminated const 8 bit char array */

// Utility type trait to strip all qualifiers including pointers
template <typename T> struct remove_all_qualifiers { using type = typename std::remove_cv<typename std::remove_pointer<typename std::remove_reference<T>::type>::type>::type; };
template <typename T> using remove_all_qualifiers_t = typename remove_all_qualifiers<T>::type;

// Primary template definition for _CXType
template<typename TT, typename Enable = void>
struct _CXType {
	using T = void;
	static _CX const char* Name() { return "UNDEFINED"; }
	static _CX TT Max() { return TT(); }
	static _CX TT Min() { return TT(); }
};

// Type by blocksize. use as _CXNType<-2>::T to emit S16, 4=U32
template<int Sz, typename Enable = void> struct _CXNType;

// Macro to define specializations for each type
#define _mRegNType(mN, mT, mTT, LMax, LMin) using mT = mTT; \
	template<> struct _CXType<mT> { using T = mT; \
		static _CX const char* Name() { return #mT; } \
		static _CX const mTT Max() { return LMax; } \
		static _CX const mTT Min() { return LMin; } }; \
	template<> struct _CXNType<mN> : _CXType<mT> { };

_mRegNType(1, U8, unsigned char, 255, 0)
_mRegNType(-1, S8, signed char, 127, -128)
_mRegNType(2, U16, unsigned short, 65535, 0)
_mRegNType(-2, S16, signed short, 32767, -32768)
_mRegNType(4, U32, unsigned long, 4294967295U, 0)
_mRegNType(-4, S32, signed long, 2147483647, -2147483648)
_mRegNType(8, U64, unsigned long long, 18446744073709551615ULL, 0)
_mRegNType(-8, S64, signed long long, 9223372036854775807LL, -(1LL << 63))
_mRegNType(-3, F32, float, 3.402823466e+38F, -3.402823466e+38F)
_mRegNType(-7, F64, double, 1.7976931348623158e+308, -1.7976931348623158e+308)


template<int Blk>
using _NType = typename _CXNType<Blk>::T;

// Utility functions for block operations
template<int Blk = 1> // count the length of a null-terminated string
_CX S32 _CXSCntZ(STR_ hL, S64 N = -1) { S32 cnt = 0;
	if (hL) for (N = N < 0 ? (U32)~0 : N; N > 0; --N, ++cnt, ++hL) if (!*hL) break;
	return cnt; // Return the length of the string
}

template<int Blk = 1> // byte copy R to L until zero value in R (0 is copied)
_CX S32 _CXSCpyZ(STR L, STR_ R, S64 N=-1) { S32 cnt = 0; 
	if( L && R ) for( N = N<0 ? (U32)~0 : N; N > 0; --N, ++cnt, ++R) { *L++ = *R; if( !*R ) break; }
	return cnt;
}

template<int Blk = 4,typename LType, typename RType> //copy R to L
_CX void _CXMCpy(LType& L, const RType& R, S64 N=-1) { //copy R to L, N is bytes not Blks
	auto* hL = reinterpret_cast<U32*>(&L); auto * hR = reinterpret_cast<const U32*>(&R);
	N = N<0 ? min(sizeof(L),sizeof(R))/4: N/4; for (S64 i = 0; i < N; ++i) hL[i] = hR[i];
}

template <typename T> _CX void _CXSwap_(T& L, T& R) { T temp = L; L = R; R = temp; }
auto _CXSwap = [](auto& a, auto& b) _CX { std::tie(a, b) = std::make_tuple(b, a); };


template<int TN> class Gdz;
template< typename T> class GdzTxt;
extern GdzTxt<Gdz<16>> GdzT;
template<typename T> // Flags: NumOrd:1 BlkOrd:2 Zeros:4 SignPos:8 Plus:16
STR _CXToStr(const T& Num, U32 Radix = 10, S32 DBlk = 0, U32 Flags = 0, const char* iSep = ":", const char* iWrap = "", STR hOut = nullptr);
extern S8 gTblBit[256][4];
template< typename T>
class GdzTxt { public:
	enum { NumOrd = 1, BlkOrd = 2, Zeros = 4, SignPos = 8, Plus = 16, Default = 256 };
	U32 _Radix, _Blk, _Flags;
	CHR _Sep[32], _Wrap[32]; // \a delimited
	STR hOut; T * hNum;
	GdzTxt() : _Radix(16), _Blk(8), _Flags(0), _Sep(":"), _Wrap("{\a}"), hOut(nullptr), hNum(0) {};
	GdzTxt( T& iNum ) : _Radix(16), _Blk(8), _Flags(0), _Sep(":"), _Wrap("{\a}"), hOut(nullptr), hNum(&iNum) {}
	template<typename U> GdzTxt( GdzTxt<U> from ) { _CXMcpy(*this, from); }
	operator STR() const { if( !hNum ) return (STR)""; return _CXToStr( *hNum, _Radix, _Blk, _Flags, _Sep, _Wrap, hOut); }
	GdzTxt& Set( U32 iRadix = 16, U32 iBlk = 8, U32 iFlags = 0, STR_ iSep = ":", STR_ iWrap = "{\a}" )
	{ _Radix=iRadix; _Blk=iBlk; _Flags=iFlags; _CXSCpyZ(_Sep,iSep,32); _CXSCpyZ(_Wrap,iWrap,32); 
	   if(_Flags &Default){ _Flags &= ~Default; _CXMCpy(GdzT, *this); }; return *this; }
	GdzTxt& Radix(U32 I=10) { _Radix=I; return *this; } 
	GdzTxt& Blk(U32 I=~0) { _Blk=I; return *this; }
	GdzTxt& Flags(U32 I=0) { _Flags=I; if(_Flags &Default){ _Flags &= ~Default; _CXMCpy(GdzT, *this); }; return *this; }
	GdzTxt& Sep(STR hI=0) { _CXSCpyZ(_Sep,hI,32); if(!hI) *_Sep = 0; return *this; }
	GdzTxt& Wrap(STR hI=0) { _CXSCpyZ(_Wrap,hI,32); if(!hI) *_Wrap = 0; return *this; }
	GdzTxt& Num(T& iNum) { hNum = &iNum; return *this; }
	// Overloaded << operator for chaining
	template<typename T>
	auto& operator<<(const T& value) { return std::cout << (STR) value; }
};

template<int TN>
class Gdz {
public:
	static constexpr size_t Sz = TN > 0 ? TN : -TN; // Determine size based on TN (absolute value of TN)
	static constexpr bool is_signed_v = TN < 0;     // Determine if TN is signed (if TN is negative, then it must be signed)
	U32 data[Sz];

#if 1 // constructor fold
	// Default constructor
	_CX Gdz() : data{} {}
	// forwarding copy constructor
	_CX Gdz(const Gdz& R) { for(S32 i = 0; i < Sz; ++i) data[i]=R.data[i]; }

	// Copy constructor for different sizes
	template<U32 M>
	_CX Gdz(const Gdz<M>& other) { U32 i = 0;
		for (; i < min((U32)Sz, M); ++i) data[i] = other.data[i];
		U32 sgn = other.sign() ? 0xFFFFFFFF : 0;
		for (; i < Sz; ++i)	data[i] = sgn;
	}

	_CX Gdz(std::initializer_list<U32> init) : data{} { U32 i = 0, sgn;
		for (auto it = init.begin(); it != init.end() && i < Sz; ++it, ++i)	data[i] = *it;
		for (sgn = 0xFFFFFFFF * ((0x80000000 & data[i?i-1:0]) != 0); i < Sz; data[i++] = sgn); // -|0 extend
	}

	//construct from S64
	template <typename T, typename = std::enable_if_t<std::is_integral<T>::value && !std::is_floating_point<T>::value>>
	_CX Gdz(const T V) : data{} {
		U32 sgn = 0xFFFFFFFF * (V < 0), i = 0;  // Determine sign for extension
		data[i++] = (U32)(V % 0x100000000ull);
		if(V / 0x100000000ull != 0) data[i++] = (U32) (V/0x100000000ull);
		for(; i < Sz; data[i++] = sgn); //-|0 extend
	}

	// Construct from F64
	Gdz(const F64 iV) : data{} { //not _CX due to float ops
		const F64 Rx = (F64) 0x100000000ull;
		bool sgn = (iV < 0.0), carry = 1; U32 i = 0;
		F64 V = std::fabs(iV), intPart = 0.0; // Work with the absolute V for digit extraction
		for (; i < Sz; ++i)
			if (V < 1.0) data[i] = 0;
			else {
				intPart = std::floor(V);
				data[i] = static_cast<U32>(intPart - std::floor(intPart / Rx) * Rx);
				V = (intPart - data[i]) / Rx;
			}
		if (V != 0.0) std::cerr << "Warning: Value truncated during Gdz construction\n";    
		if (sgn) neg(); //must put the number in 2~ form because it is negative
		sign(sgn);		// force top bit to correct sign 
	}

	// Multiple homogeneous arguments constructor
	template <typename T, typename... Args,
			  typename = std::enable_if_t<(sizeof...(Args) > 0) && std::conjunction_v<std::is_same<T, Args>...>>>
	_CX Gdz(const T first, const Args... args) {
		// Cast all args to U32, Last arg determines sign
		U32 temp[] = { (U32)first, (U32)args... };
		U32 sgn = 0xFFFFFFFF * ((0x80000000 & temp[sizeof(temp) / 4 - 1]) != 0);
		U32 i = 0;
		for (; i < std::min(sizeof(temp) / 4, Sz); ++i) data[i] = temp[i];
		for(; i < Sz; data[i++] = sgn);  //-|0 extend
	}

	// Multiple heterogeneous arguments constructor
	template <typename... Args, typename = std::enable_if_t<(sizeof...(Args) > 1 && !std::conjunction_v<std::is_same<U32, Args>...>)>>
	_CX Gdz(const Args... args) {
		static_assert(sizeof...(args) <= Sz, "Too many arguments");

		size_t i = 0;
		auto process_arg = [&](auto arg) {
			using ArgType = std::decay_t<decltype(arg)>;
			if constexpr (std::is_same_v<ArgType, Gdz>) {
				// Copy data from another Gdz, handling different sizes
				size_t len = std::min(Sz - i, ArgType::Sz);
				for (size_t j = 0; j < len; ++j) {
					data[i++] = arg.data[j];
				}
			} else if constexpr (std::is_floating_point_v<ArgType>) {
				// Convert floating-point to U32 representation
				F64 intPart;
				F64 fracPart = std::modf(static_cast<F64>(arg), &intPart);
				data[i++] = static_cast<U32>(intPart);
				while (i < Sz && fracPart != 0.0) {
					fracPart = std::modf(fracPart * std::numeric_limits<U32>::max(), &intPart);
					data[i++] = static_cast<U32>(intPart);
				}
			} else {
				// Handle other types (e.g., integral types)
				data[i++] = static_cast<U32>(arg);
			}
		};
		(process_arg(args), ...); // Process each argument

		// Determine the sign based on the most significant bit of the last initialized element
		U32 sgn = 0xFFFFFFFF * ((0x80000000 & data[i - 1]) != 0);

		// Extend the remaining elements with the sign
		for (; i < Sz; ++i) {
			data[i] = sgn;
		}
	}

	// Sequence constructor
	template <std::size_t tN>
	_CX Gdz(const U32(&init)[tN]) {	static_assert(tN <= Sz, "Initializer list too large");
		U32 sgn = 0xFFFFFFFF * ((0x80000000 & init[tN - 1]) != 0), i = 0; // Determine sign for extension
				for (; i < tN && i < Sz; ++i) data[i] = init[i]; // Copy initializer list elements to data array
		for (; i < Sz; ++i) data[i] = sgn; // Sign extend or zero-fill the remaining elements
	}
#endif //constructor fold

// Sequence iterator interface
	using value_type = U32;
	using iterator = value_type*;
	using const_iterator = const value_type*;
	using size_type = size_t;
	using reference = value_type&;
	using const_reference = const value_type&;
	_CX iterator begin() { return data; }
	_CX const_iterator begin() const { return data; }
	_CX const_iterator cbegin() const { return data; }
	_CX iterator end() { return data + Sz; }
	_CX const_iterator end() const { return data + Sz; }
	_CX const_iterator cend() const { return data + Sz; }
	_CX size_type size() const { return Sz; }
	_CX size_type max_size() const { return Sz; }
	_CX reference operator[](size_type Idx) { return data[Idx]; }
	_CX const_reference operator[](size_type Idx) const { return data[Idx]; }
	_CX reference at(size_type Idx) { if (Idx >= Sz) throw std::out_of_range("Index out of range"); return data[Idx]; }
	_CX const_reference at(size_type Idx) const { if (Idx >= Sz) throw std::out_of_range("Index out of range"); return data[Idx]; }

// cast operators
	_CX operator U32*() { return data; }
	_CX operator const U32*() const { return data; }
	_CX explicit operator bool() const { return _CXMsd(*this) >= 0; }
	_CX explicit operator U64() const { if _CX(sizeof(Gdz)<8) return (S64) data[0]; else return *(U64*)this; }
	_CX explicit operator S64() const { if _CX(sizeof(Gdz)<8) return (S64) data[0]; else return *(S64*)this; }
	_CX explicit operator U32() const { return data[0]; }
	_CX explicit operator S32() const { return (S32)data[0]; }
	_CX explicit operator U16() const { return (U16)data[0]; }
	_CX explicit operator S16() const { return (S16)data[0]; }
	_CX explicit operator U8() const { return (U8)data[0]; }
	_CX explicit operator S8() const { return (S8)data[0]; }
	_CX operator STR() const { return _CXToStr( *this, GdzT._Radix, GdzT._Blk, GdzT._Flags, GdzT._Sep, GdzT._Wrap, GdzT.hOut); }
	template<int tN> _CX explicit operator Gdz<tN>() const { Gdz<tN> RVal; _CXMCpy(&RVal, this, min(sizeof(RVal), sizeof(Gdz))); return RVal; }
	operator F32() const { return (F32)(F64)*this; }
	operator F64() const { const F64 Rx = (F64) 0x100000000ull;	F64 result = 0.0, factor = 1.0; S32 i = 0;
		if (data[Sz - 1] & 0x80000000) {
			Gdz Tmp(*this); Tmp.neg();
			for (; i < Sz; ++i, factor *= Rx) result += static_cast<F64>(Tmp.data[i]) * factor;
			result = -result;
		} else for (; i < Sz; ++i, factor *= Rx) result += static_cast<F64>(data[i]) * factor;

		return result;
	}	
	_CX explicit operator unsigned int() const { return (unsigned int)data[0]; }	// these resolve gay ambiguity
	_CX explicit operator int() const { return (int)data[0]; }			

// comparison operations
	template<typename T>
	_CX S64 cmp(const T& RV) const {
		if constexpr (std::is_floating_point_v<T>) {
			return cmp(Gdz<Sz+1>(RV));
		} else if constexpr (sizeof(T) == 8) {
			return cmp(Gdz<2>(RV));
		} else {
			U32 LSgn = sign(), RSgn = RV < 0; S64 LSz = Sz; U32 R = (U32) RV;

			if (LSgn && RV >= 0) return -1; // L is negative, R is positive
			if (!LSgn && RV < 0) return 1;  // L is positive, R is negative

			LSgn *= ~0LL;

			while (--LSz > 0) if (data[LSz] != LSgn) return LSgn ? -++LSz : ++LSz; // check for value in excess
			return *data == R ? 0 : *data < R ? -1 : 1;
		}
	}

	template<int RN>
	_CX S64 cmp(const Gdz<RN>& R) const { //returns <=-,>=+,==0 the magnitude the position of the first difference +1
		U32 LSgn =  ~0LL * sign(), RSgn = ~0LL*R.sign(); S64 LSz = Sz, RSz = R.Sz;

		if (LSgn && !RSgn) return -max(LSz,RSz); // L is negative, R is positive
		if (!LSgn && RSgn) return max(LSz, RSz); // L is positive, R is negative

		while (--RSz > LSz) if (R.data[RSz] != RSgn) return RSgn ? ++RSz : -++RSz;
		while (--LSz > RSz) if (data[LSz] != LSgn) return LSgn ? -++LSz : ++LSz;
		//if this point is reached, LSz==RSz and are pointing to the MSD and LSgn==RSgn
		for (; LSz >= 0;--LSz) if (data[LSz] != R.data[LSz]) return data[LSz] < R.data[LSz] ? -++LSz : ++LSz;

		return 0;
	}

#define _Op(Cmp) template<typename T> _CX bool operator##Cmp (const T& Rhs) const { return cmp(Rhs) Cmp 0; }
	_Op(==) _Op(!=) _Op(<) _Op(>) _Op(<=) _Op(>=)
#undef _Op

// utility operations
	_CX Gdz& set(U32 V=0) { for(S32 i=0; i<Sz; data[i++]=V); return *this; }
	_CX bool sign() const { return (data[Sz-1] & 0x80000000) != 0; }
	_CX Gdz& sign(bool sbit) { if(sbit) data[Sz-1] |= 0x80000000; else data[Sz-1] &= ~0x80000000; return *this; }
	_CX Gdz& neg() { bool carry = true;
		for (S32 i = 0; i < Sz; ++i) {
			data[i] = ~data[i];
			if (carry) {
				carry = (data[i] == 0xFFFFFFFF);
				data[i]++;
			}
		}
		return *this;
	}
	_CX Gdz& abs() { if (sign()) neg(); return *this; };
	_CX Gdz operator-() const { bool carry = true; Gdz RVal;
		for (S32 i = 0; i < Sz; ++i) {
			RVal.data[i] = ~data[i];
			if (carry) {
				carry = (RVal.data[i] == 0xFFFFFFFF);
				RVal.data[i]++;
			}
		}
		return RVal;
	}
	_CX Gdz operator+() const { return (data[Sz - 1] & 0x80000000) != 0 ? -*this : *this; }
	_CX Gdz& Not(){ for(U32 i=0; i<Sz;++i) data[i] = ~data[i]; return *this; } 
	_CX Gdz& Rev(const bool& fBits = 0){ U32 i = 0, T=0;	// rev blk ordering. if fBits then bits are also reversed per blk
	if (!fBits) for (; i < Sz / 2; ++i){ T = data[i]; data[i] = data[Sz - 1 - i]; data[Sz - 1 - i] = T; }
		else for (; i < Sz / 2; ++i) { T = RevV32(data[i]); data[i] = RevV32(data[Sz - 1 - i]); data[Sz - 1 - i] = T; };
		return *this;
	}

	//return the Blk index of the lowest non-zero element of data
	_CX S64 Lsd() const {	for (S64 i = 0; i < Sz; ++i) if (data[i] != 0) return i; return -1; }
	//return the bit position of the lowest non-zero bit
	_CX S64 Lsb() const { S64 R = Lsd(); return R == -1 ? -1 : R * 32 + LsbV32(data[R]); }
	//return the Blk index of the uppper most non-zero element of data
	_CX S64 Msd() const {	for (S64 i = Sz - 1; i >= 0; --i) if (data[i] != 0) return i; return -1; }
	//return the bit position of the uppper most non-zero bit for a Blk array
	_CX S64 Msb() const { S64 R = Msd(); return R == -1 ? -1 : R * 32 + MsbV32(data[R]); }
	_CX U64 CntBits() const { U64 Cnt = 0; for (S64 i = 0; i < Sz; ++i) Cnt += CntBitsV32(data[i]); return Cnt; }

	_CX static U32 CntBitsV32(U32 V) { return gTblBit[(U8)(V & 0xFF)][0] + gTblBit[(U8)((V >> 8) & 0xFF)][0] + gTblBit[(U8)((V >> 16) & 0xFF)][0] + gTblBit[(U8)((V >> 24) & 0xFF)][0]; }	
	_CX static S32 LsbV32(U32 V) { return (V & 0xFF ? gTblBit[(U8)(V & 0xFF)][1] : ((V >> 8) & 0xFF ? 8 + gTblBit[(U8)((V >> 8) & 0xFF)][1] : 
			   ((V >> 16) & 0xFF ? 16 + gTblBit[(U8)((V >> 16) & 0xFF)][1] : ((V >> 24) & 0xFF ? 24 + gTblBit[(U8)((V >> 24) & 0xFF)][1] : -1))));
	}
	_CX static S32 MsbV32(U32 V) { return ((V >> 24) & 0xFF ? 24 + gTblBit[(U8)((V >> 24) & 0xFF)][2] : ((V >> 16) & 0xFF ? 16 + gTblBit[(U8)((V >> 16) & 0xFF)][2] : 
			   ((V >> 8) & 0xFF ? 8 + gTblBit[(U8)((V >> 8) & 0xFF)][2] : (V & 0xFF ? gTblBit[(U8)(V & 0xFF)][2] : -1))));
	}
	_CX static U32 RevV32(U32 V) { U32 R = (U32)(U8)gTblBit[(U8)(V & 0xFF)][3] << 24;
	R |= (U32)(U8)gTblBit[(U8)((V >> 8) & 0xFF)][3] << 16;
		R |= (U32)(U8)gTblBit[(U8)((V >> 16) & 0xFF)][3] << 8; 
		R |= (U32)(U8)gTblBit[(U8)((V >> 24) & 0xFF)][3]; return R;
	}


// bit operations
#define mBitOps(Name, Sym) template<typename T> _CX void Name(const T& R, Gdz& oL) const { \
		if _CX (requires { T::Sz; }) { S32 i = 0; \
			for (; i < min(Sz, R.Sz); ++i) oL.data[i] = data[i] Sym R.data[i]; \
			for (U32 sgnExt = R.sign() ? 0xFFFFFFFF : 0; i < Sz; ++i) oL.data[i] = data[i] Sym sgnExt; \
		} else { S32 i = 0;  oL.data[i++] = data[i] Sym (U32)R; \
			if _CX (sizeof(T) > 4) oL.data[i++] = data[i] Sym (U32)(R >> 32); \
			for (U32 sgnExt = R < 0 ? 0xFFFFFFFF : 0; i < Sz; ++i) oL.data[i] = data[i] Sym sgnExt; } } \
	template<typename T>  Gdz operator Sym(const T& Rhs) const { Gdz RVal; Name(Rhs, RVal); return RVal; } \
	template<typename T>  Gdz& operator Sym##=(const T& Rhs) { Name(Rhs, *this); return *this; }
	mBitOps(And, &) mBitOps(Or, |) mBitOps(Xor, ^)

	Gdz operator~() const { Gdz RVal; for(U32 i=0; i<Sz;++i) RVal.data[i] = ~data[i]; return RVal; }
#undef mBitOps

// shift operations
	_CX Gdz& LShift(const U32 V, Gdz& oL) const {
		if (!V) return oL = *this; S32 VBlk = V / 32, VBit = V % 32, i = 0;
		if (VBlk >= Sz) { while (i < Sz) oL.data[i++] = 0; return oL; }

		// 32-bit optimization
		if (!VBit) { for (i = Sz - 1; i >= VBlk; --i) oL.data[i] = data[i - VBlk];
			while (i >= 0) oL.data[i--] = 0; return oL; }

		// Perform the shift in place from end to beginning
		for (i = Sz - 1; i >  VBlk; --i) 
			oL.data[i] = (U32)((((U64)data[i - VBlk - 1 ]) | (((U64)data[i - VBlk]) << 32)) >> (32 - VBit));		

		oL.data[i] = data[i - VBlk] << VBit;
		while (--i >= 0) oL.data[i] = 0;
		return oL;
	}

	_CX Gdz& RShift(U32 V, Gdz& oL, bool fSigned = 1) const {   // sign preserving rshift - force Sgn = 0 for non-arithmatic 
		if (!V) return oL=*this; U32 Sgn = (sign() && fSigned) ? ~0U : 0; S32 VBlk = V / 32, VBit = V % 32, i = 0;
		if (VBlk >= Sz) { while (i < Sz) oL.data[i++] = Sgn; return oL; }

		// 32-bit optimization
		if (!VBit) { for (i = 0; i < (S64)Sz - VBlk; ++i) oL.data[i] = data[i + VBlk];
			while (i < Sz) oL.data[i++] = Sgn; return oL; }

		// Perform the shift in place from beginning to end
		for( i = 0; i < Sz - VBlk - 1; i++)
			oL.data[i] = (U32)((((U64)data[i + VBlk]) | (((U64)data[i + VBlk + 1]) << 32)) >> VBit );

		oL.data[i] = (U32)(((((U64) Sgn) <<32 ) | data[i+VBlk]) >> VBit);
		while (++i < Sz) oL.data[i] = Sgn;
		return oL;
	}

	_CX Gdz& Shift(S32 V, Gdz& oL, bool fSigned = 1 ) { return V<0 ? RShift(-V,oL, fSigned) : LShift(V,oL); }

	_CX Gdz& Rotate(S32 V) {if (V == 0) return *this;
		V = V % (32 * Sz); if (V < 0) V += (32 * Sz);  // Normalize the rotation and convert - right rotate to equivalent left rotate
		S64 VBlk = V / 32, VBit = V % 32, i = 0; U32 temp[Sz] = {0};
		if (VBlk >= Sz) VBlk %= Sz; // Block rotation

		// 32-bit optimization
		if (VBit == 0) {
			for (i = 0; i < Sz; ++i) temp[(i + VBlk) % Sz] = data[i];
			for (i = 0; i < Sz; ++i) data[i] = temp[i];
			return *this;
		}

		// Perform the shift in place
		for (i = 0; i < Sz; ++i) {
			S64 idx = (i + VBlk) % Sz;
			temp[idx] = (data[i] << VBit) | (data[(i + 1) % Sz] >> (32 - VBit));
		}
		for (i = 0; i < Sz; ++i) data[i] = temp[i];
		return *this;
	}

#define mShiftOps(Name, Sym) \
	_CX Gdz operator Sym(U32 V) const { Gdz RVal; return Name(V,RVal); } \
	_CX Gdz& operator Sym##=(U32 V) { return Name(V,*this); }
	mShiftOps(LShift, << ) mShiftOps(RShift, >> )
#undef mShiftOps

	// Add function
	_CX Gdz& operator++()   { Add(1, *this); return *this; }                    // Pre-increment
	_CX Gdz  operator++(int){ Gdz temp = *this; Add(1, *this); return temp; }   // Post-increment
	_CX Gdz& operator--()   { Sub(1, *this); return *this; }                    // Pre-decrement
	_CX Gdz  operator--(int){ Gdz temp = *this; Sub(1, *this); return temp; }   // Post-decrement

	template<typename T> 
	_CX Gdz& Add(const T& R, Gdz& oL) const { S32 i = 0; U64 Acc = 0;
		if _CX(requires { T::Sz; }) {
			for (; i < min(Sz, R.Sz); ++i) {
				Acc += (U64)data[i] + (U64)R.data[i];
				oL.data[i] = (U32)Acc;
				Acc >>= 32;
			}
			for (; i < Sz; ++i) {
				Acc += (U64)data[i];
				oL.data[i] = (U32)Acc;
				Acc >>= 32;
			}
		} else {
			Acc = (U64)R;
			for (i = 0; i < Sz; ++i) {
				Acc += (U64)data[i];
				oL.data[i] = (U32)Acc;
				Acc >>= 32;
			}
		}
		return oL;
	}

	template<typename T> 
	_CX Gdz& Sub(const T& R, Gdz& oL) const { S32 i = 0;S64 Acc = 0;
		if _CX(requires { T::Sz; }) {
			for (; i < min(Sz, R.Sz); ++i) {
				Acc += (S64)data[i] - (S64)R.data[i];
				oL.data[i] = (U32)Acc;
				Acc >>= 32;
			}
			for (; i < Sz; ++i) {
				Acc += (S64)data[i];
				oL.data[i] = (U32)Acc;
				Acc >>= 32;
			}
		} else {
			Acc = -(S64)R;
			for (i = 0; i < Sz; ++i) {
				Acc += (S64)data[i];
				oL.data[i] = (U32)Acc;
				Acc >>= 32;
			}
		}
		return oL;
	}

template<typename T>
_CX Gdz& Mul(const T& R, Gdz& oL) const {
	bool resultNegative = (this->sign() != 0) ^ (R < 0); // Calculate result sign
	Gdz TmpL = +*this; // Copy the current object
	Gdz TmpR = +R; // Copy the R value

	S32 iL = 0, iR = 0;
	U64 Acc = 0; // Use U64 to handle overflow properly

	// Zero out the result
	for (iL = 0; iL < Sz; ++iL) oL.data[iL] = 0;

	if _CX(requires { T::Sz; }) {
		// Perform multiplication for Gdz type
		for (iL = 0; iL < Sz; ++iL) {
			for (Acc = 0, iR = 0; iR < TmpR.Sz && ((S64)iL + iR) < Sz; ++iR) {
				Acc += (U64)TmpL.data[iL] * (U64)TmpR.data[iR] + (U64)oL.data[iL + iR];
				oL.data[iL + iR] = (U32)Acc;
				Acc >>= 32;
			}
			while (((S64)iL + iR) < Sz && Acc) {
				Acc += (U64)oL.data[iL + iR];
				oL.data[iL + iR++] = (U32)Acc;
				Acc >>= 32;
			}
		}
	} else {
		// Perform multiplication for integral type
		U32 rData = (U32)TmpR; // Use the absolute value of R
		for (iL = 0; iL < Sz; ++iL) {
			Acc += (U64)TmpL.data[iL] * rData + (U64)oL.data[iL];
			oL.data[iL] = (U32)Acc;
			Acc >>= 32;
		}
	}

	if (resultNegative) oL.neg(); // Apply the sign to the result if necessary
	return oL;
}

template<int TN>
_CX U32 DivV32(U32 Rv, Gdz<TN>& O) {
	Gdz<TN> LL(+*this); U64 Acc = 0; 
	 
	for (S32 i = LL.Sz - 1; i >= 0; --i) { 
		Acc = (Acc << 32) + LL.data[i]; 
		O.data[i] =(U32)(Acc / Rv); 
		Acc %= Rv; 
	}

	return (U32)Acc;
}

template<int TN>
_CX Gdz& DivL( const Gdz<TN>& _div, Gdz& q, Gdz& dvd ) const { // dvd hold modulus/remainder
	bool sgn = sign() + _div.sign(); Gdz<Sz+2> mdiv, div = _div; dvd = *this; div.abs(); dvd.abs(); q.set(0);
	S32 dvdM, divM=(S32)div.Msd(), shift, i = 0; U64 dvdA, divA = div[divM]; 

	if( div == 0 ) { if(dvd==0) q = 1; else{ dvd.set(0); q.set(~(U32)!sgn).sign(sgn); } return q; } //handle div 0 as +-max or 1 if 0/0

	while( dvd > div ){
		if(divA > (dvdA = dvd[dvdM = (S32)dvd.Msd()]) && dvdM-1 >= divM )
			{ dvdA <<=32; dvdA += dvd[--dvdM]; }
		mdiv = div; //possible overflow
		mdiv <<= 32*(dvdM-divM);
		dvd -= mdiv * (q[dvdM-divM] = (U32)(dvdA/divA));

		if( dvd.sign() && q[dvdM-divM] > 0 )
		{
			Gdz<10> mDiv1 =mdiv;
			for( shift = 32, mdiv <<= 5; dvd.sign() && q[dvdM-divM] > shift; shift <<= 1, mdiv <<=1 )	// fix over subtraction
				{	dvd += mdiv; q[dvdM-divM]-=shift;  cout << "+";}
			for( shift >>= 2, mdiv >>= 2; !dvd.sign(); )												// fix over addition
				{	dvd -= mdiv; q[dvdM-divM]+=shift; if( shift >1 ) { shift >>= 1; mdiv >>=1; } cout << "-"; }
		
			for( mdiv >>=LsbV32(shift); dvd.sign() && q[dvdM-divM] > 0; q[dvdM-divM]--)					// fix final subtraction
				{dvd += mdiv; cout << "*"; }
		}
	}
	if( sgn ) { dvd.neg(); q.neg(); } // use this for standard -mod handling instead of remainder: { dvd.neg() += div; q.neg()++; }

	return q;
}

template<typename T> _CX Gdz& Div(const T& R, Gdz& oL) const { Gdz Rem; if _CX(requires { T::Sz; }) 
	return DivL( R, oL, Rem); else { Gdz<2> tR = R; return DivL( tR, oL, Rem); } }
template<typename T> _CX Gdz& Mod(const T& R, Gdz& oRem) const { Gdz oL; if _CX(requires { T::Sz; }) DivL( R, oL, oRem);
	else { Gdz<2> tR = R; DivL( tR,oL, oRem); } return oRem; }


#define mArithOps(Name, Sym) \
	template<typename T> Gdz operator Sym(const T& Rhs) const { Gdz RVal(*this); return Name(Rhs, RVal); } \
	template<typename T> Gdz& operator Sym##=(const T& Rhs) { return Name(Rhs, *this); }
	mArithOps(Add, +) mArithOps(Sub, -) mArithOps(Mul, *) mArithOps(Div, /) mArithOps(Mod, %)
#undef mArithOps

};

using S128 = Gdz<4>;
using S256 = Gdz<8>;
using S512 = Gdz<16>;
using S1K  = Gdz<32>;
using GdzL = Gdz<16>;

// General specialization for Gdz types, only accepts N > 4
template<int tN>
struct _CXType<Gdz<tN>, typename std::enable_if<(tN > 2)>::type> { // enable for typ above 64 bits to avoid collision with builtin
	using T = Gdz<tN>;
	static _CX const char* Name() { static char name[12]; snprintf(name, sizeof(name), "Gdz<%d>", tN); return name; }
	static _CX T Max() { return T(-1).sign(0); }
	static _CX T Min() { return T(0).sign(1); }
};

template<int Size>
struct _CXNType<Size, typename std::enable_if<(Size < -8 || Size > 8)>::type> // create integer by size
	: _CXType<Gdz<Size>> {
	using T = typename _CXType<Gdz<Size>>::T;
};

// Specialize type traits for Gdz
namespace std {
	template<int tN> struct std::is_integral<Gdz<tN>> :     true_type {};
	template<int tN> struct std::is_arithmetic<Gdz<tN>> :   true_type {};
	template<int tN> struct std::is_signed<Gdz<tN>> :       true_type {};
	template<int tN> struct std::is_unsigned<Gdz<tN>> :     false_type {};
	template<int tN> struct std::is_constructible<Gdz<tN>> :true_type {};
	template<int tN> struct std::is_default_constructible<Gdz<tN>> :true_type {};
	template<int tN> struct std::is_copy_constructible<Gdz<tN>> :   true_type {};
	template<int tN> struct std::is_move_constructible<Gdz<tN>> :   true_type {};
};

#ifdef DO_LIMITS
namespace std {
	template<int tN>
	class numeric_limits<Gdz<tN>> {
	public:
		static constexpr bool is_specialized = true;
		static constexpr Gdz<tN> min() noexcept { return CXType < Gdz<tN>::Min(); }
		static constexpr Gdz<tN> max() noexcept { return CXType < Gdz<tN>::Max(); }
		static constexpr Gdz<tN> lowest() noexcept { return CXType < Gdz<tN>::Min(); }
		static constexpr int digits = std::numeric_limits<int32_t>::digits * tN;
		static constexpr int digits10 = static_cast<int>(digits * 0.30102999566); // log10(2) ~ 0.30102999566
		static constexpr bool is_signed = true;
		static constexpr bool is_integer = true;
		static constexpr bool is_exact = true;
		static constexpr int radix = 2;
		static constexpr bool is_bounded = true;
		//most floating point+ specific number traits do not apply to Gdz
		static constexpr Gdz<tN> epsilon() noexcept { return Gdz<tN>(0); }
		static constexpr Gdz<tN> round_error() noexcept { return Gdz<tN>(0); }
		static constexpr int min_exponent = 0;
		static constexpr int min_exponent10 = 0;
		static constexpr int max_exponent = 0;
		static constexpr int max_exponent10 = 0;
		static constexpr bool has_infinity = true;
		static constexpr bool has_quiet_NaN = true;
		static constexpr bool has_signaling_NaN = true;
		static constexpr float_denorm_style has_denorm = std::denorm_absent;
		static constexpr bool has_denorm_loss = false;
		static constexpr Gdz<tN> infinity() noexcept { return _CXType<Gdz<tN>>::Max(); }
		static constexpr Gdz<tN> quiet_NaN() noexcept { return _CXType<Gdz<tN>>::Max(); }
		static constexpr Gdz<tN> signaling_NaN() noexcept { return _CXType<Gdz<tN>>::Max(); }
		static constexpr Gdz<tN> denorm_min() noexcept { return Gdz<tN>(0); }
		static constexpr bool is_iec559 = false;
		static constexpr bool is_modulo = false;
		static constexpr bool traps = std::numeric_limits<int32_t>::traps;
		static constexpr bool tinyness_before = false;
	};
};
#endif //DO_LIMITS

// global scope comparison operators
#define _OP(Op) template<typename T, int TN, typename std::enable_if_t<!std::is_same_v<Gdz<TN>, std::decay_t<T>> && (std::is_integral_v<T> || std::is_floating_point_v<T>), int> = 0> \
_CX bool operator Op (const T& L, const Gdz<TN>& R) { return L Op (const T)R; }
_OP(==) _OP(!=) _OP(<) _OP(>) _OP(<=) _OP(>=)
#undef _OP

// global bitwise and arithmetic operators - including assignment
#define _OP(Op) \
template<typename T, int TN, typename std::enable_if_t<!std::is_same_v<Gdz<TN>, std::decay_t<T>> && std::is_integral_v<T>, int> = 0> \
_CX T operator Op (const T& L, const Gdz<TN>& R) { return L Op static_cast<T>(R); } \
template<typename T, int TN, typename std::enable_if_t<!std::is_same_v<Gdz<TN>, std::decay_t<T>> && std::is_integral_v<T>, int> = 0> \
_CX T& operator Op##= (T& L, const Gdz<TN>& R) { L = L Op static_cast<T>(R); return L; }
_OP(&) _OP(|) _OP(^) _OP(+) _OP(-) _OP(*) _OP(/) _OP(%)
#undef _OP

template<int N> std::ostream& operator<<(std::ostream& os, const Gdz<N>& G) { return os << (STR) G; } // ostream integration

// this version has MSD partial block which is standard convention.  it writes each digit 2x and requires an extra digit buffer
template<int TN> // Flags: NumOrd:1 BlkOrd:2 Zeros:4 SignPos:8 Plus:16
STR _CXToStr(const Gdz<TN>& Num, U32 Radix, S32 DBlk, U32 Flags, const char* iSep, const char* iWrap, STR hOut) {
	_CX U32 BufMax = 512; DBlk = DBlk < 1 ? 1 : DBlk; Radix = min(max(Radix, (U32)2), (U32)64); // Blk > 0, Radix: 2-64, //if hOut, size>=256
	static char buf[BufMax * 2], DigBuf[BufMax]; static const char DigChr[] = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_@";
	S32 NDir = (Flags & 1) ? -1 : 1, BDir = (Flags & 2) ? -1 : 1, Dir = NDir * -1, DigCnt = 0, Blk, i;  // Number/Block direction: 0 for forward, 1 for backward

	Gdz<TN> Acc(+Num); // Collect digits into DigBuf
	for( *DigBuf = '0'; Acc != 0 && DigCnt < BufMax * 2; DigBuf[DigCnt++] = DigChr[Acc.DivV32(Radix, Acc)]); //DigChr[_CXDivV(Acc, Acc, Radix)]);
	if( DigCnt < 1 ) { DigBuf[0] = '0'; DigCnt = 1; }
	if( !( Flags & 4) ) DBlk = DBlk>DigCnt ? DigCnt : DBlk;
	//for(i=0;i<BufMax;++i) buf[i]='.';  // debugging 
	// Prepare wrap and sep for simple use - if NDir=-1, each must be reversed
	S32 PartCnt = DigCnt % DBlk, DigIdx, SepSz = 0, LWSz = 0, RWSz = 0;
	CHR Sep[16]{0}, LWrap[16]{0}, RWrap[16]{0}, * hC, * hCur; 
	for (hC = (STR)iSep; hC && hC[SepSz]; ++SepSz) Sep[SepSz] = hC[SepSz]; Sep[SepSz] = 0;
	for (hC = (STR)iWrap; hC && *hC > 7; ++LWSz) LWrap[LWSz] = *hC++; LWrap[i = LWSz] = 0;
	// add sign to LWrap or RWrap
	if (Num.sign() || (Flags & 16) ) //is there anything to write?
		if( ( NDir + !(Flags & 8) ) >0 ){ //sign left - ad it to LWrap and re terminate
			LWrap[LWSz++] = (Flags & 16) ? '+' : '-'; LWrap[LWSz] = 0;
		} else { // sign right - prepend it to RWrap
			RWrap[RWSz++] = (Flags & 16) ? '+' : '-';
		}
	for (++hC; hC && *hC; ++RWSz) RWrap[RWSz] = *hC++; RWrap[RWSz] = 0;
	if( Dir < 0 ) {// partial will be on left so process descending and flip wraps and sep
		for (i = 0; i < SepSz/2; ++i) _CXSwap(Sep[i], Sep[SepSz-1-i]);
		for (i = 0; i < LWSz/2;  ++i) _CXSwap(LWrap[i], LWrap[LWSz-1-i]);
		for (i = 0; i < RWSz/2;  ++i) _CXSwap(RWrap[i], RWrap[RWSz-1-i]);
	}

	Blk = DigCnt/DBlk + !!PartCnt -1;// Blk = Blk >= 0 ? Blk : 0;
	S32 TotalSz = LWSz + RWSz + DigCnt + SepSz * Blk + ((Flags & 4)&&PartCnt ? DBlk - PartCnt : 0);
	STR hX = hOut ? hOut + BufMax - 1 : buf + sizeof(buf) - 1, hO = hOut ? hOut : buf;

	if( ! hOut && ( TotalSz >= BufMax * 4 || DigCnt > BufMax * 2 ) ) { hC = DigBuf + min( (S32)BufMax - 5, DigCnt + 5);
		*hC++ = 'O'; *hC++ = 'V'; *hC++ = 'E'; *hC++ = 'R'; *hC = 0;  return DigBuf; } // over buffer, no hOut

	if( Dir < 0 ){ // set for descend, RWrap and terminate
		hO = (hOut ? hOut : buf) + TotalSz; 
		hX = hOut ? hOut : buf; hC = RWrap; *hO-- = 0; 
	} else {        // set for ascend, LWrap and terminate
		hO = hOut ? hOut : buf;                        
		hX = (hOut ? hOut : buf) + TotalSz; hC = LWrap; *hX-- = 0;
	}         
	for( hCur = hO; *hC; hCur += Dir ) *hCur = *hC++; // write start L|R wrap to hCur
	
	// place all digits into output fom lsd to msd, orient blocks and digits to NDir, BDir
	Blk = NDir * BDir > 0 ? 0: BDir * (DBlk-1);     // NB,-N-B = 0; -NB,N-B = Bdir*(DBlk-1)
	for( DigIdx = 0; DigIdx < DigCnt; ++DigIdx )
	{
		if (!(DigIdx % DBlk) && DigIdx ){  // don't add a leading or trailing sep
			for( hCur+=Dir*DBlk, hC=Sep; *hC; hCur+=Dir ) *hCur = *hC++;
			if( DigIdx >= DigCnt - PartCnt) 
				Blk += ( Flags & 4 || NDir*BDir > 0 ) ? 0 : NDir * (DBlk - PartCnt);  //add in part shift
		}
		hCur[ Blk - (BDir*DigIdx%DBlk) ] = DigBuf[DigIdx];                       //write the digit LSD to MSD order
	}

	if(PartCnt && (Flags & 4)) // write out required zeros to inner or outer pad
	do { hCur[ Blk - (BDir*DigIdx%DBlk)] = '0'; // no shift for zero fill
	} while( ++DigIdx < DigCnt + DBlk - PartCnt );


	// write the end L|R wrap to hCur
	hCur +=  Dir*(((Flags &4) || !PartCnt) ? DBlk : PartCnt);
	for( hC = Dir < 0 ? LWrap: RWrap; *hC; hCur += Dir ) *hCur = *hC++; //write start L|R wrap to hCur
	
	return Dir < 0 ? hX : hO;
}

// UDT for ""_g[radix] literals - binary, octal, decimal and hex
_CX GdzL operator"" _g2(const char* str, std::size_t len) {	
	GdzL RVal = 0; S32 skip = 0, sgn = (len > 0 && str[0] == '-') ? ++skip : 0, i = skip; U32 v = 0;
	for (; i < len; ++i) {
		v = (str[i] == '1') ? 1 : (str[i] == '0') ? 0 : ~0U;
		if (~0U == v) ++skip;
		else { RVal <<= 1; RVal.data[0] += v; }
	}
	return sgn ? RVal.neg() : RVal;
}
_CX GdzL operator"" _g8(const char* str, std::size_t len) { 
	GdzL RVal = 0; S32 skip = 0, sgn = (len > 0 && str[0] == '-') ? ++skip : 0, i = skip; U32 v = 0;
	for (; i < len; ++i) {
		v = (str[i] >= '0' && str[i] <= '7') ? (str[i] - '0') : ~0U;
		if (~0U == v) ++skip;
		else { RVal <<= 3; RVal.data[0] += v; }
	}
	return sgn ? RVal.neg() : RVal;
}
_CX GdzL operator"" _g10(const char* str, std::size_t len) { 
	GdzL RVal = 0; S32 skip = 0, sgn = (len > 0 && str[0] == '-') ? ++skip : 0, i = skip; U32 v = 0;
	for (; i < len; ++i) {
		v = (str[i] >= '0' && str[i] <= '9') ? (str[i] - '0') : ~0U;
		if (~0U == v) ++skip;
		else { RVal *= 10; RVal.data[0] += v; }
	}
	return sgn ? RVal.neg() : RVal;
}
_CX GdzL operator"" _g16(const char* str, std::size_t len) { 
	GdzL RVal = 0; S32 skip = 0, sgn = (len > 0 && str[0] == '-') ? ++skip : 0, i = skip; U32 v = 0;
	for (; i < len; ++i) {
		v = (str[i] >= '0' && str[i] <= '9') ? (str[i] - '0') :
			(str[i] >= 'A' && str[i] <= 'F') ? (10 + (str[i] - 'A')) :
			(str[i] >= 'a' && str[i] <= 'f') ? (10 + (str[i] - 'a')) : ~0U;
		if (~0U == v) ++skip;
		else { RVal <<= 4; RVal.data[0] += v; }
	}
	return sgn ? RVal.neg() : RVal;
}

mNotCreate(extern) GdzTxt<Gdz<16>> GdzT;
mNotCreate(extern) S8 gTblBit[256][4] mCreate( = { //0:cnt 1:lsb 2:msb 3:flp
{0, -1, -1, 0}, {1, 0, 0, -128}, {1, 1, 1, 64}, {2, 0, 1, -64}, {1, 2, 2, 32}, {2, 0, 2, -96}, {2, 1, 2, 96}, {3, 0, 2, -32},
{1, 3, 3, 16}, {2, 0, 3, -112}, {2, 1, 3, 80}, {3, 0, 3, -48}, {2, 2, 3, 48}, {3, 0, 3, -80}, {3, 1, 3, 112}, {4, 0, 3, -16},
{1, 4, 4, 8}, {2, 0, 4, -120}, {2, 1, 4, 72}, {3, 0, 4, -56}, {2, 2, 4, 40}, {3, 0, 4, -88}, {3, 1, 4, 104}, {4, 0, 4, -24},
{2, 3, 4, 24}, {3, 0, 4, -104}, {3, 1, 4, 88}, {4, 0, 4, -40}, {3, 2, 4, 56}, {4, 0, 4, -72}, {4, 1, 4, 120}, {5, 0, 4, -8},
{1, 5, 5, 4}, {2, 0, 5, -124}, {2, 1, 5, 68}, {3, 0, 5, -60}, {2, 2, 5, 36}, {3, 0, 5, -92}, {3, 1, 5, 100}, {4, 0, 5, -28},
{2, 3, 5, 20}, {3, 0, 5, -108}, {3, 1, 5, 84}, {4, 0, 5, -44}, {3, 2, 5, 52}, {4, 0, 5, -76}, {4, 1, 5, 116}, {5, 0, 5, -12},
{2, 4, 5, 12}, {3, 0, 5, -116}, {3, 1, 5, 76}, {4, 0, 5, -52}, {3, 2, 5, 44}, {4, 0, 5, -84}, {4, 1, 5, 108}, {5, 0, 5, -20},
{3, 3, 5, 28}, {4, 0, 5, -100}, {4, 1, 5, 92}, {5, 0, 5, -36}, {4, 2, 5, 60}, {5, 0, 5, -68}, {5, 1, 5, 124}, {6, 0, 5, -4},
{1, 6, 6, 2}, {2, 0, 6, -126}, {2, 1, 6, 66}, {3, 0, 6, -62}, {2, 2, 6, 34}, {3, 0, 6, -94}, {3, 1, 6, 98}, {4, 0, 6, -30},
{2, 3, 6, 18}, {3, 0, 6, -110}, {3, 1, 6, 82}, {4, 0, 6, -46}, {3, 2, 6, 50}, {4, 0, 6, -78}, {4, 1, 6, 114}, {5, 0, 6, -14},
{2, 4, 6, 10}, {3, 0, 6, -118}, {3, 1, 6, 74}, {4, 0, 6, -54}, {3, 2, 6, 42}, {4, 0, 6, -86}, {4, 1, 6, 106}, {5, 0, 6, -22},
{3, 3, 6, 26}, {4, 0, 6, -102}, {4, 1, 6, 90}, {5, 0, 6, -38}, {4, 2, 6, 58}, {5, 0, 6, -70}, {5, 1, 6, 122}, {6, 0, 6, -6},
{2, 5, 6, 6}, {3, 0, 6, -122}, {3, 1, 6, 70}, {4, 0, 6, -58}, {3, 2, 6, 38}, {4, 0, 6, -90}, {4, 1, 6, 102}, {5, 0, 6, -26},
{3, 3, 6, 22}, {4, 0, 6, -106}, {4, 1, 6, 86}, {5, 0, 6, -42}, {4, 2, 6, 54}, {5, 0, 6, -74}, {5, 1, 6, 118}, {6, 0, 6, -10},
{3, 4, 6, 14}, {4, 0, 6, -114}, {4, 1, 6, 78}, {5, 0, 6, -50}, {4, 2, 6, 46}, {5, 0, 6, -82}, {5, 1, 6, 110}, {6, 0, 6, -18},
{4, 3, 6, 30}, {5, 0, 6, -98}, {5, 1, 6, 94}, {6, 0, 6, -34}, {5, 2, 6, 62}, {6, 0, 6, -66}, {6, 1, 6, 126}, {7, 0, 6, -2},
{1, 7, 7, 1}, {2, 0, 7, -127}, {2, 1, 7, 65}, {3, 0, 7, -63}, {2, 2, 7, 33}, {3, 0, 7, -95}, {3, 1, 7, 97}, {4, 0, 7, -31},
{2, 3, 7, 17}, {3, 0, 7, -111}, {3, 1, 7, 81}, {4, 0, 7, -47}, {3, 2, 7, 49}, {4, 0, 7, -79}, {4, 1, 7, 113}, {5, 0, 7, -15},
{2, 4, 7, 9}, {3, 0, 7, -119}, {3, 1, 7, 73}, {4, 0, 7, -55}, {3, 2, 7, 41}, {4, 0, 7, -87}, {4, 1, 7, 105}, {5, 0, 7, -23},
{3, 3, 7, 25}, {4, 0, 7, -103}, {4, 1, 7, 89}, {5, 0, 7, -39}, {4, 2, 7, 57}, {5, 0, 7, -71}, {5, 1, 7, 121}, {6, 0, 7, -7},
{2, 5, 7, 5}, {3, 0, 7, -123}, {3, 1, 7, 69}, {4, 0, 7, -59}, {3, 2, 7, 37}, {4, 0, 7, -91}, {4, 1, 7, 101}, {5, 0, 7, -27},
{3, 3, 7, 21}, {4, 0, 7, -107}, {4, 1, 7, 85}, {5, 0, 7, -43}, {4, 2, 7, 53}, {5, 0, 7, -75}, {5, 1, 7, 117}, {6, 0, 7, -11},
{3, 4, 7, 13}, {4, 0, 7, -115}, {4, 1, 7, 77}, {5, 0, 7, -51}, {4, 2, 7, 45}, {5, 0, 7, -83}, {5, 1, 7, 109}, {6, 0, 7, -19},
{4, 3, 7, 29}, {5, 0, 7, -99}, {5, 1, 7, 93}, {6, 0, 7, -35}, {5, 2, 7, 61}, {6, 0, 7, -67}, {6, 1, 7, 125}, {7, 0, 7, -3},
{2, 6, 7, 3}, {3, 0, 7, -125}, {3, 1, 7, 67}, {4, 0, 7, -61}, {3, 2, 7, 35}, {4, 0, 7, -93}, {4, 1, 7, 99}, {5, 0, 7, -29},
{3, 3, 7, 19}, {4, 0, 7, -109}, {4, 1, 7, 83}, {5, 0, 7, -45}, {4, 2, 7, 51}, {5, 0, 7, -77}, {5, 1, 7, 115}, {6, 0, 7, -13},
{3, 4, 7, 11}, {4, 0, 7, -117}, {4, 1, 7, 75}, {5, 0, 7, -53}, {4, 2, 7, 43}, {5, 0, 7, -85}, {5, 1, 7, 107}, {6, 0, 7, -21},
{4, 3, 7, 27}, {5, 0, 7, -101}, {5, 1, 7, 91}, {6, 0, 7, -37}, {5, 2, 7, 59}, {6, 0, 7, -69}, {6, 1, 7, 123}, {7, 0, 7, -5},
{3, 5, 7, 7}, {4, 0, 7, -121}, {4, 1, 7, 71}, {5, 0, 7, -57}, {4, 2, 7, 39}, {5, 0, 7, -89}, {5, 1, 7, 103}, {6, 0, 7, -25},
{4, 3, 7, 23}, {5, 0, 7, -105}, {5, 1, 7, 87}, {6, 0, 7, -41}, {5, 2, 7, 55}, {6, 0, 7, -73}, {6, 1, 7, 119}, {7, 0, 7, -9},
{4, 4, 7, 15}, {5, 0, 7, -113}, {5, 1, 7, 79}, {6, 0, 7, -49}, {5, 2, 7, 47}, {6, 0, 7, -81}, {6, 1, 7, 111}, {7, 0, 7, -17},
{5, 3, 7, 31}, {6, 0, 7, -97}, {6, 1, 7, 95}, {7, 0, 7, -33}, {6, 2, 7, 63}, {7, 0, 7, -65}, {7, 1, 7, 127}, {8, 0, 7, -1}
} );

#endif //GDZ_HPP
