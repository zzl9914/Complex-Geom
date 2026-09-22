/**
 * @file scomplex.hpp
 * @brief Complex / split-complex / dual-number arithmetic (C++ templates).
 *
 * Lives in std::complex::experimental::zhangzl. The header is self-contained
 * via the standard library headers it includes; it does not include
 * <complex>, so it does not collide with class template std::complex.
 *
 * Algebra matches the original scomplex_no_throw.hpp. Unlike the canvas
 * no-throw build, this port actually throws on error instead of swallowing
 * exceptions.
 *
 * Non-template helpers are defined in scomplex.cpp. Functor templates
 * (differential, invert, integral) and the remaining templates must stay
 * in this header. Link scomplex.cpp when building.
 */

#pragma once
#ifndef SCOMPLEX_HPP_
#define SCOMPLEX_HPP_

#if defined(__GNUC__)||defined(__clang__)
    #define COMPLEX_ALWAYS_INLINE __attribute__((always_inline))
    #define COMPLEX_CONST __attribute__((const))
    #define COMPLEX_PURE __attribute__((pure))
    #define COMPLEX_LIKELY(x) __builtin_expect(!!(x),1)
    #define COMPLEX_UNLIKELY(x) __builtin_expect(!!(x),0)
    #define COMPLEX_RESTRICT __restrict__
#elif defined(_MSC_VER)
    #define COMPLEX_ALWAYS_INLINE __forceinline
    #define COMPLEX_CONST
    #define COMPLEX_PURE
    #define COMPLEX_LIKELY(x) (x)
    #define COMPLEX_UNLIKELY(x) (x)
    #define COMPLEX_RESTRICT __restrict
#else
    #define COMPLEX_ALWAYS_INLINE
    #define COMPLEX_CONST
    #define COMPLEX_PURE
    #define COMPLEX_LIKELY(x) (x)
    #define COMPLEX_UNLIKELY(x) (x)
    #define COMPLEX_RESTRICT
#endif

#if defined(_MSC_VER)
    #define COMPLEX_ALIGNAS(x) __declspec(align(x))
#else
    #define COMPLEX_ALIGNAS(x)
#endif

#if defined(__AVX__)||defined(__SSE2__)||defined(_M_AMD64)||defined(_M_X64)||defined(_M_IX86_FP)
#define COMPLEX_SIMD_ENABLED 1
#else
#define COMPLEX_SIMD_ENABLED 0
#endif

#if defined(__GNUC__)||defined(__clang__)
    #define COMPLEX_PREFETCH(addr,rw,locality) __builtin_prefetch(addr,rw,locality)
#elif defined(_MSC_VER)&&defined(_M_AMD64)
    #define COMPLEX_PREFETCH__GET_HINT(locality) ((locality)?(((locality)==1)?_MM_HINT_T1:(((locality)==2)?_MM_HINT_T2:_MM_HINT_NTA)):_MM_HINT_T0)
    #define COMPLEX_PREFETCH(addr,rw,locality) _mm_prefetch((const char*)(addr),COMPLEX_PREFETCH__GET_HINT(locality))
#else
    #define COMPLEX_PREFETCH(addr,rw,locality) ((void)0)
#endif

#if defined(__GNUC__) || defined(__clang__)
    #define COMPLEX_VECTORIZE _Pragma("GCC ivdep")
#elif defined(_MSC_VER)
    #define COMPLEX_VECTORIZE __pragma(loop(ivdep))
#else
    #define COMPLEX_VECTORIZE
#endif

#if defined(__GNUC__) || defined(__clang__)
    #define COMPLEX_EXPECT_BRANCH(cond, expected) __builtin_expect(!!(cond), expected)
#else
    #define COMPLEX_EXPECT_BRANCH(cond, expected) (cond)
#endif

#if defined(__GNUC__) || defined(__clang__)
    #define COMPLEX_HOT __attribute__((hot))
    #define COMPLEX_COLD __attribute__((cold))
#elif defined(_MSC_VER)
    #define COMPLEX_HOT __declspec(noinline)
    #define COMPLEX_COLD __declspec(noinline)
#else
    #define COMPLEX_HOT
    #define COMPLEX_COLD
#endif

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <ios>
#include <istream>
#include <limits>
#include <math.h>
#include <ostream>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#if defined(_MSC_VER)
    #include <intrin.h>
#else
    #include <emmintrin.h>
#endif

namespace std{
namespace complex{
namespace experimental{
namespace zhangzl{

/**
 * @brief Custom exception class for complex number operations
 */
struct complex_error:public runtime_error{
    /**
     * @brief Constructor for complex_error
     * @param what_arg Error message string
     */
    complex_error(const std::string& what_arg):std::runtime_error(what_arg){}
};

/**
 * @brief Exception for division by zero in complex numbers
 */
struct complex_division_by_zero:public complex_error{
    complex_division_by_zero():complex_error("Division by zero in complex number"){}
};

/**
 * @brief Exception for overflow in complex number operations
 */
struct complex_overflow:public complex_error{
    complex_overflow():complex_error("Overflow in complex number operation"){}
};

/**
 * @brief Exception for domain errors in complex number operations
 */
struct complex_domain_error:public complex_error{
    complex_domain_error(const std::string& m):complex_error(m){}
};

/**
 * @brief Exception for convergence errors in complex number operations
 */
struct complex_convergence_error:public complex_error{
    complex_convergence_error():complex_error("Failed to converge in complex operation"){}
};

/**
 * @brief Rounds a number to the left of the decimal point
 * @param a Number to round
 * @return Rounded number
 */
double round_left(double a);

/**
 * @brief Performs modular arithmetic on floating point numbers
 * @param a Number to be modified
 * @param b Modulus
 * @return Reference to modified number
 */
double& medium_mod(double& a,const double b);

/**
 * @brief Performs left shift operation on floating point numbers
 * @param a Number to shift
 * @param b Number of positions to shift
 * @return Shifted number
 */
double shift_left(double a,uint64_t b)noexcept;

/**
 * @brief Performs right shift operation on floating point numbers
 * @param a Number to shift
 * @param b Number of positions to shift
 * @return Shifted number
 */
double shift_right(double a,uint64_t b)noexcept;

/**
 * @brief Performs bitwise NOT operation on floating point numbers
 * @param a Number to operate on
 * @return Result of NOT operation
 */
template<typename _Tp>_Tp digital_not(_Tp a)noexcept{
    union{_Tp a;uint64_t ap;}_a;
    _a.a=a;
    _a.ap=~_a.ap;
    return _a.a;
}

/**
 * @brief Performs bitwise OR operation on floating point numbers
 * @param a First operand
 * @param b Second operand
 * @return Result of OR operation
 */
template<typename _Tp>_Tp digital_or(_Tp a,_Tp b)noexcept{
    union{_Tp a;uint64_t ap;}_a,_b,ans;
    _a.a=a;
    _b.a=b;
    ans.ap=_a.ap|_b.ap;
	if(!(ans.ap&((1<<52)-1)))return 0;
    return ans.a;
}

/**
 * @brief Performs bitwise AND operation on floating point numbers
 * @param a First operand
 * @param b Second operand
 * @return Result of AND operation
 */
template<typename _Tp>_Tp digital_and(_Tp a,_Tp b)noexcept{
    return digital_not(digital_or(digital_not(a),digital_not(b)));
}

/**
 * @brief Performs bitwise XOR operation on floating point numbers
 * @param a First operand
 * @param b Second operand
 * @return Result of XOR operation
 */
template<typename _Tp>_Tp digital_xor(_Tp a,_Tp b)noexcept{
    return digital_and(digital_or(a,b),digital_not(digital_and(a,b)));
}

/**
 * @brief Complex number template class with advanced mathematical operations
 * @tparam _Tp Type of the real and imaginary parts (default: double)
 */
template<typename _Tp=double>struct COMPLEX_ALIGNAS(16) alignas(16) complex
{
    typedef _Tp value_type;

    /**
     * @brief Sets the number system for complex operations
     * @param m System parameter (default: -1)
     */
    static void set_system(_Tp m=_Tp(-1))noexcept{_system=m;}

    /**
     * @brief Gets the current number system
     * @return Current system parameter
     */
    static _Tp get_system()noexcept{return _system;}

    /**
     * @brief Default constructor
     * @param real Real part (default: 0)
     * @param imag Imaginary part (default: 0)
     */
    constexpr complex(const _Tp& real=_Tp(),const _Tp& imag=_Tp()):_real(real),_imag(imag){
        check_overflow(_real);
        check_overflow(_imag);
    }

    /**
     * @brief Copy constructor
     * @param a Complex number to copy
     */
    constexpr complex(const complex& a):_real(a._real),_imag(a._imag){
        check_overflow(_real);
        check_overflow(_imag);
    }

    /**
	 * @brief Improved string constructor with better error handling
	 * @param str String representation of complex number
	 */
	constexpr complex(const std::string& str){
	    std::regex pattern(R"(^\s*([+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)\s*([+-])\s*([+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)?i\s*$|^\s*([+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)\s*$|^\s*([+-])?i\s*$)");
	    std::smatch matches;
	    
	    if(!std::regex_match(str,matches,pattern))throw std::invalid_argument("Invalid complex number format: "+str);
	    _real = _Tp();
	    _imag = _Tp();
	    // Case 1: a+bi or a-bi
	    if(matches[1].matched){
	        _real=std::stod(matches[1].str());
	        std::string imag_sign=matches[2].str();
	        if(matches[3].matched)_imag=std::stod(matches[3].str());
			else _imag=_Tp(1);
	        if(imag_sign=="-")_imag = -_imag;
	    }
	    // Case 2: real only
	    else if(matches[4].matched)_real=std::stod(matches[4].str());
	    // Case 3: imaginary only
	    else if(matches[5].matched||(matches[0].str().find('i')!=std::string::npos)){
			std::string imag_part=matches[5].matched?matches[5].str():"";
	        if(imag_part=="-")_imag=_Tp(-1);
			else _imag=_Tp(1);
	    }
	    check_overflow(_real);
	    check_overflow(_imag);
	}

    /**
	 * @brief Move constructor
	 * @param other Complex number to move from
	 */
	constexpr complex(complex&& other) noexcept
	    : _real(std::move(other._real)), _imag(std::move(other._imag)) {
	    other._real = _Tp();
	    other._imag = _Tp();
	}

    /**
     * @brief Assignment operator
     * @param a Complex number to assign
     * @return Reference to this complex number
     */
    complex& operator=(const complex& a){
        if(this!=&a){
            complex temp=a;
            _Tp tempr=temp._real,tempi=temp._imag;
            temp._real=_real;
            temp._imag=_imag;
            _real=tempr;
            _imag=tempi;
        }
        return *this;
    }

	/**
	 * @brief Move assignment operator
	 * @param other Complex number to move from
	 * @return Reference to this complex number
	 */
	complex& operator=(complex&& other) noexcept {
	    if (this != &other) {
	        _real = std::move(other._real);
	        _imag = std::move(other._imag);
	        other._real = _Tp();
	        other._imag = _Tp();
	    }
	    return *this;
	}

    /**
     * @brief Subscript operator
     * @param n Index (0 for real, 1 for imaginary)
     * @return Reference to the selected component
     */
    COMPLEX_CONST _Tp& operator[](bool n)noexcept{return n?_imag:_real;}

    /**
     * @brief Const subscript operator
     * @param n Index (0 for real, 1 for imaginary)
     * @return Const reference to the selected component
     */
    COMPLEX_CONST const _Tp& operator[](bool n)const noexcept{return n?_imag:_real;}

	/**
	 * @brief Gets the real part
	 * @return Real part
	 */
	COMPLEX_CONST _Tp real()const noexcept{return _real;}

	/**
	 * @brief Gets the imaginary part
	 * @return Imaginary part
	 */
	COMPLEX_CONST _Tp imag()const noexcept{return _imag;}

	/**
	 * @brief Sets the real part
	 * @param r New real value
	 */
	void real(_Tp r){
	    _real=r;
	    check_overflow(r);
	}

	/**
	 * @brief Sets the imaginary part
	 * @param i New imaginary value
	 */
	void imag(_Tp i){
	    _imag=i;
	    check_overflow(i);
	}

    /**
     * @brief Boolean conversion operator
     * @return True if either component is non-zero
     */
    explicit operator bool()const noexcept{return _real||_imag;}

    /**
     * @brief Unary plus operator
     * @return Copy of this complex number
     */
    COMPLEX_ALWAYS_INLINE constexpr complex operator+()const noexcept{return *this;}

    /**
     * @brief Unary minus operator
     * @return Negated complex number
     */
    COMPLEX_ALWAYS_INLINE constexpr complex operator-()const noexcept{return complex(-_real,-_imag);}

    /**
     * @brief Equality comparison with real number
     * @param rhs Real number to compare
     * @return True if equal
     */
    COMPLEX_ALWAYS_INLINE constexpr bool operator==(const _Tp rhs)const noexcept{return _real==rhs&&!_imag;}

    /**
     * @brief Inequality comparison with real number
     * @param rhs Real number to compare
     * @return True if not equal
     */
    COMPLEX_ALWAYS_INLINE constexpr bool operator!=(const _Tp rhs)const noexcept{return !(*this==rhs);}

    /**
     * @brief Friend equality operator
     * @param lhs Real number
     * @param rhs Complex number
     * @return True if equal
     */
    friend inline constexpr bool operator==(_Tp lhs,complex rhs)noexcept{return rhs==lhs;}

    /**
     * @brief Friend inequality operator
     * @param lhs Real number
     * @param rhs Complex number
     * @return True if not equal
     */
    friend inline constexpr bool operator!=(_Tp lhs,complex rhs)noexcept{return rhs!=lhs;}

    /**
     * @brief Equality comparison with complex number
     * @param rhs Complex number to compare
     * @return True if equal
     */
    COMPLEX_ALWAYS_INLINE constexpr bool operator==(const complex rhs)const noexcept{
        return _real==rhs._real&&_imag==rhs._imag;
    }

    /**
     * @brief Inequality comparison with complex number
     * @param rhs Complex number to compare
     * @return True if not equal
     */
    COMPLEX_ALWAYS_INLINE constexpr bool operator!=(const complex rhs)const noexcept{
        return !(*this==rhs);
    }

    /**
     * @brief Less than comparison
     * @param rhs Complex number to compare
     * @return True if less
     */
    COMPLEX_ALWAYS_INLINE constexpr bool operator<(const complex rhs)const noexcept{
        return (_real==rhs._real)?(_imag<rhs._imag):(_real<rhs._real);
    }

    /**
     * @brief Greater than comparison
     * @param rhs Complex number to compare
     * @return True if greater
     */
    COMPLEX_ALWAYS_INLINE bool operator>(const complex rhs)const noexcept{
        return !(*this<=rhs);
    }

    /**
     * @brief Less than or equal comparison
     * @param rhs Complex number to compare
     * @return True if less or equal
     */
    COMPLEX_ALWAYS_INLINE bool operator<=(const complex rhs)const noexcept{
        return (*this<rhs)||(*this==rhs);
    }

    /**
     * @brief Greater than or equal comparison
     * @param rhs Complex number to compare
     * @return True if greater or equal
     */
    COMPLEX_ALWAYS_INLINE bool operator>=(const complex rhs)const noexcept{
        return !(*this<rhs);
    }

    /**
     * @brief Pre-increment operator
     * @return Reference to incremented complex number
     */
    COMPLEX_ALWAYS_INLINE complex& operator++(){
        ++_real;
        check_overflow(_real);
        return *this;
    }

    /**
     * @brief Pre-decrement operator
     * @return Reference to decremented complex number
     */
    COMPLEX_ALWAYS_INLINE complex& operator--(){
        --_real;
        check_overflow(_real);
        return *this;
    }

    /**
     * @brief Post-increment operator
     * @return Copy of complex number before increment
     */
    COMPLEX_ALWAYS_INLINE constexpr complex operator++(int){
        complex b=*this;
        ++*this;
        return b;
    }

    /**
     * @brief Post-decrement operator
     * @return Copy of complex number before decrement
     */
    COMPLEX_ALWAYS_INLINE constexpr complex operator--(int){
        complex b=*this;
        --*this;
        return b;
    }

    /**
     * @brief Bitwise NOT operator
     * @return Complex number with bitwise NOT applied to components
     */
    COMPLEX_ALWAYS_INLINE constexpr complex operator~()const noexcept{
        return complex(digital_not(_real),digital_not(_imag));
    }

    /**
     * @brief Addition assignment with real number
     * @param rhs Real number to add
     * @return Reference to this complex number
     */
    COMPLEX_ALWAYS_INLINE complex& operator+=(const _Tp rhs){
        _real+=rhs;
        check_overflow(_real);
        return *this;
    }

    /**
     * @brief Subtraction assignment with real number
     * @param rhs Real number to subtract
     * @return Reference to this complex number
     */
    COMPLEX_ALWAYS_INLINE complex& operator-=(const _Tp rhs){
        _real-=rhs;
        check_overflow(_real);
        return *this;
    }

    /**
     * @brief Multiplication assignment with real number
     * @param rhs Real number to multiply
     * @return Reference to this complex number
     */
    COMPLEX_ALWAYS_INLINE complex& operator*=(const _Tp rhs){
        _real*=rhs;
        _imag*=rhs;
        check_overflow(_real);
        check_overflow(_imag);
        return *this;
    }

    /**
     * @brief Division assignment with real number
     * @param rhs Real number to divide by
     * @return Reference to this complex number
     */
    COMPLEX_ALWAYS_INLINE complex& operator/=(const _Tp rhs){
        check_division_by_zero(rhs);
        _real/=rhs;
        _imag/=rhs;
        return *this;
    }

    /**
     * @brief Modulo assignment with real number
     * @param rhs Real number to modulo by
     * @return Reference to this complex number
     */
    COMPLEX_ALWAYS_INLINE complex& operator%=(const _Tp rhs){
        check_division_by_zero(rhs);
        _real=medium_mod(_real,rhs);
        _imag=medium_mod(_imag,rhs);
        return *this;
    }

    /**
     * @brief Addition with real number
     * @param rhs Real number to add
     * @return Result of addition
     */
    COMPLEX_ALWAYS_INLINE inline constexpr complex operator+(const _Tp rhs){
        complex b=*this;
        b+=rhs;
        return b;
    }

    /**
     * @brief Subtraction with real number
     * @param rhs Real number to subtract
     * @return Result of subtraction
     */
    COMPLEX_ALWAYS_INLINE inline constexpr complex operator-(const _Tp rhs){
        complex b=*this;
        b-=rhs;
        return b;
    }

    /**
     * @brief Multiplication with real number
     * @param rhs Real number to multiply
     * @return Result of multiplication
     */
    COMPLEX_ALWAYS_INLINE inline constexpr complex operator*(const _Tp rhs){
        complex b=*this;
        b*=rhs;
        return b;
    }

    /**
     * @brief Division with real number
     * @param rhs Real number to divide by
     * @return Result of division
     */
    COMPLEX_ALWAYS_INLINE inline constexpr complex operator/(const _Tp rhs){
        complex b=*this;
        b/=rhs;
        return b;
    }

    /**
     * @brief Modulo with real number
     * @param rhs Real number to modulo by
     * @return Result of modulo
     */
    COMPLEX_ALWAYS_INLINE inline constexpr complex operator%(const _Tp rhs){
        complex b=*this;
        b%=rhs;
        return b;
    }

    /**
     * @brief Friend addition operator
     * @param a Real number
     * @param b Complex number
     * @return Result of addition
     */
    friend inline constexpr complex operator+(_Tp a,complex b){return b+a;}

    /**
     * @brief Friend subtraction operator
     * @param a Real number
     * @param b Complex number
     * @return Result of subtraction
     */
    friend inline constexpr complex operator-(_Tp a,complex b){return -b+a;}

    /**
     * @brief Friend multiplication operator
     * @param a Real number
     * @param b Complex number
     * @return Result of multiplication
     */
    friend inline constexpr complex operator*(_Tp a,complex b){return b*a;}

    /**
     * @brief Friend division operator
     * @param a Real number
     * @param b Complex number
     * @return Result of division
     */
    friend inline constexpr complex operator/(_Tp a,complex b){return complex(a)/b;}

    /**
     * @brief Friend modulo operator
     * @param a Real number
     * @param b Complex number
     * @return Result of modulo
     */
    friend inline constexpr complex operator%(_Tp a,complex b){return complex(a)%b;}

    /**
     * @brief Addition assignment with complex number
     * @param rhs Complex number to add
     * @return Reference to this complex number
     */
    COMPLEX_ALWAYS_INLINE complex& operator+=(const complex rhs){
        _real+=rhs._real;
        _imag+=rhs._imag;
        check_overflow(_real);
        check_overflow(_imag);
        return *this;
    }

    /**
     * @brief Subtraction assignment with complex number
     * @param rhs Complex number to subtract
     * @return Reference to this complex number
     */
    COMPLEX_ALWAYS_INLINE complex& operator-=(const complex rhs){
        _real-=rhs._real;
        _imag-=rhs._imag;
        check_overflow(_real);
        check_overflow(_imag);
        return *this;
    }

    /**
     * @brief Multiplication assignment with complex number
     * @param rhs Complex number to multiply
     * @return Reference to this complex number
     */
    COMPLEX_ALWAYS_INLINE complex& operator*=(const complex rhs){
        _Tp _temp=_real*rhs._imag+_imag*rhs._real;
        _real=_real*rhs._real+_system*_imag*rhs._imag;
		_imag=_temp;
        check_overflow(_real);
        check_overflow(_imag);
        return *this;
    }

    /**
     * @brief Division assignment with complex number
     * @param rhs Complex number to divide by
     * @return Reference to this complex number
     */
    COMPLEX_ALWAYS_INLINE complex& operator/=(const complex rhs){
    	_Tp _temp=(_imag*rhs._real-_real*rhs._imag);
        _Tp denom=rhs._real*rhs._real-_system*rhs._imag*rhs._imag;
        check_division_by_zero(denom);
        _real=(_real*rhs._real-_system*_imag*rhs._imag)/denom;
		_imag=_temp/denom;
        check_overflow(_real);
        check_overflow(_imag);
        return *this;
    }

    /**
     * @brief Modulo assignment with complex number
     * @param rhs Complex number to modulo by
     * @return Reference to this complex number
     */
    COMPLEX_ALWAYS_INLINE complex& operator%=(const complex rhs){
        *this/=rhs;
        _real=round_left(_real);
        _imag=round_left(_imag);
        *this*=rhs;
        return *this;
    }

    /**
     * @brief Addition with complex number
     * @param rhs Complex number to add
     * @return Result of addition
     */
    COMPLEX_ALWAYS_INLINE inline constexpr complex operator+(const complex rhs){
        complex b=*this;
        b+=rhs;
        return b;
    }

    /**
     * @brief Subtraction with complex number
     * @param rhs Complex number to subtract
     * @return Result of subtraction
     */
    COMPLEX_ALWAYS_INLINE inline constexpr complex operator-(const complex rhs){
        complex b=*this;
        b-=rhs;
        return b;
    }

    /**
     * @brief Multiplication with complex number
     * @param rhs Complex number to multiply
     * @return Result of multiplication
     */
    COMPLEX_ALWAYS_INLINE inline constexpr complex operator*(const complex rhs){
        complex b=*this;
        b*=rhs;
        return b;
    }

    /**
     * @brief Division with complex number
     * @param rhs Complex number to divide by
     * @return Result of division
     */
    COMPLEX_ALWAYS_INLINE inline constexpr complex operator/(const complex rhs){
        complex b=*this;
        b/=rhs;
        return b;
    }

    /**
     * @brief Modulo with complex number
     * @param rhs Complex number to modulo by
     * @return Result of modulo
     */
    COMPLEX_ALWAYS_INLINE inline constexpr complex operator%(const complex rhs){
        complex b=*this;
        b%=rhs;
        return b;
    }

    /**
     * @brief Left shift assignment
     * @param rhs Number of positions to shift
     * @return Reference to this complex number
     */
    COMPLEX_ALWAYS_INLINE complex& operator<<=(const int rhs){
        _real=shift_left(_real,rhs);
        _imag=shift_left(_imag,rhs);
        check_overflow(_real);
        check_overflow(_imag);
        return *this;
    }

    /**
     * @brief Right shift assignment
     * @param rhs Number of positions to shift
     * @return Reference to this complex number
     */
    COMPLEX_ALWAYS_INLINE complex& operator>>=(const int rhs){
        _real=shift_right(_real,rhs);
        _imag=shift_right(_imag,rhs);
        check_overflow(_real);
        check_overflow(_imag);
        return *this;
    }

    /**
     * @brief Bitwise AND assignment with real number
     * @param rhs Real number to AND with
     * @return Reference to this complex number
     */
    COMPLEX_ALWAYS_INLINE complex& operator&=(const _Tp rhs)noexcept{
        _real=digital_and(_real,rhs);
        _imag=digital_and(_imag,0);
        return *this;
    }

    /**
     * @brief Bitwise OR assignment with real number
     * @param rhs Real number to OR with
     * @return Reference to this complex number
     */
    COMPLEX_ALWAYS_INLINE complex& operator|=(const _Tp rhs)noexcept{
        _real=digital_or(_real,rhs);
        _imag=digital_or(_imag,0);
        return *this;
    }

    /**
     * @brief Bitwise XOR assignment with real number
     * @param rhs Real number to XOR with
     * @return Reference to this complex number
     */
    COMPLEX_ALWAYS_INLINE complex& operator^=(const _Tp rhs)noexcept{
        _real=digital_xor(_real,rhs);
        _imag=digital_xor(_imag,0);
        return *this;
    }

    /**
     * @brief Left shift operator
     * @param rhs Number of positions to shift
     * @return Result of shift
     */
    COMPLEX_ALWAYS_INLINE constexpr complex operator<<(const int rhs){
        complex b=*this;
        b<<=rhs;
        return b;
    }

    /**
     * @brief Right shift operator
     * @param rhs Number of positions to shift
     * @return Result of shift
     */
    COMPLEX_ALWAYS_INLINE constexpr complex operator>>(const int rhs){
        complex b=*this;
        b>>=rhs;
        return b;
    }

    /**
     * @brief Bitwise AND operator with real number
     * @param rhs Real number to AND with
     * @return Result of AND
     */
    COMPLEX_ALWAYS_INLINE constexpr complex operator&(const _Tp rhs){
        complex b=*this;
        b&=rhs;
        return b;
    }

    /**
     * @brief Bitwise OR operator with real number
     * @param rhs Real number to OR with
     * @return Result of OR
     */
    COMPLEX_ALWAYS_INLINE constexpr complex operator|(const _Tp rhs){
        complex b=*this;
        b|=rhs;
        return b;
    }

    /**
     * @brief Bitwise XOR operator with real number
     * @param rhs Real number to XOR with
     * @return Result of XOR
     */
    COMPLEX_ALWAYS_INLINE constexpr complex operator^(const _Tp rhs){
        complex b=*this;
        b^=rhs;
        return b;
    }

    /**
     * @brief Friend bitwise AND operator
     * @param a Real number
     * @param b Complex number
     * @return Result of AND
     */
    friend inline constexpr complex operator&(_Tp a,complex b){return b&a;}

    /**
     * @brief Friend bitwise OR operator
     * @param a Real number
     * @param b Complex number
     * @return Result of OR
     */
    friend inline constexpr complex operator|(_Tp a,complex b){return b|a;}

    /**
     * @brief Friend bitwise XOR operator
     * @param a Real number
     * @param b Complex number
     * @return Result of XOR
     */
    friend inline constexpr complex operator^(_Tp a,complex b){return b^a;}

    /**
     * @brief Bitwise AND assignment with complex number
     * @param rhs Complex number to AND with
     * @return Reference to this complex number
     */
    COMPLEX_ALWAYS_INLINE complex& operator&=(const complex rhs){
        _real=digital_and(_real,rhs._real);
        _imag=digital_and(_imag,rhs._imag);
        return *this;
    }

    /**
     * @brief Bitwise OR assignment with complex number
     * @param rhs Complex number to OR with
     * @return Reference to this complex number
     */
    COMPLEX_ALWAYS_INLINE complex& operator|=(const complex rhs){
        _real=digital_or(_real,rhs._real);
        _imag=digital_or(_imag,rhs._imag);
        return *this;
    }

    /**
     * @brief Bitwise XOR assignment with complex number
     * @param rhs Complex number to XOR with
     * @return Reference to this complex number
     */
    COMPLEX_ALWAYS_INLINE complex& operator^=(const complex rhs){
        _real=digital_xor(_real,rhs._real);
        _imag=digital_xor(_imag,rhs._imag);
        return *this;
    }

    /**
     * @brief Bitwise AND operator with complex number
     * @param rhs Complex number to AND with
     * @return Result of AND
     */
    COMPLEX_ALWAYS_INLINE constexpr complex operator&(const complex rhs){
        complex b=*this;
        b&=rhs;
        return b;
    }

    /**
     * @brief Bitwise OR operator with complex number
     * @param rhs Complex number to OR with
     * @return Result of OR
     */
    COMPLEX_ALWAYS_INLINE constexpr complex operator|(const complex rhs){
        complex b=*this;
        b|=rhs;
        return b;
    }

    /**
     * @brief Bitwise XOR operator with complex number
     * @param rhs Complex number to XOR with
     * @return Result of XOR
     */
    COMPLEX_ALWAYS_INLINE constexpr complex operator^(const complex rhs){
        complex b=*this;
        b^=rhs;
        return b;
    }

    /**
     * @brief Prefetches complex number data
     * @param addr Address of data to prefetch
     */
    static void prefetch(const complex *addr)noexcept{COMPLEX_PREFETCH(addr,0,3);}

    /**
     * @brief Applies function to array of complex numbers
     * @param data Array of complex numbers
     * @param n Size of array
     * @param func Function to apply
     */
    template<typename _Func>static void for_each(complex *data,size_t n,const _Func func){
        for(size_t i=0;i<n;++i){
            if(i+4<n)prefetch(&data[i+4]);
            data[i]=func(data[i]);
        }
    }

	static bool print_number;///< Iostream style

    static constexpr _Tp zero()noexcept{return _Tp();}
    static constexpr _Tp one()noexcept{return _Tp(1);}
    static constexpr _Tp sqrt2()noexcept{return _Tp(1.41421356237309504880);}
    static constexpr _Tp sqrt1_2()noexcept{return _Tp(.70710678118654752440);}
    static constexpr _Tp cbrt2()noexcept{return _Tp(1.25992104989487316476);}
    static constexpr _Tp cbrt1_2()noexcept{return _Tp(.79370052598409973737);}
    static constexpr _Tp pi()noexcept{return _Tp(3.14159265358979323846);}
    static constexpr _Tp pi2()noexcept{return _Tp(6.28318530717958647692);}
    static constexpr _Tp pi1_2()noexcept{return _Tp(1.57079632679489661923);}
    static constexpr _Tp pi1_3()noexcept{return _Tp(1.04719755119659774615);}
    static constexpr _Tp e()noexcept{return _Tp(2.71828182845904523536);}
    static constexpr complex i()noexcept{return complex(_Tp(0),_Tp(1));}
    static constexpr complex exp_i_pi_1_6()noexcept{return complex(_Tp(.5),_Tp(.86602540378443864676));}
    static constexpr complex exp_i_pi_1_9()noexcept{return complex(_Tp(.93969262078590838405),_Tp(.34202014332566873304));}

private:
    struct{
        union{
            struct{_Tp _real,_imag;};  ///< Real and imaginary parts
            __m128d simd_complex;      ///< SIMD representation
        };
    };
    static _Tp _system;  ///< Number system parameter

	/**
	 * @brief Checks if the complex number is finite
	 * @return True if both components are finite
	 */
	bool is_finite()const noexcept{return std::isfinite(_real)&&std::isfinite(_imag);}

	/**
	 * @brief Checks if the complex number is NaN
	 * @return True if either component is NaN
	 */
	bool is_nan()const noexcept{return std::isnan(_real)||std::isnan(_imag);}

	/**
	 * @brief Checks if the complex number is zero
	 * @return True if both components are zero
	 */
	bool is_zero()const noexcept{return _real==_Tp()&&_imag==_Tp();}

    /**
     * @brief Checks for overflow in complex number operations
     * @param v Value to check
     */
    void check_overflow(_Tp v)const{
    	if(!std::isfinite(v))throw complex_overflow();
        if(v<0)v=-v;
        if(v>std::numeric_limits<_Tp>::max()/2)throw complex_overflow();
    }

    /**
     * @brief Checks for division by zero
     * @param d Divisor to check
     */
    void check_division_by_zero(_Tp d)const{
        if(d<0)d=-d;
        if (d<std::numeric_limits<_Tp>::epsilon())throw complex_division_by_zero();
    }

    /**
     * @brief Creates backup of complex number
     * @return Backup copy
     */
    complex create_backup()const noexcept{return complex(_real,_imag);}

    /**
     * @brief Restores complex number from backup
     * @param backup Backup to restore from
     */
    void restore_backup(const complex& backup){
        _real=backup._real;
        _imag=backup._imag;
    }
};

// Initialize static member
template<typename _Tp>_Tp complex<_Tp>::_system=_Tp(-1);
template<typename _Tp>bool complex<_Tp>::print_number=0;

/**
 * @brief Gets real part of complex number
 * @param a Complex number
 * @return Real part
 */
template<typename _Tp>inline constexpr _Tp real(complex<_Tp> a)noexcept{return a[0];}

/**
 * @brief Gets imaginary part of complex number
 * @param a Complex number
 * @return Imaginary part
 */
template<typename _Tp>inline constexpr _Tp imag(complex<_Tp> a)noexcept{return a[1];}

/**
 * @brief Calculates norm of complex number
 * @param a Complex number
 * @return Norm value
 */
template<typename _Tp>inline constexpr _Tp norm(complex<_Tp> a)noexcept{
    return real(a)*real(a)-complex<_Tp>::get_system()*imag(a)*imag(a);
}

/**
 * @brief Calculates squared norm of complex number
 * @param a Complex number
 * @return Squared norm value
 */
template<typename _Tp>inline constexpr _Tp norm_norm(complex<_Tp> a)noexcept{
    return real(a)*real(a)+imag(a)*imag(a);
}

/**
 * @brief Calculates complex conjugate
 * @param a Complex number
 * @return Complex conjugate
 */
template<typename _Tp>inline constexpr complex<_Tp> conj(complex<_Tp> a)noexcept{
    return complex<_Tp>(real(a),-imag(a));
}

/**
 * @brief Calculates absolute value
 * @param a Complex number
 * @return Absolute value
 */
template<typename _Tp>inline constexpr complex<_Tp> abs(complex<_Tp> a)noexcept{
    _Tp _real=real(a),_imag=imag(a);
    if(_real<0)_real=-_real;
    if(_imag<0)_imag=-_imag;
    return complex<_Tp>(_real,_imag);
}

/**
 * @brief Raises complex number to integer power
 * @param base Base complex number
 * @param exponent Integer exponent
 * @return Result of power operation
 */
template<typename _Tp>complex<_Tp> pow(complex<_Tp> base,int exponent){
    if(!exponent)return 1;
    if(!base)return 0;
    if(exponent<0)return 1/pow(base,-exponent);
    if(exponent%2)return base*pow(base,exponent-1);
    complex<_Tp> ans=pow(base,exponent/2);
    return ans*ans;
}

/**
 * @brief Calculates square root of complex number
 * @param a Complex number
 * @param max_iterations Maximum iterations for convergence
 * @return Square root
 */
template<typename _Tp>complex<_Tp> sqrt(complex<_Tp> a,int max_iterations){
    if(!imag(a)&&real(a)>=0)return std::sqrt(real(a));
    _Tp L=norm(a);
	if(L<0)throw domain_error("Square root of non-positive number");
	L=norm_norm(a);
    if(L>2)return sqrt(a>>1,max_iterations)*complex<_Tp>::sqrt2();
    if(L<.5)return sqrt(a<<1,max_iterations)*complex<_Tp>::sqrt1_2();
    auto sqrt_with_init=[](complex<_Tp> a,complex<_Tp> b,int max_iterations)->complex<_Tp>{
        complex<_Tp> bp=b;
        for(int i=0;i<max_iterations;++i){
            try{
                b=(bp+a/bp)/2;
                complex<_Tp> d(real(b)-real(bp),imag(b)-imag(bp));
                d=abs(d);
                if(real(d)<std::numeric_limits<_Tp>::epsilon()&&imag(d)<std::numeric_limits<_Tp>::epsilon())break;
            }catch(const complex_error&){throw complex_convergence_error();}
            bp=b;
        }
        return b;
    };
    if(real(a)>0)return sqrt_with_init(a,1,max_iterations);
    else if(imag(a)>0)return sqrt_with_init(a,complex<_Tp>::exp_i_pi_1_6(),max_iterations);
    else return sqrt_with_init(a,conj(complex<_Tp>::exp_i_pi_1_6()),max_iterations);
}

/**
 * @brief Calculates square root with default iterations
 * @param a Complex number
 * @return Square root
 */
template<typename _Tp>complex<_Tp> sqrt(complex<_Tp> a){return sqrt(a,100);}

/**
 * @brief Calculates cube root of complex number
 * @param a Complex number
 * @param max_iterations Maximum iterations for convergence
 * @return Cube root
 */
template<typename _Tp>complex<_Tp> cbrt(complex<_Tp> a,int max_iterations){
    if(!imag(a))return std::cbrt(real(a));
    _Tp L=norm_norm(a);
    if(L>2)return cbrt(a>>1,max_iterations)*complex<_Tp>::cbrt2();
    if(L<.5)return cbrt(a<<1,max_iterations)*complex<_Tp>::cbrt1_2();
    auto cbrt_with_init=[](complex<_Tp> a,complex<_Tp> b,int max_iterations)->complex<_Tp>{
        complex<_Tp> bp=b;
        for(int i=0;i<max_iterations;++i){
            try{
                b=(bp*2+a/(bp*bp))/3;
                complex<_Tp> d(real(b)-real(bp),imag(b)-imag(bp));
                d=abs(d);
                if(real(d)<std::numeric_limits<_Tp>::epsilon()&&imag(d)<std::numeric_limits<_Tp>::epsilon())break;
            }catch(const complex_error&){throw complex_convergence_error();}
            bp=b;
        }
        return b;
    };
    if(real(a)>0)return cbrt_with_init(a,1,max_iterations);
    else if(imag(a)>0)return cbrt_with_init(a,complex<_Tp>::exp_i_pi_1_9(),max_iterations);
    else return cbrt_with_init(a,conj(complex<_Tp>::exp_i_pi_1_9()),max_iterations);
}

/**
 * @brief Calculates cube root with default iterations
 * @param a Complex number
 * @return Cube root
 */
template<typename _Tp>complex<_Tp> cbrt(complex<_Tp> a){return cbrt(a,100);}

/**
 * @brief Calculates nth root of complex number
 * @param a Complex number
 * @param n Root degree
 * @param max_iterations Maximum iterations for convergence
 * @return nth root
 */
template<typename _Tp>complex<_Tp> root(complex<_Tp> a,int n,int max_iterations){
    if(n==1)return a;
    if(n==2)return sqrt(a);
    if(n==3)return cbrt(a);
    if(!imag(a))return std::pow(real(a),1.0/n);
    _Tp L=norm_norm(a);
    if(L>2)return root(a>>1,n,max_iterations)*_Tp(std::pow(2.0,1.0/n));
    if(L<.5)return root(a<<1,n,max_iterations)*_Tp(std::pow(0.5,1.0/n));
    const complex<_Tp> exp_i_pi_1_3_n(_Tp(std::cos(complex<_Tp>::pi1_3()/n)),_Tp(std::sin(complex<_Tp>::pi1_3()/n)));
    auto root_with_init=[](complex<_Tp> a,complex<_Tp> b,int n,int max_iterations)->complex<_Tp>{
        complex<_Tp> bp=b;
        for(int i=0;i<max_iterations;++i){
            try{
                b=(bp*(n-1)+a/pow(bp,n-1))/n;
                complex<_Tp> d(real(b)-real(bp),imag(b)-imag(bp));
                d=abs(d);
                if(real(d)<std::numeric_limits<_Tp>::epsilon()&&imag(d)<std::numeric_limits<_Tp>::epsilon())break;
            }catch(const complex_error&){throw complex_convergence_error();}
            bp=b;
        }
        return b;
    };
    if(real(a)>0)return root_with_init(a,1,n,max_iterations);
    else if(imag(a)>0)return root_with_init(a,exp_i_pi_1_3_n,n,max_iterations);
    else return root_with_init(a,conj(exp_i_pi_1_3_n),n,max_iterations);
}

/**
 * @brief Calculates nth root with default iterations
 * @param a Complex number
 * @param n Root degree
 * @return nth root
 */
template<typename _Tp>complex<_Tp> root(complex<_Tp> a,int n){return root(a,n,100);}

/**
 * @brief Calculates exponential of complex number
 * @param a Complex number
 * @param max_iterations Maximum iterations for convergence
 * @return Exponential value
 */
template<typename _Tp>complex<_Tp> exp(complex<_Tp> a,int max_iterations){
    if(!imag(a))return std::exp(real(a));
    _Tp L=norm_norm(a);
    if(L>1){
        complex<_Tp> b=exp(a>>1,max_iterations);
        return b*b;
    }
    complex<_Tp> b=1,_temp=complex<_Tp>(0,imag(a));
    for(int i=max_iterations;i>=0;i--)b=(_temp*b)/(i+1)+1;
    return b*std::exp(real(a));
}

/**
 * @brief Calculates exponential with default iterations
 * @param a Complex number
 * @return Exponential value
 */
template<typename _Tp>complex<_Tp> exp(complex<_Tp> a){return exp(a,20);}

/**
 * @brief Creates complex number from polar coordinates
 * @param norm Magnitude
 * @param arg Argument
 * @return Complex number
 */
template<typename _Tp>complex<_Tp> polar(_Tp norm,_Tp arg){
    return exp(complex<_Tp>(std::log(norm),arg));
}

/**
 * @brief Calculates natural logarithm of complex number
 * @param a Complex number
 * @param max_iterations Maximum iterations for convergence
 * @return Natural logarithm
 */
template<typename _Tp>complex<_Tp> log(complex<_Tp> a,int max_iterations){
    if(!a)throw domain_error("Logarithm of zero is undefined");
    if(!imag(a)&&real(a)>0)return std::log(real(a));
    _Tp L=norm(a);
    if(L<0)throw domain_error("Logarithm of non-positive number");
    L=norm_norm(a);
    try{
        if(L>1)return log(a>>1,max_iterations)+_Tp(.69314718055994530942);
        if(L<0.25)return log(a<<1,max_iterations)-_Tp(.69314718055994530942);
        complex<_Tp> Q=1,q=1;
        for(int i=1;i<max_iterations;++i)q=q/2;
        auto AGM=[](complex<_Tp> a,complex<_Tp> b,int max_iterations)->complex<_Tp>{
            complex<_Tp> s=a+b,p=a*b;
            for(int i=0;i<max_iterations;++i){
                a=s/2;
                b=sqrt(p);
                s=a+b;
                p=a*b;
            }
            return a;
        };
        complex<_Tp> ans1=AGM(Q,q,8),ans2=AGM(Q,a*q,8);
        _Tp Lans1=norm_norm(ans1),Lans2=norm_norm(ans2);
        if(Lans1<std::numeric_limits<_Tp>::epsilon()
        ||Lans2<std::numeric_limits<_Tp>::epsilon())throw complex_convergence_error();
        return (Q/ans1-Q/ans2)*complex<_Tp>::pi1_2();
    }catch(const complex_error&){throw complex_convergence_error();}
}

/**
 * @brief Calculates natural logarithm with default iterations
 * @param a Complex number
 * @return Natural logarithm
 */
template<typename _Tp>complex<_Tp> log(complex<_Tp> a){return log(a,40);}

/**
 * @brief Calculates logarithm with custom base
 * @param a Complex number
 * @param b Base
 * @return Logarithm value
 */
template<typename _Tp>complex<_Tp> log(complex<_Tp> a,complex<_Tp> b){
    return log(a)/log(b);
}

/**
 * @brief Calculates argument of complex number
 * @param a Complex number
 * @return Argument value
 */
template<typename _Tp>_Tp arg(complex<_Tp> a){return imag(log(a));}

/**
 * @brief Calculates phase of complex number
 * @param a Complex number
 * @return Phase value
 */
template<typename _Tp>_Tp phase(complex<_Tp> a){return arg(a);}

/**
 * @brief Calculates cosine of complex number
 * @param a Complex number
 * @param max_iterations Maximum iterations for convergence
 * @return Cosine value
 */
template<typename _Tp>complex<_Tp> cos(complex<_Tp> a,int max_iterations){
    if(!imag(a))return std::cos(real(a));
    while(real(a)>complex<_Tp>::pi2())a-=complex<_Tp>::pi2();
    while(real(a)<-complex<_Tp>::pi2())a+=complex<_Tp>::pi2();
    _Tp L=norm_norm(a);
    if(L>1){
        complex<_Tp> b=cos(a>>1,max_iterations);
        return 2*b*b-1;
    }
    complex<_Tp> b=1;
    for(int i=max_iterations;i>=0;i--){
        if(i%2)b=(a*b)/(i+1);
        else if(i%4)b=(a*b)/(i+1)-1;
        else b=(a*b)/(i+1)+1;
    }
    return b;
}

/**
 * @brief Calculates cosine with default iterations
 * @param a Complex number
 * @return Cosine value
 */
template<typename _Tp>complex<_Tp> cos(complex<_Tp> a){return cos(a,20);}

/**
 * @brief Calculates hyperbolic cosine
 * @param a Complex number
 * @return Hyperbolic cosine value
 */
template<typename _Tp>complex<_Tp> cosh(complex<_Tp> a){
    return (exp(a)+exp(-a))/2;
}

/**
 * @brief Calculates hyperbolic sine
 * @param a Complex number
 * @return Hyperbolic sine value
 */
template<typename _Tp>complex<_Tp> sinh(complex<_Tp> a){
    return (exp(a)-exp(-a))/2;
}

/**
 * @brief Calculates hyperbolic tangent
 * @param a Complex number
 * @return Hyperbolic tangent value
 */
template<typename _Tp>complex<_Tp> tanh(complex<_Tp> a){
    return sinh(a)/cosh(a);
}

/**
 * @brief Calculates hyperbolic cotangent
 * @param a Complex number
 * @return Hyperbolic cotangent value
 */
template<typename _Tp>complex<_Tp> coth(complex<_Tp> a){
    return cosh(a)/sinh(a);
}

/**
 * @brief Calculates hyperbolic secant
 * @param a Complex number
 * @return Hyperbolic secant value
 */
template<typename _Tp>complex<_Tp> sech(complex<_Tp> a){
    return 1/cosh(a);
}

/**
 * @brief Calculates hyperbolic cosecant
 * @param a Complex number
 * @return Hyperbolic cosecant value
 */
template<typename _Tp>complex<_Tp> csch(complex<_Tp> a){
    return 1/sinh(a);
}

/**
 * @brief Calculates inverse hyperbolic cosine
 * @param a Complex number
 * @return Inverse hyperbolic cosine value
 */
template<typename _Tp>complex<_Tp> acosh(complex<_Tp> a){
    return log(a+sqrt(a*a-1));
}

/**
 * @brief Calculates inverse hyperbolic sine
 * @param a Complex number
 * @return Inverse hyperbolic sine value
 */
template<typename _Tp>complex<_Tp> asinh(complex<_Tp> a){
    return log(a+sqrt(a*a+1));
}

/**
 * @brief Calculates inverse hyperbolic tangent
 * @param a Complex number
 * @return Inverse hyperbolic tangent value
 */
template<typename _Tp>complex<_Tp> atanh(complex<_Tp> a){
    return log((1+a)/(1-a))/2;
}

/**
 * @brief Calculates inverse hyperbolic cotangent
 * @param a Complex number
 * @return Inverse hyperbolic cotangent value
 */
template<typename _Tp>complex<_Tp> acoth(complex<_Tp> a){
    return log((a+1)/(a-1))/2;
}

/**
 * @brief Calculates inverse hyperbolic secant
 * @param a Complex number
 * @return Inverse hyperbolic secant value
 */
template<typename _Tp>complex<_Tp> asech(complex<_Tp> a){
    return log((sqrt(1-a*a)+1)/a);
}

/**
 * @brief Calculates inverse hyperbolic cosecant
 * @param a Complex number
 * @return Inverse hyperbolic cosecant value
 */
template<typename _Tp>complex<_Tp> acsch(complex<_Tp> a){
    return log((sqrt(1+a*a)+1)/a);
}

/**
 * @brief Calculates sine
 * @param a Complex number
 * @return Sine value
 */
template<typename _Tp>complex<_Tp> sin(complex<_Tp> a){
	return cos(complex<_Tp>::pi1_2()-a);
}

/**
 * @brief Calculates tangent
 * @param a Complex number
 * @return Tangent value
 */
template<typename _Tp>complex<_Tp> tan(complex<_Tp> a){
    return sin(a)/cos(a);
}

/**
 * @brief Calculates cotangent
 * @param a Complex number
 * @return Cotangent value
 */
template<typename _Tp>complex<_Tp> cot(complex<_Tp> a){
    return cos(a)/sin(a);
}

/**
 * @brief Calculates secant
 * @param a Complex number
 * @return Secant value
 */
template<typename _Tp>complex<_Tp> sec(complex<_Tp> a){
    return 1/cos(a);
}

/**
 * @brief Calculates cosecant
 * @param a Complex number
 * @return Cosecant value
 */
template<typename _Tp>complex<_Tp> csc(complex<_Tp> a){
    return 1/sin(a);
}

/**
 * @brief Calculates hyper-exponential function
 * @param a Complex number
 * @param n Number of iterations
 * @return Hyper-exponential value
 */
template<typename _Tp>complex<_Tp> hyper_xexp(complex<_Tp> a,int n){
    complex<_Tp> ans=a;
    for(int i=1;i<=n;i++)ans=exp(ans);
    for(int i=1;i<=-n;i++)ans=log(ans);
    return ans*a;
}

/**
 * @brief Raises complex number to complex power
 * @param base Base complex number
 * @param exponent Complex exponent
 * @return Result of power operation
 */
template<typename _Tp>complex<_Tp> pow(complex<_Tp> base,complex<_Tp> exponent){
    if(!exponent)return 1;
    if(!base)return 0;
    if(real(exponent)==std::round(real(exponent))&&(!imag(exponent)))
        return pow(base,int(real(exponent)));
    return exp(log(base)*exponent);
}

/**
 * @brief Calculates complex root
 * @param base Base complex number
 * @param exponent Complex exponent
 * @return Complex root value
 */
template<typename _Tp>complex<_Tp> root(complex<_Tp> base,complex<_Tp> exponent){
    if(exponent==1)return base;
    if(!base)return 0;
    if(real(exponent)==std::round(real(exponent))&&(!imag(exponent)))
        return root(base,int(real(exponent)));
    return exp(log(base)/exponent);
}

/**
 * @brief Calculates definite integral
 * @param a Lower bound
 * @param b Upper bound
 * @param func Function to integrate
 * @param delta Step size
 * @return Integral value
 */
template<typename _Tp,typename _Func>complex<_Tp> integral(complex<_Tp> a,complex<_Tp> b,const _Func func,_Tp delta){
    complex<_Tp> d=delta*(b-a),ans,ins=func(a)*d;
    for(int i=0;delta*_Tp(i)<=1;++i){
        ans+=ins;
        a+=d;
        ins=func(a)*d;
    }
    return ans;
}

/**
 * @brief Calculates definite integral with default step size
 * @param a Lower bound
 * @param b Upper bound
 * @param func Function to integrate
 * @return Integral value
 */
template<typename _Tp,typename _Func>complex<_Tp> integral(complex<_Tp> a,complex<_Tp> b,const _Func func){
    return integral(a,b,func,1e-6);
}

/**
 * @brief Calculates derivative
 * @param a Point of differentiation
 * @param func Function to differentiate
 * @param delta Step size
 * @return Derivative value
 */
template<typename _Tp,typename _Func>complex<_Tp> differential(complex<_Tp> a,const _Func func,_Tp delta){
    complex<_Tp> d_real=(func(a+delta)-func(a))/delta;
    complex<_Tp> d_imag=(func(a+complex<_Tp>(0,delta))-func(a))/complex<_Tp>(0,delta);
    complex<_Tp> d=abs((d_real-d_imag)/func(a));
    if(real(d)>1e-3||imag(d)>1e-3)throw domain_error("The function is not differentiable");
    return (d_real+d_imag)/2;
}

/**
 * @brief Calculates derivative with default step size
 * @param a Point of differentiation
 * @param func Function to differentiate
 * @return Derivative value
 */
template<typename _Tp,typename _Func>complex<_Tp> differential(complex<_Tp> a,const _Func func){
    return differential(a,func,1e-6);
}

/**
 * @brief Finds inverse function value
 * @param a Target value
 * @param func Function to invert
 * @param max_iterations Maximum iterations
 * @return Inverse function value
 */
template<typename _Tp,typename _Func>complex<_Tp> invert(complex<_Tp> a,const _Func func,int max_iterations){
    complex<_Tp> b=1,bp=b;
    for(int i=0;i<max_iterations;++i){
        try{
            b=bp-func(bp)/differential(bp,func)+a/differential(bp,func);
            complex<_Tp> d(real(b)-real(bp),imag(b)-imag(bp));
            d=abs(d);
            if(real(d)<1e-3&&imag(d)<1e-3)break;
        }catch(const complex_error&){throw complex_convergence_error();}
        bp=b;
    }
    return b;
}

/**
 * @brief Finds inverse function value with default iterations
 * @param a Target value
 * @param func Function to invert
 * @return Inverse function value
 */
template<typename _Tp,typename _Func>complex<_Tp> invert(complex<_Tp> a,const _Func func){
    return invert(a,func,100);
}

/**
 * @brief Calculates inverse cosine
 * @param a Complex number
 * @return Inverse cosine value
 */
template<typename _Tp>complex<_Tp> acos(complex<_Tp> a){
    return invert(a,(complex<_Tp>(*)(complex<_Tp>))cos);
}

/**
 * @brief Calculates inverse sine
 * @param a Complex number
 * @return Inverse sine value
 */
template<typename _Tp>complex<_Tp> asin(complex<_Tp> a){
    return complex<_Tp>::pi1_2()-acos(a);
}

/**
 * @brief Calculates inverse tangent
 * @param a Complex number
 * @return Inverse tangent value
 */
template<typename _Tp>complex<_Tp> atan(complex<_Tp> a){
    return acos(_Tp(1.)/sqrt(pow(a,2)+_Tp(1.)));
}

/**
 * @brief Calculates inverse cotangent
 * @param a Complex number
 * @return Inverse cotangent value
 */
template<typename _Tp>complex<_Tp> acot(complex<_Tp> a){
    return complex<_Tp>::pi1_2()-atan(a);
}

/**
 * @brief Calculates inverse secant
 * @param a Complex number
 * @return Inverse secant value
 */
template<typename _Tp>complex<_Tp> asec(complex<_Tp> a){
    return acos(_Tp(1.)/a);
}

/**
 * @brief Calculates inverse cosecant
 * @param a Complex number
 * @return Inverse cosecant value
 */
template<typename _Tp>complex<_Tp> acsc(complex<_Tp> a){
    return complex<_Tp>::pi1_2()-asec(a);
}

/**
 * @brief Calculates hyper-omega function
 * @param a Complex number
 * @param n Number of iterations
 * @return Hyper-omega value
 */
template<typename _Tp>complex<_Tp> hyper_omega(complex<_Tp> a,int n){
    if(n>=0)return invert(a,[n](complex<_Tp> x)->complex<_Tp>{return hyper_xexp(x,n);});
    complex<_Tp> ans=hyper_omega(a,-n);
    for(int i=1;i<=-n;++i)ans=exp(ans);
    return ans;
}

/**
 * @brief Calculates gamma function
 * @param a Complex number
 * @param max_iterations Maximum iterations
 * @return Gamma value
 */
template<typename _Tp>complex<_Tp> gamma(complex<_Tp> a,int max_iterations){
    if(a==1)return 1;
    if(!imag(a)&&real(a)<0&&round(real(a))==real(a))
        throw domain_error("Gamma of non-positive integer is undefined");
    if(real(a)<.5)
        return complex<_Tp>::pi()/(sin(a*complex<_Tp>::pi())*gamma(1-a,max_iterations));
    if(real(a)>1)return gamma(a-1,max_iterations)*(a-1);
    if(real(a)<.5)
        return complex<_Tp>::pi()/(sin(a*complex<_Tp>::pi())*gamma(1-a,max_iterations));
    complex<_Tp> b=-log(a);
    for(int i=1;i<=max_iterations;i++){
        b-=log(1+a/i);
        b+=a*std::log(1+1./i);
    }
    return exp(b);
}

/**
 * @brief Calculates gamma function with default iterations
 * @param a Complex number
 * @return Gamma value
 */
template<typename _Tp>complex<_Tp> gamma(complex<_Tp> a){
    return gamma(a,100000);
}

/**
 * @brief Calculates beta function
 * @param a First parameter
 * @param b Second parameter
 * @return Beta value
 */
template<typename _Tp>complex<_Tp> beta(complex<_Tp> a,complex<_Tp> b){
    return gamma(a)*gamma(b)/gamma(a+b);
}

template<typename _Tp,typename Char_Tp,typename Traits>
std::basic_ostream<Char_Tp, Traits>& numbero(std::basic_ostream<Char_Tp, Traits>& os){complex<_Tp>::print_number=1;return os;}

template<typename _Tp,typename Char_Tp,typename Traits>
std::basic_ostream<Char_Tp, Traits>& coordinateo(std::basic_ostream<Char_Tp, Traits>& os){complex<_Tp>::print_number=0;return os;}

template<typename _Tp,typename Char_Tp,typename Traits>
std::basic_istream<Char_Tp, Traits>& numberi(std::basic_istream<Char_Tp, Traits>& is){complex<_Tp>::print_number=1;return is;}

template<typename _Tp,typename Char_Tp,typename Traits>
std::basic_istream<Char_Tp, Traits>& coordinatei(std::basic_istream<Char_Tp, Traits>& is){complex<_Tp>::print_number=0;return is;}

/**
 * @brief Output stream operator
 * @param os Output stream
 * @param z Complex number
 * @return Output stream
 */
template<typename _Tp,typename Char_Tp,typename Traits>
std::basic_ostream<Char_Tp, Traits>& operator<<(std::basic_ostream<Char_Tp,Traits>& os,const complex<_Tp>& z){
    std::basic_ostringstream<Char_Tp, Traits> s;
    s.flags(os.flags());
    s.imbue(os.getloc());
    s.precision(os.precision());
    if(complex<_Tp>::print_number){
		s<<real(z);
		if(imag(z)){
			if(imag(z)>0)s<<" + ";
			s<<imag(z)<<" i";
		}
    }
	else s<<"( "<<real(z)<<" , "<<imag(z)<<" )";
    return os<<s.str();
}

/**
 * @brief Input stream operator
 * @param is Input stream
 * @param z Complex number
 * @return Input stream
 */
template<typename _Tp,typename Char_Tp, typename Traits>
std::basic_istream<Char_Tp,Traits>& operator>>(std::basic_istream<Char_Tp,Traits>& is,complex<_Tp>& z){
    _Tp real,imag;
    Char_Tp ch;
    if(!(is>>ch)) return is;
    if(complex<_Tp>::print_number)
    {
		if(!(is>>real)) return is;
		if(!(is>>ch)) return is;
		if(ch!='+'&&ch!='-'){is.setstate(std::ios::failbit);return is;}
		if(!(is>>imag)) return is;
		if(!(is>>ch)) return is;
    	if(ch!='i'){is.setstate(std::ios::failbit);return is;}
	}
	else{
    	if(ch!='('){is.setstate(std::ios::failbit);return is;}
    	if(!(is>>real)) return is;
    	if(!(is>>ch)) return is;
    	if(ch!=','){is.setstate(std::ios::failbit);return is;}
    	if(!(is>>imag)) return is;
    	if(!(is>>ch)) return is;
    	if(ch!=')'){is.setstate(std::ios::failbit);return is;}
	}
    z=complex<_Tp>(real,imag);
    return is;
}

/**
 * @brief Calculates M function
 * @param a First parameter (complex number)
 * @param c Third parameter (complex number)
 * @param func Function pointer parameter, default is square function
 * @param L Flag parameter, default is 0
 * @return M function value (complex number)
 * 
 * This function template calculates different expressions based on flag L:
 * - When L is true: calculates (a - func(a)) + c
 * - When L is false: calculates func(a) + c
 */
template<typename _Tp>complex<_Tp> M(complex<_Tp> a, complex<_Tp> c, 
              const std::function<complex<_Tp>(complex<_Tp>)> func=[](complex<_Tp> x)->complex<_Tp>{return pow(x,2);},
              bool L = 0){
    return (L?(a-func(a)):func(a))+c;
}

/**
 * @brief Calculates T function
 * @param a First parameter (complex number)
 * @param c Third parameter (complex number)
 * @param func Function pointer parameter, default is square function
 * @param L Flag parameter, default is 0
 * @return T function value (complex number)
 * 
 * This function template is implemented by calling M function with:
 * - The complex conjugate of a as the first parameter of M
 * - Other parameters passed directly to M
 */
template<typename _Tp>complex<_Tp> T(complex<_Tp> a, complex<_Tp> c,
              const std::function<complex<_Tp>(complex<_Tp>)> func=[](complex<_Tp> x)->complex<_Tp>{return pow(x,2);},
              bool L = 0){
    return M(conj(a),c,func,L);
}

/**
 * @brief Calculates B function
 * @param a First parameter (complex number)
 * @param c Third parameter (complex number)
 * @param func Function pointer parameter, default is square function
 * @param L Flag parameter, default is 0
 * @return B function value (complex number)
 * 
 * This function template is implemented by calling M function with:
 * - The absolute value of a as the first parameter of M
 * - Other parameters passed directly to M
 */
template<typename _Tp>complex<_Tp> B(complex<_Tp> a, complex<_Tp> c,
              const std::function<complex<_Tp>(complex<_Tp>)> func=[](complex<_Tp> x)->complex<_Tp>{return pow(x,2);},
              bool L = 0){
    return M(abs(a),c,func,L);
}

/**
 * @brief Sorts complex numbers
 * @param first Start iterator
 * @param last End iterator
 * @param cmp Comparison function
 */
template<typename _Tp>void sort(complex<_Tp> *first,complex<_Tp> *last,
    bool (*cmp)(const complex<_Tp>,const complex<_Tp>)=
    [](complex<_Tp> a,complex<_Tp> b)->bool{
        return (norm(a)==norm(b))?(a<b):(norm(a)<norm(b));
    }){
    std::sort(first,last,cmp);
}

namespace batch{

/**
 * @brief Batch positive operation
 * @param dest Destination array
 * @param n Size of array
 */
template<typename _Tp>void pos(complex<_Tp> *dest,size_t n){
    complex<_Tp>::for_each(dest,n,[](complex<_Tp> x)->complex<_Tp>{return +x;});
}

/**
 * @brief Batch negative operation
 * @param dest Destination array
 * @param n Size of array
 */
template<typename _Tp>void neg(complex<_Tp> *dest,size_t n){
    complex<_Tp>::for_each(dest,n,[](complex<_Tp> x)->complex<_Tp>{return -x;});
}

/**
 * @brief Batch increment operation
 * @param dest Destination array
 * @param n Size of array
 */
template<typename _Tp>void inc(complex<_Tp> *dest,size_t n){
    complex<_Tp>::for_each(dest,n,[](complex<_Tp> x)->complex<_Tp>{return ++x;});
}

/**
 * @brief Batch decrement operation
 * @param dest Destination array
 * @param n Size of array
 */
template<typename _Tp>void dec(complex<_Tp> *dest,size_t n){
    complex<_Tp>::for_each(dest,n,[](complex<_Tp> x)->complex<_Tp>{return --x;});
}

/**
 * @brief Batch bitwise NOT operation
 * @param dest Destination array
 * @param n Size of array
 */
template<typename _Tp>void nott(complex<_Tp> *dest,size_t n){
    complex<_Tp>::for_each(dest,n,[](complex<_Tp> x)->complex<_Tp>{return ~x;});
}

/**
 * @brief Batch addition with scalar
 * @param dest Destination array
 * @param src Scalar value
 * @param n Size of array
 */
template<typename _Tp>void add(complex<_Tp> *dest,const _Tp src,size_t n){
    for(size_t i=0;i<n;++i)dest[i]+=src;
}

/**
 * @brief Batch subtraction with scalar
 * @param dest Destination array
 * @param src Scalar value
 * @param n Size of array
 */
template<typename _Tp>void sub(complex<_Tp> *dest,const _Tp src,size_t n){
    for(size_t i=0;i<n;++i)dest[i]-=src;
}

/**
 * @brief Batch multiplication with scalar
 * @param dest Destination array
 * @param src Scalar value
 * @param n Size of array
 */
template<typename _Tp>void mul(complex<_Tp> *dest,const _Tp src,size_t n){
    for(size_t i=0;i<n;++i)dest[i]*=src;
}

/**
 * @brief Batch division with scalar
 * @param dest Destination array
 * @param src Scalar value
 * @param n Size of array
 */
template<typename _Tp>void div(complex<_Tp> *dest,const _Tp src,size_t n){
    for(size_t i=0;i<n;++i)dest[i]/=src;
}

/**
 * @brief Batch modulo with scalar
 * @param dest Destination array
 * @param src Scalar value
 * @param n Size of array
 */
template<typename _Tp>void mod(complex<_Tp> *dest,const _Tp src,size_t n){
    for(size_t i=0;i<n;++i)dest[i]%=src;
}

/**
 * @brief Batch addition with array
 * @param dest Destination array
 * @param src Source array
 * @param n Size of array
 */
template<typename _Tp>void add(complex<_Tp> *dest,const _Tp *src,size_t n){
    for(size_t i=0;i<n;++i)dest[i]+=src[i];
}

/**
 * @brief Batch subtraction with array
 * @param dest Destination array
 * @param src Source array
 * @param n Size of array
 */
template<typename _Tp>void sub(complex<_Tp> *dest,const _Tp *src,size_t n){
    for(size_t i=0;i<n;++i)dest[i]-=src[i];
}

/**
 * @brief Batch multiplication with array
 * @param dest Destination array
 * @param src Source array
 * @param n Size of array
 */
template<typename _Tp>void mul(complex<_Tp> *dest,const _Tp *src,size_t n){
    for(size_t i=0;i<n;++i)dest[i]*=src[i];
}

/**
 * @brief Batch division with array
 * @param dest Destination array
 * @param src Source array
 * @param n Size of array
 */
template<typename _Tp>void div(complex<_Tp> *dest,const _Tp *src,size_t n){
    for(size_t i=0;i<n;++i)dest[i]/=src[i];
}

/**
 * @brief Batch modulo with array
 * @param dest Destination array
 * @param src Source array
 * @param n Size of array
 */
template<typename _Tp>void mod(complex<_Tp> *dest,const _Tp *src,size_t n){
    for(size_t i=0;i<n;++i)dest[i]%=src[i];
}

/**
 * @brief Batch addition with complex scalar
 * @param dest Destination array
 * @param src Complex scalar
 * @param n Size of array
 */
template<typename _Tp>void add(complex<_Tp> *dest,const complex<_Tp> src,size_t n){
    for(size_t i=0;i<n;++i)dest[i]+=src;
}

/**
 * @brief Batch subtraction with complex scalar
 * @param dest Destination array
 * @param src Complex scalar
 * @param n Size of array
 */
template<typename _Tp>void sub(complex<_Tp> *dest,const complex<_Tp> src,size_t n){
    for(size_t i=0;i<n;++i)dest[i]-=src;
}

/**
 * @brief Batch multiplication with complex scalar
 * @param dest Destination array
 * @param src Complex scalar
 * @param n Size of array
 */
template<typename _Tp>void mul(complex<_Tp> *dest,const complex<_Tp> src,size_t n){
    for(size_t i=0;i<n;++i)dest[i]*=src;
}

/**
 * @brief Batch division with complex scalar
 * @param dest Destination array
 * @param src Complex scalar
 * @param n Size of array
 */
template<typename _Tp>void div(complex<_Tp> *dest,const complex<_Tp> src,size_t n){
    for(size_t i=0;i<n;++i)dest[i]/=src;
}

/**
 * @brief Batch modulo with complex scalar
 * @param dest Destination array
 * @param src Complex scalar
 * @param n Size of array
 */
template<typename _Tp>void mod(complex<_Tp> *dest,const complex<_Tp> src,size_t n){
    for(size_t i=0;i<n;++i)dest[i]%=src;
}

/**
 * @brief Batch addition with complex array
 * @param dest Destination array
 * @param src Source array
 * @param n Size of array
 */
template<typename _Tp>void add(complex<_Tp> *dest,const complex<_Tp> *src,size_t n){
    for(size_t i=0;i<n;++i)dest[i]+=src[i];
}

/**
 * @brief Batch subtraction with complex array
 * @param dest Destination array
 * @param src Source array
 * @param n Size of array
 */
template<typename _Tp>void sub(complex<_Tp> *dest,const complex<_Tp> *src,size_t n){
    for(size_t i=0;i<n;++i)dest[i]-=src[i];
}

/**
 * @brief Batch multiplication with complex array
 * @param dest Destination array
 * @param src Source array
 * @param n Size of array
 */
template<typename _Tp>void mul(complex<_Tp> *dest,const complex<_Tp> *src,size_t n){
    for(size_t i=0;i<n;++i)dest[i]*=src[i];
}

/**
 * @brief Batch division with complex array
 * @param dest Destination array
 * @param src Source array
 * @param n Size of array
 */
template<typename _Tp>void div(complex<_Tp> *dest,const complex<_Tp> *src,size_t n){
    for(size_t i=0;i<n;++i)dest[i]/=src[i];
}

/**
 * @brief Batch modulo with complex array
 * @param dest Destination array
 * @param src Source array
 * @param n Size of array
 */
template<typename _Tp>void mod(complex<_Tp> *dest,const complex<_Tp> *src,size_t n){
    for(size_t i=0;i<n;++i)dest[i]%=src[i];
}

/**
 * @brief Batch left shift
 * @param dest Destination array
 * @param src Shift amount
 * @param n Size of array
 */
template<typename _Tp>void lshift(complex<_Tp> *dest,const int src,size_t n){
    for(size_t i=0;i<n;++i)dest[i]<<=src;
}

/**
 * @brief Batch right shift
 * @param dest Destination array
 * @param src Shift amount
 * @param n Size of array
 */
template<typename _Tp>void rshift(complex<_Tp> *dest,const int src,size_t n){
    for(size_t i=0;i<n;++i)dest[i]>>=src;
}

/**
 * @brief Batch bitwise AND with scalar
 * @param dest Destination array
 * @param src Scalar value
 * @param n Size of array
 */
template<typename _Tp>void andd(complex<_Tp> *dest,const _Tp src,size_t n){
    for(size_t i=0;i<n;++i)dest[i]&=src;
}

/**
 * @brief Batch bitwise OR with scalar
 * @param dest Destination array
 * @param src Scalar value
 * @param n Size of array
 */
template<typename _Tp>void orr(complex<_Tp> *dest,const _Tp src,size_t n){
    for(size_t i=0;i<n;++i)dest[i]|=src;
}

/**
 * @brief Batch bitwise XOR with scalar
 * @param dest Destination array
 * @param src Scalar value
 * @param n Size of array
 */
template<typename _Tp>void xorr(complex<_Tp> *dest,const _Tp src,size_t n){
    for(size_t i=0;i<n;++i)dest[i]^=src;
}

/**
 * @brief Batch left shift with array
 * @param dest Destination array
 * @param src Shift amounts array
 * @param n Size of array
 */
template<typename _Tp>void lshift(complex<_Tp> *dest,const int *src,size_t n){
    for(size_t i=0;i<n;++i)dest[i]<<=src[i];
}

/**
 * @brief Batch right shift with array
 * @param dest Destination array
 * @param src Shift amounts array
 * @param n Size of array
 */
template<typename _Tp>void rshift(complex<_Tp> *dest,const int *src,size_t n){
    for(size_t i=0;i<n;++i)dest[i]>>=src[i];
}

/**
 * @brief Batch bitwise AND with array
 * @param dest Destination array
 * @param src Source array
 * @param n Size of array
 */
template<typename _Tp>void andd(complex<_Tp> *dest,const _Tp *src,size_t n){
    for(size_t i=0;i<n;++i)dest[i]&=src[i];
}

/**
 * @brief Batch bitwise OR with array
 * @param dest Destination array
 * @param src Source array
 * @param n Size of array
 */
template<typename _Tp>void orr(complex<_Tp> *dest,const _Tp *src,size_t n){
    for(size_t i=0;i<n;++i)dest[i]|=src[i];
}

/**
 * @brief Batch bitwise XOR with array
 * @param dest Destination array
 * @param src Source array
 * @param n Size of array
 */
template<typename _Tp>void xorr(complex<_Tp> *dest,const _Tp *src,size_t n){
    for(size_t i=0;i<n;++i)dest[i]^=src[i];
}

/**
 * @brief Batch bitwise AND with complex scalar
 * @param dest Destination array
 * @param src Complex scalar
 * @param n Size of array
 */
template<typename _Tp>void andd(complex<_Tp> *dest,const complex<_Tp> src,size_t n){
    for(size_t i=0;i<n;++i)dest[i]&=src;
}

/**
 * @brief Batch bitwise OR with complex scalar
 * @param dest Destination array
 * @param src Complex scalar
 * @param n Size of array
 */
template<typename _Tp>void orr(complex<_Tp> *dest,const complex<_Tp> src,size_t n){
    for(size_t i=0;i<n;++i)dest[i]|=src;
}

/**
 * @brief Batch bitwise XOR with complex scalar
 * @param dest Destination array
 * @param src Complex scalar
 * @param n Size of array
 */
template<typename _Tp>void xorr(complex<_Tp> *dest,const complex<_Tp> src,size_t n){
    for(size_t i=0;i<n;++i)dest[i]^=src;
}

/**
 * @brief Batch bitwise AND with complex array
 * @param dest Destination array
 * @param src Source array
 * @param n Size of array
 */
template<typename _Tp>void andd(complex<_Tp> *dest,const complex<_Tp> *src,size_t n){
    for(size_t i=0;i<n;++i)dest[i]&=src[i];
}

/**
 * @brief Batch bitwise OR with complex array
 * @param dest Destination array
 * @param src Source array
 * @param n Size of array
 */
template<typename _Tp>void orr(complex<_Tp> *dest,const complex<_Tp> *src,size_t n){
    for(size_t i=0;i<n;++i)dest[i]|=src[i];
}

/**
 * @brief Batch bitwise XOR with complex array
 * @param dest Destination array
 * @param src Source array
 * @param n Size of array
 */
template<typename _Tp>void xorr(complex<_Tp> *dest,const complex<_Tp> *src,size_t n){
    for(size_t i=0;i<n;++i)dest[i]^=src[i];
}

/**
 * @brief Batch complex conjugate operation
 * @param data Array of complex numbers
 * @param n Size of array
 */
template<typename _Tp>void conj(complex<_Tp> *data,size_t n)noexcept{
    complex<_Tp>::for_each(data,n,zhangzl::conj);
}

/**
 * @brief Batch absolute value operation
 * @param data Array of complex numbers
 * @param n Size of array
 */
template<typename _Tp>void abs(complex<_Tp> *data,size_t n)noexcept{
    complex<_Tp>::for_each(data,n,zhangzl::abs);
}

/**
 * @brief Batch square root operation
 * @param data Array of complex numbers
 * @param n Size of array
 */
template<typename _Tp>void sqrt(complex<_Tp> *data,size_t n){
    complex<_Tp>::for_each(data,n,zhangzl::sqrt);
}

/**
 * @brief Batch cube root operation
 * @param data Array of complex numbers
 * @param n Size of array
 */
template<typename _Tp>void cbrt(complex<_Tp> *data,size_t n){
    complex<_Tp>::for_each(data,n,zhangzl::cbrt);
}

/**
 * @brief Batch nth root operation
 * @param dest Destination array
 * @param src Root degree
 * @param n Size of array
 */
template<typename _Tp>void root(complex<_Tp> *dest,const int src,size_t n){
    for(size_t i=0;i<n;++i)dest[i]=root(dest[i],src);
}

/**
 * @brief Batch nth root operation
 * @param dest Destination array
 * @param src Root degrees array
 * @param n Size of array
 */
template<typename _Tp>void root(complex<_Tp> *dest,const int *src,size_t n){
    for(size_t i=0;i<n;++i)dest[i]=root(dest[i],src[i]);
}

/**
 * @brief Batch exponential operation
 * @param data Array of complex numbers
 * @param n Size of array
 */
template<typename _Tp>void exp(complex<_Tp> *data,size_t n){
    complex<_Tp>::for_each(data,n,zhangzl::exp);
}

/**
 * @brief Batch logarithm operation
 * @param data Array of complex numbers
 * @param n Size of array
 */
template<typename _Tp>void log(complex<_Tp> *data,size_t n){
    complex<_Tp>::for_each(data,n,zhangzl::log);
}

/**
 * @brief Batch cosine operation
 * @param data Array of complex numbers
 * @param n Size of array
 */
template<typename _Tp>void cos(complex<_Tp> *data,size_t n){
    complex<_Tp>::for_each(data,n,cos);
}

/**
 * @brief Batch hyperbolic cosine operation
 * @param data Array of complex numbers
 * @param n Size of array
 */
template<typename _Tp>void cosh(complex<_Tp> *data,size_t n){
    complex<_Tp>::for_each(data,n,cosh);
}

/**
 * @brief Batch hyperbolic sine operation
 * @param data Array of complex numbers
 * @param n Size of array
 */
template<typename _Tp>void sinh(complex<_Tp> *data,size_t n){
    complex<_Tp>::for_each(data,n,sinh);
}

/**
 * @brief Batch hyperbolic tangent operation
 * @param data Array of complex numbers
 * @param n Size of array
 */
template<typename _Tp>void tanh(complex<_Tp> *data,size_t n){
    complex<_Tp>::for_each(data,n,tanh);
}

/**
 * @brief Batch hyperbolic cotangent operation
 * @param data Array of complex numbers
 * @param n Size of array
 */
template<typename _Tp>void coth(complex<_Tp> *data,size_t n){
    complex<_Tp>::for_each(data,n,coth);
}

/**
 * @brief Batch hyperbolic secant operation
 * @param data Array of complex numbers
 * @param n Size of array
 */
template<typename _Tp>void sech(complex<_Tp> *data,size_t n){
    complex<_Tp>::for_each(data,n,sech);
}

/**
 * @brief Batch hyperbolic cosecant operation
 * @param data Array of complex numbers
 * @param n Size of array
 */
template<typename _Tp>void csch(complex<_Tp> *data,size_t n){
    complex<_Tp>::for_each(data,n,csch);
}

/**
 * @brief Batch inverse hyperbolic cosine operation
 * @param data Array of complex numbers
 * @param n Size of array
 */
template<typename _Tp>void acosh(complex<_Tp> *data,size_t n){
    complex<_Tp>::for_each(data,n,acosh);
}

/**
 * @brief Batch inverse hyperbolic sine operation
 * @param data Array of complex numbers
 * @param n Size of array
 */
template<typename _Tp>void asinh(complex<_Tp> *data,size_t n){
    complex<_Tp>::for_each(data,n,asinh);
}

/**
 * @brief Batch inverse hyperbolic tangent operation
 * @param data Array of complex numbers
 * @param n Size of array
 */
template<typename _Tp>void atanh(complex<_Tp> *data,size_t n){
    complex<_Tp>::for_each(data,n,atanh);
}

/**
 * @brief Batch inverse hyperbolic cotangent operation
 * @param data Array of complex numbers
 * @param n Size of array
 */
template<typename _Tp>void acoth(complex<_Tp> *data,size_t n){
    complex<_Tp>::for_each(data,n,acoth);
}

/**
 * @brief Batch inverse hyperbolic secant operation
 * @param data Array of complex numbers
 * @param n Size of array
 */
template<typename _Tp>void asech(complex<_Tp> *data,size_t n){
    complex<_Tp>::for_each(data,n,asech);
}

/**
 * @brief Batch inverse hyperbolic cosecant operation
 * @param data Array of complex numbers
 * @param n Size of array
 */
template<typename _Tp>void acsch(complex<_Tp> *data,size_t n){
    complex<_Tp>::for_each(data,n,acsch);
}

/**
 * @brief Batch sine operation
 * @param data Array of complex numbers
 * @param n Size of array
 */
template<typename _Tp>void sin(complex<_Tp> *data,size_t n){
    complex<_Tp>::for_each(data,n,sin);
}

/**
 * @brief Batch tangent operation
 * @param data Array of complex numbers
 * @param n Size of array
 */
template<typename _Tp>void tan(complex<_Tp> *data,size_t n){
    complex<_Tp>::for_each(data,n,tan);
}

/**
 * @brief Batch cotangent operation
 * @param data Array of complex numbers
 * @param n Size of array
 */
template<typename _Tp>void cot(complex<_Tp> *data,size_t n){
    complex<_Tp>::for_each(data,n,cot);
}

/**
 * @brief Batch secant operation
 * @param data Array of complex numbers
 * @param n Size of array
 */
template<typename _Tp>void sec(complex<_Tp> *data,size_t n){
    complex<_Tp>::for_each(data,n,sec);
}

/**
 * @brief Batch cosecant operation
 * @param data Array of complex numbers
 * @param n Size of array
 */
template<typename _Tp>void csc(complex<_Tp> *data,size_t n){
    complex<_Tp>::for_each(data,n,csc);
}

/**
 * @brief Batch power operation
 * @param dest Destination array
 * @param src Complex exponent
 * @param n Size of array
 */
template<typename _Tp>void pow(complex<_Tp> *dest,const complex<_Tp> src,size_t n){
    for(size_t i=0;i<n;++i)dest[i]=pow(dest[i],src);
}

/**
 * @brief Batch power operation with array
 * @param dest Destination array
 * @param src Array of complex exponents
 * @param n Size of array
 */
template<typename _Tp>void pow(complex<_Tp> *dest,const complex<_Tp> *src,size_t n){
    for(size_t i=0;i<n;++i)dest[i]=pow(dest[i],src[i]);
}

/**
 * @brief Batch inverse cosine operation
 * @param data Array of complex numbers
 * @param n Size of array
 */
template<typename _Tp>void acos(complex<_Tp> *data,size_t n){
    complex<_Tp>::for_each(data,n,acos);
}

/**
 * @brief Batch inverse sine operation
 * @param data Array of complex numbers
 * @param n Size of array
 */
template<typename _Tp>void asin(complex<_Tp> *data,size_t n){
    complex<_Tp>::for_each(data,n,asin);
}

/**
 * @brief Batch inverse tangent operation
 * @param data Array of complex numbers
 * @param n Size of array
 */
template<typename _Tp>void atan(complex<_Tp> *data,size_t n){
    complex<_Tp>::for_each(data,n,atan);
}

/**
 * @brief Batch inverse cotangent operation
 * @param data Array of complex numbers
 * @param n Size of array
 */
template<typename _Tp>void acot(complex<_Tp> *data,size_t n){
    complex<_Tp>::for_each(data,n,acot);
}

/**
 * @brief Batch inverse secant operation
 * @param data Array of complex numbers
 * @param n Size of array
 */
template<typename _Tp>void asec(complex<_Tp> *data,size_t n){
    complex<_Tp>::for_each(data,n,asec);
}

/**
 * @brief Batch inverse cosecant operation
 * @param data Array of complex numbers
 * @param n Size of array
 */
template<typename _Tp>void acsc(complex<_Tp> *data,size_t n){
    complex<_Tp>::for_each(data,n,acsc);
}

/**
 * @brief Batch gamma function operation
 * @param data Array of complex numbers
 * @param n Size of array
 */
template<typename _Tp>void gamma(complex<_Tp> *data,size_t n){
    complex<_Tp>::for_each(data,n,gamma);
}

/**
 * @brief Batch beta function operation
 * @param dest Destination array
 * @param src Complex parameter
 * @param n Size of array
 */
template<typename _Tp>void beta(complex<_Tp> *dest,const complex<_Tp> src,size_t n){
    for(size_t i=0;i<n;++i)dest[i]=beta(dest[i],src);
}

/**
 * @brief Batch beta function operation with array
 * @param dest Destination array
 * @param src Array of complex parameters
 * @param n Size of array
 */
template<typename _Tp>void beta(complex<_Tp> *dest,const complex<_Tp> *src,size_t n){
    for(size_t i=0;i<n;++i)dest[i]=beta(dest[i],src[i]);
}

/**
 * @brief Perform transformation operations on complex number arrays
 * @tparam _Tp Data type of the complex numbers
 * @param dest Pointer to destination array for storing results
 * @param srcc Pointer to source array providing input data
 * @param n Length of the arrays
 * @param func Function pointer with default value of square function, applied to each element
 * @param L Boolean flag with default value 0. When true, modifies func to x - func(x)
 */
template<typename _Tp,typename _Func>
void M(complex<_Tp> *dest, const complex<_Tp> *srcc, size_t n, 
       _Func func=[](complex<_Tp> x)->complex<_Tp>{return pow(x,2);}, 
       const bool L = 0){
    // If L is true, modify the function to x - func(x)
    if(L)func=[func](complex<_Tp> x)->complex<_Tp>{return x-func(x);};
    // Apply func to each element in dest array
    complex<_Tp>::for_each(dest, n, func);
    // Add contents of srcc array to dest array
    add(dest,srcc,n);
}

/**
 * @brief Perform conjugate transformation on complex number arrays
 * @tparam _Tp Data type of the complex numbers
 * @param dest Pointer to destination array
 * @param srcc Pointer to source array
 * @param n Length of the arrays
 * @param func Function pointer with default square function
 * @param L Boolean flag with default value 0
 */
template<typename _Tp,typename _Func>
void T(complex<_Tp> *dest, const complex<_Tp> *srcc, size_t n,
       _Func func=[](complex<_Tp> x)->complex<_Tp>{return pow(x,2);},
       const bool L = 0){
    // Take conjugate of each element in dest array
    conj(dest, n);
    // Call M function for further processing
    M(dest, srcc, n, func, L);
}

/**
 * @brief Perform magnitude transformation on complex number arrays
 * @tparam _Tp Data type of the complex numbers
 * @param dest Pointer to destination array
 * @param srcc Pointer to source array
 * @param n Length of the arrays
 * @param func Function pointer with default square function
 * @param L Boolean flag with default value 0
 */
template<typename _Tp,typename _Func>
void B(complex<_Tp> *dest, const complex<_Tp> *srcc, size_t n,
       _Func func=[](complex<_Tp> x)->complex<_Tp>{return pow(x,2);},
       const bool L = 0){
    // Calculate magnitude of each element in dest array
    abs(dest, n);
    // Call M function for further processing
    M(dest, srcc, n, func, L);
}
}
}
}
}
}

#endif
