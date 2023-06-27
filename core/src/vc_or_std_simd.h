#pragma once
// This header files either includes <Vc/Vc> or it includes <experimental/simd> and defines a compatibility namespace Vc.
// In either way, the including code can access Vc functions.

#if __cplusplus >= 201703L
#define HAVE_STDSIMD
#endif

#ifdef HAVE_STDSIMD
// if there is c++17, we use std::experimental::simd

// define compatibility layer such that std::simd appears like Vc
#include <experimental/simd>
namespace Vc
{
  // define types
  using double_v = std::experimental::native_simd<double>;
  using int_v = std::experimental::native_simd<int>;
  
  // import math functions
  using std::experimental::abs;
  using std::experimental::exp;
  using std::experimental::log;
  using std::experimental::where;
  using std::experimental::any_of;
  using std::experimental::all_of;
  
  // define constants
  constexpr int One = 1;
  constexpr int Zero = 0;
  
  // define inline if functions
  template<typename T>
  constexpr T iif(const typename T::mask_type& mask, const T& trueValue, const T& falseValue)
  {
    T result(falseValue);
    where(mask, result) = trueValue;
    return result;
  }
  
  template<typename T> 
  constexpr T iif (bool condition, const T &trueValue, const T &falseValue)
  {
    return condition ? trueValue : falseValue;
  }

  template<typename T>
  typename T::mask_type isnegative(const T &x)
  {
    return x < 0;
  }
  
  using std::experimental::isfinite;
  
  // define reduction functions
  template<typename Mask>
  int count(const Mask &mask)
  {
    return popcount(mask);
  }
  
  template<typename T>
  typename T::value_type min(const T &x)
  {
    return hmin(x);
  }
  
  template<typename T>
  typename T::value_type max(const T &x)
  {
    return hmax(x);
  }
}
#else
// if there is no c++17, we are use Vc

#include <Vc/Vc>

namespace Vc_1
{
  // define reduction functions
  template<typename Mask>
  int count(const Mask &mask)
  {
    return mask.count();
  }
  
  template<typename T>
  typename T::value_type min(const T &x)
  {
    return x.min();
  }
  
  template<typename T>
  typename T::value_type max(const T &x)
  {
    return x.max();
  }
}
#endif
