#include "quadrature/gauss.h"

#include <array>
#include <cmath>

namespace Quadrature
{

// 1 Gauss point
template<>
std::array<double, 1> Gauss<1>::
samplingPoints()
{
  return std::array<double, 1>{0.5};
}

// 2 Gauss points
template<>
std::array<double, 2> Gauss<2>::
samplingPoints()
{
  return std::array<double, 2>{
    (-1./sqrt(3.)+1)/2.,
    (+1./sqrt(3.)+1)/2.,
  };
}

// 3 Gauss points
template<>
std::array<double, 3> Gauss<3>::
samplingPoints()
{
  return std::array<double, 3>{
    (-sqrt(3./5)+1)/2.,
    1./2,
    (+sqrt(3./5)+1)/2.,
  };
}

// 4 Gauss points
template<>
std::array<double, 4> Gauss<4>::
samplingPoints()
{
  return std::array<double, 4>{
    (-sqrt(3./7+2./7*sqrt(6./5))+1)/2.,
    (-sqrt(3./7-2./7*sqrt(6./5))+1)/2.,
    (+sqrt(3./7-2./7*sqrt(6./5))+1)/2.,
    (+sqrt(3./7+2./7*sqrt(6./5))+1)/2.
  };
}

// 5 Gauss points
template<>
std::array<double, 5> Gauss<5>::
samplingPoints()
{
  return std::array<double, 5>{
    (-1./3*sqrt(5+2.*sqrt(10./7))+1)/2.,
    (-1./3*sqrt(5-2.*sqrt(10./7))+1)/2.,
    1./2.,
    (+1./3*sqrt(5-2.*sqrt(10./7))+1)/2.,
    (+1./3*sqrt(5+2.*sqrt(10./7))+1)/2.
  };
}

// 7 Gauss points
template<>
std::array<double,6> Gauss<6>::
samplingPoints()
{
  return std::array<double,6>{
    (0.6612093864662645+1)/2.,
    (-0.6612093864662645+1)/2.,
    (-0.2386191860831969+1)/2.,
    (0.2386191860831969+1)/2.,
    (-0.9324695142031521+1)/2.,
    (0.9324695142031521+1)/2.,
  };
}

// 7 Gauss points
template<>
std::array<double,7> Gauss<7>::
samplingPoints()
{
  return std::array<double,7>{
    (0.0000000000000000+1)/2.,
    (0.4058451513773972+1)/2.,
    (-0.4058451513773972+1)/2.,
    (-0.7415311855993945+1)/2.,
    (0.7415311855993945+1)/2.,
    (-0.9491079123427585+1)/2.,
    (0.9491079123427585+1)/2.
  };
}

// 8 Gauss points
template<>
std::array<double,8> Gauss<8>::
samplingPoints()
{
  return std::array<double,8>{
    (-0.1834346424956498+1.)/2.,
    (0.1834346424956498+1.)/2.,
    (-0.5255324099163290+1.)/2.,
    (0.5255324099163290+1.)/2.,
    (-0.7966664774136267+1.)/2.,
    (0.7966664774136267+1.)/2.,
    (-0.9602898564975363+1.)/2.,
    (0.9602898564975363+1.)/2.,
  };
}

// 10 Gauss points
template<>
std::array<double,10> Gauss<10>::
samplingPoints()
{
  return std::array<double,10>{
    (-0.1488743389816312+1.)/2.,
    (0.1488743389816312+1.)/2.,
    (-0.4333953941292472+1.)/2.,
    (0.4333953941292472+1.)/2.,
    (-0.6794095682990244+1.)/2.,
    (0.6794095682990244+1.)/2.,
    (-0.8650633666889845+1.)/2.,
    (0.8650633666889845+1.)/2.,
    (-0.9739065285171717+1.)/2.,
    (0.9739065285171717+1.)/2.,
  };
}

// 12 Gauss points
template<>
std::array<double,12> Gauss<12>::
samplingPoints()
{
  return std::array<double,12>{
    (-0.1252334085114689+1.)/2.,
    (0.1252334085114689+1.)/2.,
    (-0.3678314989981802+1.)/2.,
    (0.3678314989981802+1.)/2.,
    (-0.5873179542866175+1.)/2.,
    (0.5873179542866175+1.)/2.,
    (-0.7699026741943047+1.)/2.,
    (0.7699026741943047+1.)/2.,
    (-0.9041172563704749+1.)/2.,
    (0.9041172563704749+1.)/2.,
    (-0.9815606342467192+1.)/2.,
    (0.9815606342467192+1.)/2.,
  };
}

// 16 Gauss points
template<>
std::array<double,16> Gauss<16>::
samplingPoints()
{
  return std::array<double,16>{
    (-0.0950125098376374+1.)/2.,
    (0.0950125098376374+1.)/2.,
    (-0.2816035507792589+1.)/2.,
    (0.2816035507792589+1.)/2.,
    (-0.4580167776572274+1.)/2.,
    (0.4580167776572274+1.)/2.,
    (-0.6178762444026438+1.)/2.,
    (0.6178762444026438+1.)/2.,
    (-0.7554044083550030+1.)/2.,
    (0.7554044083550030+1.)/2.,
    (-0.8656312023878318+1.)/2.,
    (0.8656312023878318+1.)/2.,
    (-0.9445750230732326+1.)/2.,
    (0.9445750230732326+1.)/2.,
    (-0.9894009349916499+1.)/2.,
    (0.9894009349916499+1.)/2.,
  };
}

// 20 Gauss points
template<>
std::array<double,20> Gauss<20>::
samplingPoints()
{
  return std::array<double,20>{
    (-0.0765265211334973+1.)/2.,
    (0.0765265211334973+1.)/2.,
    (-0.2277858511416451+1.)/2.,
    (0.2277858511416451+1.)/2.,
    (-0.3737060887154195+1.)/2.,
    (0.3737060887154195+1.)/2.,
    (-0.5108670019508271+1.)/2.,
    (0.5108670019508271+1.)/2.,
    (-0.6360536807265150+1.)/2.,
    (0.6360536807265150+1.)/2.,
    (-0.7463319064601508+1.)/2.,
    (0.7463319064601508+1.)/2.,
    (-0.8391169718222188+1.)/2.,
    (0.8391169718222188+1.)/2.,
    (-0.9122344282513259+1.)/2.,
    (0.9122344282513259+1.)/2.,
    (-0.9639719272779138+1.)/2.,
    (0.9639719272779138+1.)/2.,
    (-0.9931285991850949+1.)/2.,
    (0.9931285991850949+1.)/2.,
  };
}

// 24 Gauss points
template<>
std::array<double,24> Gauss<24>::
samplingPoints()
{
  return std::array<double,24>{
    (0.1279381953467522+1.)/2.,
    (0.0640568928626056+1.)/2.,
    (-0.1911188674736163+1.)/2.,
    (0.1911188674736163+1.)/2.,
    (-0.3150426796961634+1.)/2.,
    (0.3150426796961634+1.)/2.,
    (-0.4337935076260451+1.)/2.,
    (0.4337935076260451+1.)/2.,
    (-0.5454214713888396+1.)/2.,
    (0.5454214713888396+1.)/2.,
    (-0.6480936519369755+1.)/2.,
    (0.6480936519369755+1.)/2.,
    (-0.7401241915785544+1.)/2.,
    (0.7401241915785544+1.)/2.,
    (-0.8200019859739029+1.)/2.,
    (0.8200019859739029+1.)/2.,
    (-0.8864155270044011+1.)/2.,
    (0.8864155270044011+1.)/2.,
    (-0.9382745520027328+1.)/2.,
    (0.9382745520027328+1.)/2.,
    (-0.9747285559713095+1.)/2.,
    (0.9747285559713095+1.)/2.,
    (-0.9951872199970213+1.)/2.,
    (0.9951872199970213+1.)/2.,

  };
}

// 64 Gauss points
template<>
std::array<double,64> Gauss<64>::
samplingPoints()
{
  // source: https://pomax.github.io/bezierinfo/legendre-gauss.html
  return std::array<double,64>{
    (-0.0243502926634244+1)/2.,
    (0.0243502926634244+1)/2.,
    (-0.0729931217877990+1)/2.,
    (0.0729931217877990+1)/2.,
    (-0.1214628192961206+1)/2.,
    (0.1214628192961206+1)/2.,
    (-0.1696444204239928+1)/2.,
    (0.1696444204239928+1)/2.,
    (-0.2174236437400071+1)/2.,
    (0.2174236437400071+1)/2.,
    (-0.2646871622087674+1)/2.,
    (0.2646871622087674+1)/2.,
    (-0.3113228719902110+1)/2.,
    (0.3113228719902110+1)/2.,
    (-0.3572201583376681+1)/2.,
    (0.3572201583376681+1)/2.,
    (-0.4022701579639916+1)/2.,
    (0.4022701579639916+1)/2.,
    (-0.4463660172534641+1)/2.,
    (0.4463660172534641+1)/2.,
    (-0.4894031457070530+1)/2.,
    (0.4894031457070530+1)/2.,
    (-0.5312794640198946+1)/2.,
    (0.5312794640198946+1)/2.,
    (-0.5718956462026340+1)/2.,
    (0.5718956462026340+1)/2.,
    (-0.6111553551723933+1)/2.,
    (0.6111553551723933+1)/2.,
    (-0.6489654712546573+1)/2.,
    (0.6489654712546573+1)/2.,
    (-0.6852363130542333+1)/2.,
    (0.6852363130542333+1)/2.,
    (-0.7198818501716109+1)/2.,
    (0.7198818501716109+1)/2.,
    (-0.7528199072605319+1)/2.,
    (0.7528199072605319+1)/2.,
    (-0.7839723589433414+1)/2.,
    (0.7839723589433414+1)/2.,
    (-0.8132653151227975+1)/2.,
    (0.8132653151227975+1)/2.,
    (-0.8406292962525803+1)/2.,
    (0.8406292962525803+1)/2.,
    (-0.8659993981540928+1)/2.,
    (0.8659993981540928+1)/2.,
    (-0.8893154459951141+1)/2.,
    (0.8893154459951141+1)/2.,
    (-0.9105221370785028+1)/2.,
    (0.9105221370785028+1)/2.,
    (-0.9295691721319396+1)/2.,
    (0.9295691721319396+1)/2.,
    (-0.9464113748584028+1)/2.,
    (0.9464113748584028+1)/2.,
    (-0.9610087996520538+1)/2.,
    (0.9610087996520538+1)/2.,
    (-0.9733268277899110+1)/2.,
    (0.9733268277899110+1)/2.,
    (-0.9833362538846260+1)/2.,
    (0.9833362538846260+1)/2.,
    (-0.9910133714767443+1)/2.,
    (0.9910133714767443+1)/2.,
    (-0.9963401167719553+1)/2.,
    (0.9963401167719553+1)/2.,
    (-0.9993050417357722+1)/2.,
    (0.9993050417357722+1)/2.
  };
}

}  // namespace
