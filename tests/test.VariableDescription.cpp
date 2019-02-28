#include "fub/State.hpp"
#include "fub/VariableDescription.hpp"

using namespace fub;

namespace xx {
template <typename T, typename S> struct Cons {
  BOOST_HANA_DEFINE_STRUCT(Cons, (T, x), (S, y));
};
}
using ConsShape = xx::Cons<Scalar, Vector2d>;


struct Equation : VariableDescription<ConsShape> {};

int main() {
  using boost::hana::decltype_;
  using boost::hana::make_tuple;
  using boost::hana::to_tuple;
  constexpr Equation eq{};

  static_assert(eq.Shape(cons) == make_tuple(Scalar{}, Vector2d{}));
  static_assert(eq.ValueTypes(cons) ==
                make_tuple(type_c<double>, type_c<double>));
  static_assert(decltype_(eq.Template(cons)) == type_c<template_t<xx::Cons>>);

  // static_assert(StateBase<xx::Cons<double, double>>::Names() ==
  //               make_tuple(BOOST_HANA_STRING("x"), BOOST_HANA_STRING("y")));

  static_assert(ScalarComponentType(type_c<double>, Scalar{}) ==
                type_c<double>);
  static_assert(ScalarComponentType(type_c<double>, Vector2d{}) ==
                type_c<Array<double, 2, 1>>);

  Cons<Equation> x{};
  Cons<Equation> y{};
  Cons<Equation> z = x + y;

  static_assert(Cons<Equation>::ValueTypes() == make_tuple(type_c<double>, type_c<Array<double, 2, 1>>));

  static_assert(RefComponentType(type_c<double>) == type_c<double&>);
  static_assert(RefComponentType(type_c<Array2d>) == type_c<Array2d::MapType>);
  static_assert(ConstRefComponentType(type_c<Array2d>) == type_c<Array2d::ConstMapType>);

  static_assert(RefBaseType(type_c<Cons<Equation>>) == type_c<xx::Cons<double&, Array<double, 2, 1>::MapType>>);
  static_assert(ConstRefBaseType(type_c<Cons<Equation>>) == type_c<xx::Cons<double&, Array<double, 2, 1>::MapType>>);
}