// Copyright (c) 2019 Maikel Nadolski
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#ifndef FUB_META_HPP
#define FUB_META_HPP

#include "fub/core/type_traits.hpp"

namespace fub::meta {

/// @{
/// \brief A template typedef to detect the member function
///
/// This is a template typedef which can be used in conjunction with is_detected
/// to detect in a generic context if a given expression is valid or not.
template <typename Context, typename... Args>
using PreAdvanceHierarchy = decltype(
    std::declval<Context>().PreAdvanceHierarchy(std::declval<Args>()...));

template <typename Context, typename... Args>
using PostAdvanceHierarchy = decltype(
    std::declval<Context>().PostAdvanceHierarchy(std::declval<Args>()...));

template <typename Context, typename... Args>
using PreAdvanceLevel =
    decltype(std::declval<Context>().PreAdvanceLevel(std::declval<Args>()...));

template <typename Context, typename... Args>
using PostAdvanceLevel =
    decltype(std::declval<Context>().PostAdvanceLevel(std::declval<Args>()...));

template <typename T, typename... Args>
using ResetHierarchyConfiguration = decltype(
    std::declval<T>().ResetHierarchyConfiguration(std::declval<Args>()...));

template <typename T>
using GriddingAlgorithm =
    std::decay_t<decltype(*std::declval<T>().GetGriddingAlgorithm())>;

template <typename T>
using Equation = std::decay_t<decltype(std::declval<T>().GetEquation())>;
/// @}

} // namespace fub::meta

namespace fub {

/// \brief Invokes member function obj.ResetHierarchyConfiguration(grid)
///
/// This is a helper function which invokes ResetHierarchyConfiguration if the
/// specified object obj has such a member function. If obj does not have such
/// a member function the function body will be empty.
///
/// This functionality is used by generic algorithms in include/fub/solver/*
template <typename T, typename Grid>
void ResetHierarchyConfigurationIfDetected(T&& obj, Grid&& grid) {
  if constexpr (is_detected<meta::ResetHierarchyConfiguration, T, Grid>()) {
    std::forward<T>(obj).ResetHierarchyConfiguration(std::forward<Grid>(grid));
  }
}

} // namespace fub

#endif