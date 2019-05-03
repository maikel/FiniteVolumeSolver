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

#ifndef FUB_TIME_STEP_ERROR_HPP
#define FUB_TIME_STEP_ERROR_HPP

#include "fub/Duration.hpp"

#include <string>
#include <system_error>

namespace fub {

enum class TimeStepErrc { success = 0, time_step_too_large = 1 };

}

namespace std {
// Tell the C++ 11 STL metaprogramming that enum TimeStepErrc
// is registered with the standard error code system
template <> struct is_error_code_enum<::fub::TimeStepErrc> : true_type {};
} // namespace std

namespace fub {

// Define a custom error code category derived from std::error_category
class TimeStepErrcCategory : public std::error_category {
public:
  // Return a short descriptive name for the category
  virtual const char* name() const noexcept override final {
    return "TimeStepError";
  }

  // Return what each enum means in text
  virtual std::string message(int c) const override final {
    switch (static_cast<TimeStepErrc>(c)) {
    case TimeStepErrc::success:
      return "time step successful";
    case TimeStepErrc::time_step_too_large:
      return "time step size too large";
    default:
      return "unknown";
    }
  }
};

inline const TimeStepErrcCategory& TimeStepErrc_category() {
  static TimeStepErrcCategory c{};
  return c;
}

inline std::error_code make_error_code(TimeStepErrc error) {
  return {static_cast<int>(error), TimeStepErrc_category()};
}

struct TimeStepTooLarge {
  Duration coarse_dt;
};

inline std::error_code make_error_code(const TimeStepTooLarge&) {
  return TimeStepErrc::time_step_too_large;
}

} // namespace fub

#endif
