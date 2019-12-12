//
// Created by Patrick Denzler on 20.11.19.
//

#ifndef FUB_COUNTER_COUNTER_RESULT_HPP
#define FUB_COUNTER_COUNTER_RESULT_HPP

#include <string>

namespace fub {

struct CounterResult {
  std::string name;
  unsigned long count;
  long long min;
  long long max;
  double mean;
  double stddev;
  double variance;
};

} // namespace fub

#endif // FUB_COUNTER_COUNTER_RESULT_HPP
