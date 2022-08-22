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

#ifndef FUB_OUTPUT_FACTORY_HPP
#define FUB_OUTPUT_FACTORY_HPP

#include "fub/ext/ProgramOptions.hpp"
#include "fub/output/BasicOutput.hpp"

#include <functional>
#include <map>
#include <memory>

namespace fub {

template <typename Grid> class OutputFactory {
public:
  using ProgramOptions = std::map<std::string, pybind11::object>;
  OutputFactory() = default;

  template <typename Output, typename... Args,
            typename = std::enable_if_t<std::is_constructible<
                Output, const ProgramOptions&, remove_cvref_t<Args>...>::value>,
            typename = std::enable_if_t<
                std::is_base_of<BasicOutput<Grid>, Output>::value>>
  bool RegisterOutput(std::string name, Args&&... args) {
    return factories_
        .emplace(name,
                 [=](const ProgramOptions& opts) {
                   return std::make_unique<Output>(opts, args...);
                 })
        .second;
  }

  bool Contains(const std::string& name) { return factories_.count(name) > 0; }

  std::unique_ptr<BasicOutput<Grid>> MakeOutput(const std::string& name,
                                                const ProgramOptions& opts) {
    auto it = factories_.find(name);
    if (it == factories_.end()) {
      return nullptr;
    }
    return it->second(opts);
  }

private:
  std::map<std::string, std::function<std::unique_ptr<BasicOutput<Grid>>(
                            const ProgramOptions&)>>
      factories_;
};

} // namespace fub

#endif
