// Copyright (c) 2019 Maikel Nadolski
// Copyright (c) 2019 Patrick Denzler
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

#ifndef FUB_SAMRAI_TAGGING_CONSTANT_BOX_HPP
#define FUB_SAMRAI_TAGGING_CONSTANT_BOX_HPP

#include <SAMRAI/hier/Box.h>
#include <fub/SAMRAI/GriddingAlgorithm.hpp>

namespace fub::samrai {

class ConstantBox {
public:
  ConstantBox(const SAMRAI::hier::Box&);

  void TagCellsForRefinement(GriddingAlgorithm& gridding, int level, int tag_id,
                             Duration time_point);

private:
  SAMRAI::hier::Box box_;
};

} // namespace fub::samrai
#endif // FUB_SAMRAI_TAGGING_CONSTANT_BOX_HPP
