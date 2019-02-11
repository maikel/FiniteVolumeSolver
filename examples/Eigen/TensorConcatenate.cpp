#include <unsupported/Eigen/CXX11/Tensor>
#include <iostream>

int main()
{
  Eigen::Tensor<double, 2> a(2,1);
  Eigen::Tensor<double, 2> b(2,1);
  Eigen::Tensor<double, 2> c(2,1);
  a.setConstant(1);
  b.setConstant(2);
  c.setConstant(3);
  Eigen::Tensor<double, 2> d = a.concatenate(b, 1).eval().concatenate(c, 1);
  std::cout << d.dimension(0) << '\n';
  std::cout << d.dimension(1) << '\n';
  std::cout << d << '\n';

  Eigen::Tensor<double, 2> e = a.broadcast(std::array<int, 2>{2, 2});
  std::cout << e.dimension(0) << '\n';
  std::cout << e.dimension(1) << '\n';
  std::cout << e << '\n';
}
