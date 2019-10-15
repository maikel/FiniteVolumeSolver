//
// Created by Patrick Denzler on 2019-07-30.
//

#include "fub/Solver.hpp"
#include "fub/SAMRAI.hpp"

int main(int argc, char** argv)
{
  fub::samrai::ScopeGuard guard(argc, argv);

  //fub::Advection2d advection2D{{}};
  fub::IdealGasMix<1> idealGasMix(fub::Burke2012{});
  fub::samrai::DataDescription desc = fub::samrai::RegisterVariables(idealGasMix);

  SAMRAI::hier::VariableDatabase* vardb = SAMRAI::hier::VariableDatabase::getDatabase();
  vardb->printClassData(std::cout, false);
}