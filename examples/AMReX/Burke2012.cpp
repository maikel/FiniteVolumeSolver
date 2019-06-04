#include "fub/AMReX.hpp"
#include "fub/Solver.hpp"

#include <boost/program_options.hpp>

namespace po = boost::program_options;

po::variable_map ParseOptions(int argc, char** argv);

fub::amrex::NumericalMethod NumericalMethod(const po::variable_map& opts);

fub::amrex::GriddingAlgorithm GriddingAlgorithm(const po::variable_map& opts);

fub::RunOptions RunOptions(po::variable_map& opts);

fub::amrex::HyperbolicSplitIntegratorContext Context(po::variable_map& opts);

int Main(const po::variable_map& opts) {
  fub::HyperbolicSplitSystemSolver solver(Context(opts));
  fub::RunSimulation(solver, RunOptions(opts));
}

int main(int argc, char** argv) {
  fub::amrex::ScopeGuard _(argc, argv);
  Main(ParseOptions(argc, argv));
}