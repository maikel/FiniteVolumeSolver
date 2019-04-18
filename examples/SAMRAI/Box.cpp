#include <SAMRAI/hier/Box.h>
#include <SAMRAI/tbox/SAMRAIManager.h>
#include <SAMRAI/tbox/SAMRAI_MPI.h>

int main(int argc, char** argv) {
  SAMRAI::tbox::SAMRAI_MPI::init(&argc, &argv);
  SAMRAI::tbox::SAMRAIManager::initialize();
  SAMRAI::tbox::SAMRAIManager::startup();
  { example(); }
  SAMRAIManager::shutdown();
  SAMRAIManager::finalize();
  SAMRAI_MPI::finalize();
}

void example() {
  SAMRAI::hier::BlockId block_id(0);

  // General Box, not part of a BoxLevel.
  SAMRAI::hier::Box box{{0, 0}, {4, 3}, block_id};

  // Make Box + information of ownership, part of a BoxLevel
  SAMRAI::hier::LocalId local_id(0);
  using SAMRAI::hier::LocalId;
  const int owner_rank = 0;
  SAMRAI::hier::Box coarse_box{
      {0, 0}, {3, 2}, block_id, LocalId(0), owner_rank};
  SAMRAI::hier::Box fine_box1{{2, 2}, {5, 3}, block_id, LocalId(0), owner_rank};
  SAMRAI::hier::Box fine_box2{{6, 2}, {7, 5}, block_id, LocalId(1), owner_rank};
}