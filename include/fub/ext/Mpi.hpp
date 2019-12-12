#include <mpi.h>

#include <string>

namespace fub {

std::string ReadAndBroadcastFile(std::string filepath, MPI_Comm comm);

}