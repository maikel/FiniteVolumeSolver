#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/ParallelBuffer.h"
#include "SAMRAI/tbox/SAMRAIManager.h"

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

/// This Scope Guard calls all neccessary initialize and finalize functions for
/// SAMRAI.
///
/// If a test fails an exception will be thrown and this will automatically
/// ensure that the shutdown and finalize functions is properly called.
struct ScopeGuard {
  ScopeGuard() {
    SAMRAI::tbox::SAMRAIManager::initialize();
    SAMRAI::tbox::SAMRAIManager::startup();
  }

  /// Calls SAMRAI routines to do the clean-up.
  ~ScopeGuard() noexcept {
    SAMRAI::tbox::SAMRAIManager::shutdown();
    SAMRAI::tbox::SAMRAIManager::finalize();
  }
};

TEST_CASE("std::ostream::write works with ParallelBuffer") {
  using SAMRAI::tbox::pout;
  ScopeGuard samrai_scope;

  std::streambuf* buf = pout.rdbuf();

  // Cast underlying streambuf to PrallelBuffer and fail the test if this is
  // actually not valid.
  SAMRAI::tbox::ParallelBuffer* parallel_buffer =
      dynamic_cast<SAMRAI::tbox::ParallelBuffer*>(buf);
  REQUIRE(parallel_buffer != nullptr);

  // Set the new output stream of tbox::pout to a ostringstream.
  std::ostringstream output;
  parallel_buffer->setOutputStream1(&output);

  /// Use write() to write a string to the stream.
  const char* string = "My Test String\n";
  std::streamsize size = std::strlen(string);
  parallel_buffer->sputn(string, size);
  REQUIRE(std::string(string) == output.str());
}