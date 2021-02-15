#ifndef FUB_AMREX_MEMORY_RESOURCE_ARENA_HPP
#define FUB_AMREX_MEMORY_RESOURCE_ARENA_HPP

#include <AMReX_Arena.H>
#include <memory_resource>

namespace fub::amrex {

struct MemoryResourceArena : ::amrex::Arena {
  explicit MemoryResourceArena(std::pmr::memory_resource* resource);

  void* alloc(std::size_t bytes) override { memory_resource_->allocate(bytes); }
  void free(void* pointer) override { memory_resource_->deallocate(pointer); };

private:
  std::pmr::memory_resource* memory_resource_;
};

} // namespace fub::amrex

#endif