#ifndef FUB_POOL_ALLOCATOR_HPP
#define FUB_POOL_ALLOCATOR_HPP

#include <memory_resource>
#include <vector>
#include <unordered_map>

namespace fub {

struct monotonic_bucket_resource : std::pmr::memory_resource {
  using pointer = void*;
  using size_type = std::size_t;

  monotonic_bucket_resource(std::pmr::memory_resource& resource);
  
  monotonic_bucket_resource(const monotonic_bucket_resource&) = delete;
  monotonic_bucket_resource&
  operator=(const monotonic_bucket_resource&) = delete;

  ~monotonic_bucket_resource() noexcept;

private:
  pointer do_allocate(size_type bytes, size_type alignment) override;

  void do_deallocate(pointer p, size_type bytes, size_type alignment) override;

  bool
  do_is_equal(const std::pmr::memory_resource& other) const noexcept override;

private:
  std::pmr::memory_resource* memory_resource_;
  std::unordered_map<size_type, pointer> free_buckets_{};
  std::vector<pointer> all_pointers_{};
};

}

#endif