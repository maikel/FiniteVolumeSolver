#ifndef FUB_POOL_ALLOCATOR_HPP
#define FUB_POOL_ALLOCATOR_HPP

#if __has_include(<memory_resource>)
#include <memory_resource>
#define FUB_WITH_MEMORY_RESOURCE 1
#elif __has_include(<experimental/memory_resource>)
#include <experimental/memory_resource>
#define FUB_WITH_MEMORY_RESOURCE 1
#define FUB_WITH_EXPERIMENTAL_MEMORY_RESOURCE 1
#endif
#include <unordered_map>
#include <vector>

namespace fub {

#ifdef FUB_WITH_EXPERIMENTAL_MEMORY_RESOURCE
using std::experimental::pmr::memory_resource;
#elif FUB_WITH_MEMORY_RESOURCE
using std::pmr::memory_resource;
#endif

#ifdef FUB_WITH_MEMORY_RESOURCE

struct monotonic_bucket_resource : memory_resource {
  using pointer = void*;
  using size_type = std::size_t;

  monotonic_bucket_resource(memory_resource& resource);

  monotonic_bucket_resource(const monotonic_bucket_resource&) = delete;
  monotonic_bucket_resource&
  operator=(const monotonic_bucket_resource&) = delete;

  ~monotonic_bucket_resource() noexcept;

private:
  pointer do_allocate(size_type bytes, size_type alignment) override;

  void do_deallocate(pointer p, size_type bytes, size_type alignment) override;

  bool do_is_equal(const memory_resource& other) const noexcept override;

private:
  memory_resource* memory_resource_;
  std::unordered_map<size_type, pointer> free_buckets_{};
  std::vector<pointer> all_pointers_{};
};

#endif

} // namespace fub

#endif