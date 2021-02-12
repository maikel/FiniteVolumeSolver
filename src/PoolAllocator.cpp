#include "fub/PoolAllocator.hpp"

namespace fub {

struct bucket_node {
  void* next;
};

struct bucket_pointer {
  std::size_t bytes;
};

monotonic_bucket_resource::pointer
monotonic_bucket_resource::do_allocate(size_type bytes,
                                       size_type /* alignment */) {
  pointer& head = free_buckets_[bytes];
  // If we have a pointer hand it out
  if (head != nullptr) {
    pointer p = head;
    bucket_node* node = static_cast<bucket_node*>(head);
    pointer next = node->next;
    node->~bucket_node();
    head = next;
    return p;
  }
  // our bucket is empty... we need to allocate new pointers!
  pointer p = memory_resource_->allocate(sizeof(std::size_t) + bytes,
                                         alignof(std::size_t));
  new (p) std::size_t{bytes};
  p = static_cast<pointer>(static_cast<std::size_t*>(p) + 1);
  new (p) bucket_node{head};
  head = p;
  all_pointers_.push_back(p);
  return p;
}

void monotonic_bucket_resource::do_deallocate(pointer p, size_type /* bytes */,
                                              size_type /* salignment */) {
  std::size_t bytes = *(static_cast<std::size_t*>(p) - 1);
  pointer& head = free_buckets_[bytes];
  new (p) bucket_node{head};
  head = p;
}

monotonic_bucket_resource::~monotonic_bucket_resource() noexcept {
  for (pointer p : all_pointers_) {
    pointer original_pointer =
        static_cast<pointer>(static_cast<std::size_t*>(p) - 1);
    memory_resource_->deallocate(original_pointer, 1);
  }
}

} // namespace fub