#include <span>
#include <optional>

using namespace std;

struct BufferSetCell {
    span<byte> data;
    optional<uint8_t> valid;
}

class BufferSetIteratorBase {
    BufferCell next() {

    }
}

class VarlenBufferSetIterator : BufferSetIteratorBase {

}

class FixedlenBufferSetIterator : BufferSetIteratorBase {

}

class UntypedBufferSetIterator {
    using iterator_category = std::forward_iterator_tag;
    using difference_type   = int;
    using value_type        = BufferSetCell;
    using pointer           = BufferSetCell;
    using reference         = BufferSetCell;

public:
    UntypedBufferSetIterator& operator++() {};


private:
    unique_ptr<BufferSetIteratorBase> impl_;
}