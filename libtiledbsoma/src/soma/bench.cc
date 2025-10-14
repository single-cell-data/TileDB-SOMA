#include "bench.h"
#include "column_buffer.h"

namespace tiledbsoma {

void MemoryBench::release_schema(ArrowSchema* schema) {
    auto name = std::string(schema->name);

    if (schema->name != nullptr) {
        std::cout << name << " delete 'name'" << std::endl;
        free((void*)schema->name);
        schema->name = nullptr;
    }
    if (schema->format != nullptr) {
        std::cout << name << " delete 'format'" << std::endl;
        free((void*)schema->format);
        schema->format = nullptr;
    }
    if (schema->metadata != nullptr) {
        std::cout << name << " delete 'metadata'" << std::endl;
        free((void*)schema->metadata);
        schema->metadata = nullptr;
    }

    if (schema->children != nullptr) {
        std::cout << name << " delete 'children'" << std::endl;
        for (auto i = 0; i < schema->n_children; i++) {
            if (schema->children[i] != nullptr) {
                if (schema->children[i]->release != nullptr) {
                    schema->children[i]->release(schema->children[i]);
                }

                free(schema->children[i]);
                schema->children[i] = nullptr;
            }
        }

        free(schema->children);
        schema->children = nullptr;
    }

    if (schema->dictionary != nullptr) {
        std::cout << name << " delete 'dictionary'" << std::endl;
        if (schema->dictionary->release != nullptr) {
            schema->dictionary->release(schema->dictionary);
        }
        free(schema->dictionary);
        schema->dictionary = nullptr;
    }

    std::cout << "Schema deleted" << std::endl;
    schema->release = nullptr;
}

void MemoryBench::release_vector_array(ArrowArray* array) {
    auto data = static_cast<std::vector<int32_t>*>(array->private_data);
    if (data != nullptr) {
        std::cout << "Array delete 'private_data'" << std::endl;
        delete data;
    }

    if (array->buffers != nullptr) {
        std::cout << "Array delete 'buffers'" << std::endl;
        free(array->buffers);
        array->buffers = nullptr;
    }

    if (array->children != nullptr) {
        std::cout << "Array delete 'children'" << std::endl;
        for (auto i = 0; i < array->n_children; i++) {
            if (array->children[i] != nullptr) {
                if (array->children[i]->release != nullptr) {
                    array->children[i]->release(array->children[i]);
                }
                free(array->children[i]);
                array->children[i] = nullptr;
            }
        }
        free(array->children);
        array->children = nullptr;
    }

    if (array->dictionary != nullptr) {
        std::cout << "Array delete 'dictionary'" << std::endl;
        for (int64_t i = 0; i < array->dictionary->n_buffers; ++i) {
            if (array->dictionary->buffers[i] != nullptr) {
                free(const_cast<void*>(array->dictionary->buffers[i]));
                array->dictionary->buffers[i] = nullptr;
            }
        }

        array->dictionary->release(array->dictionary);
        free(array->dictionary);
        array->dictionary = nullptr;
    }

    std::cout << "Array deleted" << std::endl;
    array->release = nullptr;
}

void MemoryBench::release_pointer_array(ArrowArray* array) {
    auto data = static_cast<int32_t*>(array->private_data);
    if (data != nullptr) {
        std::cout << "Array delete 'data'" << std::endl;
        delete[] data;
    }

    if (array->buffers != nullptr) {
        std::cout << "Array delete 'buffers'" << std::endl;
        free(array->buffers);
        array->buffers = nullptr;
    }

    if (array->children != nullptr) {
        std::cout << "Array delete 'children'" << std::endl;
        for (auto i = 0; i < array->n_children; i++) {
            if (array->children[i] != nullptr) {
                if (array->children[i]->release != nullptr) {
                    array->children[i]->release(array->children[i]);
                }
                free(array->children[i]);
                array->children[i] = nullptr;
            }
        }
        free(array->children);
        array->children = nullptr;
    }

    if (array->dictionary != nullptr) {
        std::cout << "Array delete 'dictionary'" << std::endl;
        for (int64_t i = 0; i < array->dictionary->n_buffers; ++i) {
            if (array->dictionary->buffers[i] != nullptr) {
                free(const_cast<void*>(array->dictionary->buffers[i]));
                array->dictionary->buffers[i] = nullptr;
            }
        }

        array->dictionary->release(array->dictionary);
        free(array->dictionary);
        array->dictionary = nullptr;
    }

    std::cout << "Array deleted" << std::endl;
    array->release = nullptr;
}

void MemoryBench::release_column_buffer_array(ArrowArray* array) {
    auto arrow_buffer = static_cast<ArrowBuffer*>(array->private_data);
    if (arrow_buffer != nullptr) {
        std::cout << "Array delete 'data'" << std::endl;
        std::cout << "Buffer name: " << arrow_buffer->buffer_->name() << " Count: " << arrow_buffer->buffer_.use_count() << std::endl;

        delete arrow_buffer;
    }
    if (array->buffers != nullptr) {
        std::cout << "Array delete 'buffers'" << std::endl;
        free(array->buffers);
        array->buffers = nullptr;
    }

    if (array->children != nullptr) {
        std::cout << "Array delete 'children'" << std::endl;
        for (auto i = 0; i < array->n_children; i++) {
            if (array->children[i] != nullptr) {
                if (array->children[i]->release != nullptr) {
                    array->children[i]->release(array->children[i]);
                }
                free(array->children[i]);
                array->children[i] = nullptr;
            }
        }
        
        free(array->children);
        array->children = nullptr;
    }

    if (array->dictionary != nullptr) {
        std::cout << "Array delete 'dictionary'" << std::endl;
        for (int64_t i = 0; i < array->dictionary->n_buffers; ++i) {
            if (array->dictionary->buffers[i] != nullptr) {
                free(const_cast<void*>(array->dictionary->buffers[i]));
                array->dictionary->buffers[i] = nullptr;
            }
        }

        array->dictionary->release(array->dictionary);
        free(array->dictionary);
        array->dictionary = nullptr;
    }

    std::cout << "Array deleted" << std::endl;
    array->release = nullptr;
}

ArrowTable MemoryBench::allocate_vector(uint64_t size) {
    managed_unique_ptr<ArrowSchema> schema = make_managed_unique<ArrowSchema>();
    managed_unique_ptr<ArrowArray> array = make_managed_unique<ArrowArray>();

    auto data = new std::vector<int32_t>(size);
    std::cout << "Vector allocation memory bench - Int32" << std::endl;
    std::cout << "Capacity: " << data->capacity() << " Length: " << data->size() << std::endl;
    std::cout << "Data size: " << size * sizeof(int32_t) << std::endl;
    
    auto sch = schema.get();
    auto arr = array.get();

    ArrowSchemaInitFromType(sch, NANOARROW_TYPE_INT32);
    ArrowSchemaSetName(sch, "vector_int32");
    ArrowSchemaAllocateChildren(sch, 0);
    // After allocating and initializing via nanoarrow we
    // hook our custom release function in
    schema->release = &release_schema;

    ArrowArrayInitFromSchema(arr, sch, NULL);
    ArrowArrayAllocateChildren(arr, 0);
    array->length = static_cast<int64_t>(size);

    // After allocating and initializing via nanoarrow we
    // hook our custom release function in
    array->release = &release_vector_array;
    if (array->private_data != nullptr) {  // as we use nanoarrow's init
        free(array->private_data);         // free what was allocated before
    }  // assigning our ArrowBuffer pointer

    array->private_data = data;
    array->buffers = (const void**)malloc(sizeof(void*) * 2);
    array->buffers[0] = nullptr;                          
    array->buffers[1] = data->data();  

    schema->flags &= ~ARROW_FLAG_NULLABLE;

    return std::pair(std::move(array), std::move(schema));
}

ArrowTable MemoryBench::allocate_pointer(uint64_t size) {
    managed_unique_ptr<ArrowSchema> schema = make_managed_unique<ArrowSchema>();
    managed_unique_ptr<ArrowArray> array = make_managed_unique<ArrowArray>();

    auto data = new int32_t[size];
    data[0] = 1000;
    data[size - 1] = 10000; 
    std::cout << "Pointer allocation memory bench - Int32" << std::endl;
    std::cout << "Length: " << size << std::endl;
    std::cout << "Data size: " << size * sizeof(int32_t) << std::endl;

    auto sch = schema.get();
    auto arr = array.get();

    ArrowSchemaInitFromType(sch, NANOARROW_TYPE_INT32);
    ArrowSchemaSetName(sch, "pointer_int32");
    ArrowSchemaAllocateChildren(sch, 0);
    // After allocating and initializing via nanoarrow we
    // hook our custom release function in
    schema->release = &release_schema;

    ArrowArrayInitFromSchema(arr, sch, NULL);
    ArrowArrayAllocateChildren(arr, 0);
    array->length = static_cast<int64_t>(size);

    // After allocating and initializing via nanoarrow we
    // hook our custom release function in
    array->release = &release_pointer_array;
    if (array->private_data != nullptr) {  // as we use nanoarrow's init
        free(array->private_data);         // free what was allocated before
    } 

    array->private_data = data;
    array->buffers = (const void**)malloc(sizeof(void*) * 2);
    array->buffers[0] = nullptr;                          
    array->buffers[1] = data;  

    schema->flags &= ~ARROW_FLAG_NULLABLE;

    return std::pair(std::move(array), std::move(schema));
}

ArrowTable MemoryBench::allocate_column_buffer(uint64_t size) {
    managed_unique_ptr<ArrowSchema> schema = make_managed_unique<ArrowSchema>();
    managed_unique_ptr<ArrowArray> array = make_managed_unique<ArrowArray>();

    auto data = std::make_shared<ColumnBuffer>("column_buffer_int32", TILEDB_INT32, size, size * sizeof(int32_t));
    std::vector<int32_t> a;
    a.push_back(1);
    memcpy(data->data<void*>().data(), a.data(), a.size());

    std::cout << "Pointer allocation memory bench - Int32" << std::endl;
    std::cout << "Length: " << data->size() << std::endl;
    std::cout << "Data size: " << data->size() * sizeof(int32_t) << std::endl;

    auto sch = schema.get();
    auto arr = array.get();

    ArrowSchemaInitFromType(sch, NANOARROW_TYPE_INT32);
    ArrowSchemaSetName(sch, data->name().data());
    ArrowSchemaAllocateChildren(sch, 0);
    schema->release = &release_schema;

    auto arrow_buffer = new ArrowBuffer(data);

    ArrowArrayInitFromSchema(arr, sch, NULL);
    ArrowArrayAllocateChildren(arr, 0);
    array->length = static_cast<int64_t>(data->size());

    // After allocating and initializing via nanoarrow we
    // hook our custom release function in
    array->release = &release_column_buffer_array;
    if (array->private_data != nullptr) {  // as we use nanoarrow's init
        free(array->private_data);         // free what was allocated before
    }
    array->private_data = (void*)arrow_buffer;
    array->buffers = (const void**)malloc(sizeof(void*) * 2);
    array->buffers[0] = nullptr;
    array->buffers[1] = data->data<void*>().data();  
    schema->flags &= ~ARROW_FLAG_NULLABLE;

    return std::pair(std::move(array), std::move(schema));
}
}