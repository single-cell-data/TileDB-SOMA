/**
 * @file tiledb/common/thread_pool/test/unit_thread_pool.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * Tests the `ThreadPool` class.
 */

#include <catch2/catch_test_macros.hpp>

#include <stdint.h>
#include <atomic>
#include <iostream>
#include <random>

#include "thread_pool/thread_pool.h"

using namespace tiledbsoma;

size_t random_ms(size_t max = 3) {
    thread_local static uint64_t generator_seed = std::hash<std::thread::id>()(std::this_thread::get_id());
    thread_local static std::mt19937_64 generator(generator_seed);
    std::uniform_int_distribution<size_t> distribution(0, max);
    return distribution(generator);
}

TEST_CASE("ThreadPool: Test empty", "[threadpool]") {
    for (int i = 0; i < 10; i++) {
        ThreadPool pool{4};
    }
}

TEST_CASE("ThreadPool: Test single thread", "[threadpool]") {
    std::atomic<int> result = 0;  // needs to be atomic b/c scavenging thread
                                  // can run in addition to thread pool
    std::vector<ThreadPool::Task> results;
    ThreadPool pool{1};

    for (int i = 0; i < 100; i++) {
        ThreadPool::Task task = pool.execute([&result]() {
            result++;
            return Status::Ok();
        });

        REQUIRE(task.valid());

        results.emplace_back(std::move(task));
    }
    REQUIRE(pool.wait_all(results).ok());
    REQUIRE(result == 100);
}

TEST_CASE("ThreadPool: Test multiple threads", "[threadpool]") {
    std::atomic<int> result(0);
    std::vector<ThreadPool::Task> results;
    ThreadPool pool{4};
    for (int i = 0; i < 100; i++) {
        results.push_back(pool.execute([&result]() {
            result++;
            return Status::Ok();
        }));
    }
    REQUIRE(pool.wait_all(results).ok());
    REQUIRE(result == 100);
}

TEST_CASE("ThreadPool: Test wait status", "[threadpool]") {
    std::atomic<int> result(0);
    std::vector<ThreadPool::Task> results;
    ThreadPool pool{4};
    for (int i = 0; i < 100; i++) {
        results.push_back(pool.execute([&result, i]() {
            result++;
            return i == 50 ? Status_Error("Generic error") : Status::Ok();
        }));
    }
    REQUIRE(!pool.wait_all(results).ok());
    REQUIRE(result == 100);
}

TEST_CASE("ThreadPool: Test no wait", "[threadpool]") {
    {
        ThreadPool pool{4};
        std::atomic<int> result(0);
        for (int i = 0; i < 5; i++) {
            ThreadPool::Task task = pool.execute([&result]() {
                result++;
                std::this_thread::sleep_for(std::chrono::milliseconds(random_ms(1000)));
                return Status::Ok();
            });
            REQUIRE(task.valid());
        }
        // There may be an error logged when the pool is destroyed if there are
        // outstanding tasks, but everything should still complete.
    }
}

// Defer this test, pending final design of ThreadPool lifecycle
// TEST_CASE("ThreadPool: Test construction of empty threadpool",
// "[threadpool]") {
//  try {
//    ThreadPool pool;
//  } catch (std::exception& e) {
//    auto st = Status_ThreadPoolError(e.what());
//    REQUIRE(!st.ok());
//  }
// }

TEST_CASE("ThreadPool: Test recursion, simplest case", "[threadpool]") {
    ThreadPool pool{1};

    std::atomic<int> result(0);

    std::vector<ThreadPool::Task> tasks;
    auto a = pool.execute([&pool, &result]() {
        std::vector<ThreadPool::Task> tasks;
        auto b = pool.execute([&result]() {
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
            ++result;
            return Status::Ok();
        });
        REQUIRE(b.valid());
        tasks.emplace_back(std::move(b));
        pool.wait_all(tasks);
        return Status::Ok();
    });
    REQUIRE(a.valid());
    tasks.emplace_back(std::move(a));
    pool.wait_all(tasks);
    REQUIRE(result == 1);
}

TEST_CASE("ThreadPool: Test recursion", "[threadpool]") {
    size_t num_threads = 0;
    SECTION("- One thread") {
        num_threads = 1;
    }

    SECTION("- Two threads") {
        num_threads = 2;
    }

    SECTION("- Ten threads") {
        num_threads = 10;
    }

    ThreadPool pool{num_threads};

    // Test recursive execute-and-wait.
    std::atomic<int> result(0);
    const size_t num_tasks = 100;
    const size_t num_nested_tasks = 10;
    std::vector<ThreadPool::Task> tasks;
    for (size_t i = 0; i < num_tasks; ++i) {
        auto task = pool.execute([&]() {
            std::vector<ThreadPool::Task> inner_tasks;
            for (size_t j = 0; j < num_nested_tasks; ++j) {
                auto inner_task = pool.execute([&]() {
                    std::this_thread::sleep_for(std::chrono::milliseconds(random_ms()));
                    ++result;
                    return Status::Ok();
                });

                inner_tasks.emplace_back(std::move(inner_task));
            }

            pool.wait_all(inner_tasks).ok();
            return Status::Ok();
        });

        REQUIRE(task.valid());
        tasks.emplace_back(std::move(task));
    }
    REQUIRE(pool.wait_all(tasks).ok());
    REQUIRE(result == (num_tasks * num_nested_tasks));

    // Test a top-level execute-and-wait with async-style inner tasks.
    std::condition_variable cv;
    std::mutex cv_mutex;
    tasks.clear();
    for (size_t i = 0; i < num_tasks; ++i) {
        auto task = pool.execute([&]() {
            for (size_t j = 0; j < num_nested_tasks; ++j) {
                pool.execute([&]() {
                    std::this_thread::sleep_for(std::chrono::milliseconds(random_ms()));

                    std::unique_lock<std::mutex> ul(cv_mutex);
                    if (--result == 0) {
                        cv.notify_all();
                    }
                    return Status::Ok();
                });
            }

            return Status::Ok();
        });

        REQUIRE(task.valid());
        tasks.emplace_back(std::move(task));
    }

    REQUIRE(pool.wait_all(tasks).ok());

    // Wait all inner tasks to complete.
    std::unique_lock<std::mutex> ul(cv_mutex);
    while (result > 0)
        cv.wait(ul);
}

TEST_CASE("ThreadPool: Test recursion, two pools", "[threadpool]") {
    size_t num_threads = 0;

    // Shorten original test by reducing repeat count and skipping the 1 and 2
    // thread tests.
    int repeat_count = 10;

    /*
        SECTION("- One thread") {
            num_threads = 1;
        }

        SECTION("- Two threads") {
            num_threads = 2;
        }
    */

    SECTION("- Ten threads") {
        num_threads = 10;
    }

    SECTION("- Twenty threads") {
        num_threads = 20;
    }

    ThreadPool pool_a{num_threads};
    ThreadPool pool_b{num_threads};

    // This test logic is relatively inexpensive, run it 50 times
    // to increase the chance of encountering race conditions.
    for (int t = 0; t < repeat_count; ++t) {
        // Test recursive execute-and-wait.
        std::atomic<int> result(0);
        const size_t num_tasks_a = 10;
        const size_t num_tasks_b = 10;
        const size_t num_tasks_c = 10;
        std::vector<ThreadPool::Task> tasks_a;
        for (size_t i = 0; i < num_tasks_a; ++i) {
            auto task_a = pool_a.execute([&]() {
                std::vector<ThreadPool::Task> tasks_b;
                for (size_t j = 0; j < num_tasks_c; ++j) {
                    auto task_b = pool_b.execute([&]() {
                        std::vector<ThreadPool::Task> tasks_c;
                        for (size_t k = 0; k < num_tasks_b; ++k) {
                            auto task_c = pool_a.execute([&result]() {
                                std::this_thread::sleep_for(std::chrono::milliseconds(random_ms()));
                                ++result;
                                return Status::Ok();
                            });

                            tasks_c.emplace_back(std::move(task_c));
                        }

                        pool_a.wait_all(tasks_c);
                        return Status::Ok();
                    });

                    tasks_b.emplace_back(std::move(task_b));
                }

                pool_b.wait_all(tasks_b).ok();
                return Status::Ok();
            });

            REQUIRE(task_a.valid());
            tasks_a.emplace_back(std::move(task_a));
        }
        REQUIRE(pool_a.wait_all(tasks_a).ok());
        REQUIRE(result == (num_tasks_a * num_tasks_b * num_tasks_c));

        // Test a top-level execute-and-wait with async-style inner tasks.
        std::condition_variable cv;
        std::mutex cv_mutex;
        tasks_a.clear();
        for (size_t i = 0; i < num_tasks_a; ++i) {
            auto task_a = pool_a.execute([&]() {
                std::vector<ThreadPool::Task> tasks_b;
                for (size_t j = 0; j < num_tasks_b; ++j) {
                    auto task_b = pool_b.execute([&]() {
                        std::vector<ThreadPool::Task> tasks_c;
                        for (size_t k = 0; k < num_tasks_c; ++k) {
                            auto task_c = pool_a.execute([&]() {
                                std::this_thread::sleep_for(std::chrono::milliseconds(random_ms()));
                                if (--result == 0) {
                                    std::unique_lock<std::mutex> ul(cv_mutex);
                                    cv.notify_all();
                                }
                                return Status::Ok();
                            });

                            tasks_c.emplace_back(std::move(task_c));
                        }

                        pool_a.wait_all(tasks_c);
                        return Status::Ok();
                    });

                    tasks_b.emplace_back(std::move(task_b));
                }

                pool_b.wait_all(tasks_b).ok();
                return Status::Ok();
            });

            REQUIRE(task_a.valid());
            tasks_a.emplace_back(std::move(task_a));
        }
        REQUIRE(pool_a.wait_all(tasks_a).ok());

        // Wait all inner tasks to complete.
        std::unique_lock<std::mutex> ul(cv_mutex);
        while (result > 0)
            cv.wait(ul);
    }
}

TEST_CASE("ThreadPool: Test Exceptions", "[threadpool]") {
    std::atomic<int> result(0);
    ThreadPool pool{7};

    Status unripe_banana_status = Status_TaskError("Caught msg: Unripe banana");
    Status unbaked_potato_status = Status_TileError("Unbaked potato");

    SECTION("One task error exception") {
        std::vector<ThreadPool::Task> results;

        for (int i = 0; i < 207; ++i) {
            results.push_back(pool.execute([&result]() {
                auto tmp = result++;
                if (tmp == 13) {
                    throw(std::string("Unripe banana"));
                }
                return Status::Ok();
            }));
        }

        REQUIRE(pool.wait_all(results).to_string() == unripe_banana_status.to_string());
        REQUIRE(result == 207);
    }

    SECTION("One tile error exception") {
        std::vector<ThreadPool::Task> results;

        for (int i = 0; i < 207; ++i) {
            results.push_back(pool.execute([&result, &unbaked_potato_status]() {
                auto tmp = result++;
                if (tmp == 31) {
                    throw(unbaked_potato_status);
                }
                return Status::Ok();
            }));
        }

        REQUIRE(pool.wait_all(results).to_string() == unbaked_potato_status.to_string());
        REQUIRE(result == 207);
    }

    SECTION("Two exceptions") {
        std::vector<ThreadPool::Task> results;

        for (int i = 0; i < 207; ++i) {
            results.push_back(pool.execute([&result]() {
                auto tmp = result++;
                if (tmp == 13) {
                    throw(std::string("Unripe banana"));
                }
                if (tmp == 31) {
                    throw(Status_TileError("Unbaked potato"));
                }

                return Status::Ok();
            }));
        }

        auto pool_status = pool.wait_all(results);
        REQUIRE(
            ((pool_status.to_string() == unripe_banana_status.to_string()) ||
             (pool_status.to_string() == unbaked_potato_status.to_string())));
        REQUIRE(result == 207);
    }

    SECTION("Two exceptions reverse order") {
        std::vector<ThreadPool::Task> results;

        for (int i = 0; i < 207; ++i) {
            results.push_back(pool.execute([&result]() {
                auto tmp = result++;
                if (tmp == 31) {
                    throw(std::string("Unripe banana"));
                }
                if (tmp == 13) {
                    throw(Status_TileError("Unbaked potato"));
                }

                return Status::Ok();
            }));
        }

        auto pool_status = pool.wait_all(results);
        REQUIRE(
            ((pool_status.to_string() == unripe_banana_status.to_string()) ||
             (pool_status.to_string() == unbaked_potato_status.to_string())));
        REQUIRE(result == 207);
    }

    SECTION("Two exceptions strict order") {
        std::vector<ThreadPool::Task> results;

        for (int i = 0; i < 207; ++i) {
            results.push_back(pool.execute([i, &result]() {
                result++;
                if (i == 13) {
                    throw(std::string("Unripe banana"));
                }
                if (i == 31) {
                    throw(Status_TileError("Unbaked potato"));
                }

                return Status::Ok();
            }));
        }

        REQUIRE(pool.wait_all(results).to_string() == unripe_banana_status.to_string());
        REQUIRE(result == 207);
    }

    SECTION("Two exceptions strict reverse order") {
        std::vector<ThreadPool::Task> results;

        for (int i = 0; i < 207; ++i) {
            results.push_back(pool.execute([i, &result]() {
                ++result;
                if (i == 31) {
                    throw(std::string("Unripe banana"));
                }
                if (i == 13) {
                    throw(Status_TileError("Unbaked potato"));
                }

                return Status::Ok();
            }));
        }

        REQUIRE(pool.wait_all(results).to_string() == unbaked_potato_status.to_string());
        REQUIRE(result == 207);
    }
}
