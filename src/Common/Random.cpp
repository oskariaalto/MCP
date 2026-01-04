#include "Random.hpp"
#include <random>

thread_local std::mt19937_64 Random::engine_{std::random_device{}()};
thread_local std::uniform_real_distribution<double> Random::dist_{0.0, 1.0};
