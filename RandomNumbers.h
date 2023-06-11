#ifndef RANDOMNUMBERS_H
#define RANDOMNUMBERS_H

#include <random>

template<typename T>
class RandomGenerator;

template<>
class RandomGenerator<int> {
private:
    std::mt19937 gen;
    std::uniform_int_distribution<> dis;

public:
    // Default constructor
    RandomGenerator() : gen(std::random_device()()), dis(0, 100) {} 

    // Constructor with parameters
    RandomGenerator(int seed, int lower_bound, int upper_bound)
        : gen(seed), dis(lower_bound, upper_bound) {}

    // Generate a random number
    int getRandomNumber() {
        return dis(gen);
    }
};

template<>
class RandomGenerator<float> {
private:
    std::mt19937 gen;
    std::uniform_real_distribution<float> dis;

public:
    // Default constructor
    RandomGenerator() : gen(std::random_device()()), dis(0.0f, 1.0f) {} 

    // Constructor with parameters
    RandomGenerator(int seed, float lower_bound, float upper_bound)
        : gen(seed), dis(lower_bound, upper_bound) {}

    // Generate a random number
    float getRandomNumber() {
        return dis(gen);
    }
};
#endif
