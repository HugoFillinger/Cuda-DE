#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
// Objective function example - Sphere function
double objective_function(const std::vector<double>& vec) {
    double sum = 0.0;
    for (double x : vec) {
        sum += x * x;
    }
    return sum;
}

double levy_function(const std::vector<double>& vec) {
    double w[3];
    double sum = sin(M_PI * (1 + (vec[0] - 1) / 4)) * sin(M_PI * (1 + (vec[0] - 1) / 4));
    for (int i = 0; i < 2; ++i) {
        w[i] = 1 + (vec[i] - 1) / 4;
        sum += (w[i] - 1) * (w[i] - 1) * (1 + 10 * sin(M_PI * w[i] + 1) * sin(M_PI * w[i] + 1));
    }
    w[2] = 1 + (vec[2] - 1) / 4;
    sum += (w[2] - 1) * (w[2] - 1) * (1 + sin(2 * M_PI * w[2]) * sin(2 * M_PI * w[2]));
    return sum;
}


// Differential Evolution Parameters
struct DEParams {
    int populationSize;
    int maxGenerations;
    double mutationFactor;
    double crossoverRate;
    int dimensions;
};

// Generate random double in range [min, max]
double random_double(double min, double max) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_real_distribution<double> distribution(min, max);
    return distribution(gen);
}

// Initialize population
void initialize_population(std::vector<std::vector<double>>& population, const DEParams& params) {
    for (int i = 0; i < params.populationSize; ++i) {
        std::vector<double> individual(params.dimensions);
        for (int j = 0; j < params.dimensions; ++j) {
            individual[j] = random_double(-10, 10); // Example range, adjust as needed
        }
        population.push_back(individual);
    }
}

// Differential Evolution main loop
void differential_evolution(const DEParams& params) {
    std::vector<std::vector<double>> population;
    initialize_population(population, params);

    for (int gen = 0; gen < params.maxGenerations; ++gen) {
        std::vector<std::vector<double>> new_population;

        for (int i = 0; i < params.populationSize; ++i) {
            // Mutation
            // Select three random indices different from i
            int a, b, c;
            do {
                a = rand() % params.populationSize;
            } while (a == i);

            do {
                b = rand() % params.populationSize;
            } while (b == i || b == a);

            do {
                c = rand() % params.populationSize;
            } while (c == i || c == a || c == b);

            std::vector<double> mutant(params.dimensions);
            for (int j = 0; j < params.dimensions; ++j) {
                mutant[j] = population[a][j] + params.mutationFactor * (population[b][j] - population[c][j]);
            }

            // Crossover
            std::vector<double> trial(params.dimensions);
            for (int j = 0; j < params.dimensions; ++j) {
                if (random_double(0, 1) < params.crossoverRate) {
                    trial[j] = mutant[j];
                }
                else {
                    trial[j] = population[i][j];
                }
            }

            // Selection
            if (levy_function(trial) < levy_function(population[i])) {
                new_population.push_back(trial);
            }
            else {
                new_population.push_back(population[i]);
            }
        }

        population = new_population;
    }

    // Find the best solution in the final population
    std::vector<double> best = population[0];
    for (const std::vector<double>& individual : population) {
        if (levy_function(individual) < levy_function(best)) {
            best = individual;
        }
    }

    // Output the best solution
    std::cout << "Best solution: ";
    for (double val : best) {
        std::cout << val << " ";
    }
    std::cout << "\nObjective function value: " << levy_function(best) << std::endl;
}

int main() {
    DEParams params;
    params.populationSize = 500;
    params.maxGenerations = 10000;
    params.mutationFactor = 0.5;
    params.crossoverRate = 0.7;
    params.dimensions = 3; // Ici, on définit la dimension à 3 pour la fonction de Lévy

    clock_t begin = clock();
    differential_evolution(params);
    clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    std::cout << "Time: " << elapsed_secs << std::endl;

    return 0;
}
