/*
* Code partly provided by Pablo Gonzalez Tamames
* https://github.com/tamames/toyLIFE_simulator
*/

#include "func.h"

#include <algorithm>
#include <bitset>
#include <ctime>
#include <iostream>
#include <iterator>
#include <map>
#include <random>
#include <stdexcept>
#include <string>




std::string binaryGenerator(int length, float control) {
    /**
     * Generates a random binary number as a string.
     * We can control the distribution of 1s and 0s with the control parameter.
     * @param amount The length of our binary number. Defaults to 20.
     * @param control The probability of a 1 appearing in the binary number.
     * Defaults to 0.5.
     * @return The binary number as a string.
     */
    std::string binary_string = "";
    for (int i = 0; i < length; i++) {
        float numero = (float)rand() / (float)RAND_MAX;
        binary_string += numero < control ? '1' : '0';
    }
    return binary_string;
}

// std::vector<float> createRandomArray(int size, int max /*= 8*/) {
//     /**
//      * Creates a vector and fill it with random numbers.
//      * @param size The size of the array.
//      * @param max The maximum value of the random integers. Defaults to 8.
//      */
//     // std::random_device rd;
//     // std::mt19937 gen(rd());
//     // std::uniform_int_distribution<> dis(1, max);

//     std::uniform_real_distribution<double> distribution(0.0, max);

//     std::vector<float> energies(size);
//     for (int i = 0; i < size; ++i) {
//         energies[i] = distribution(GENERATOR);
//     }
//     return energies;
// }

std::map<std::string, int> fromList2Map(std::vector<std::string> list) {
    /**
     * Converts a list of binary strings to a map of binary strings and their
     *  respective counts.
     * @param list The list of binary strings.
     * @return The map of binary strings and their respective indexes.
     */

    std::map<std::string, int> map;
    for (int i = 0; i < list.size(); ++i)
        map[list[i]] += 1;

    return map;
}

std::vector<std::string> decimal2Binary(int numbers) {
    /**
     * Creates a vector with the binary representation of the numbers from 0 to
     * max_number - 1.
     * @param numbers The amount of numbers to be converted.
     * @return The vector with the binary representation of the numbers.
     */
    std::vector<std::string> binary_numbers(numbers);
    for (int i = 0; i < numbers; ++i) {
        binary_numbers[i] = std::bitset<8>(i).to_string();
    }

    return binary_numbers;
}

void printMap(mapa_prot& mapa) {
    /**
     * Prints a map.
     * @param mapa The map to be printed.
     */
    for (auto const& pair : mapa) {
        std::cout << "{" << pair.first << ": " << pair.second << "}\n";
    }
}

void printMap(mapa_met& mapa) {
    /**
     * Prints a map.
     * @param mapa The map to be printed.
     */
    for (auto const& pair : mapa) {
        std::cout << "{" << pair.first << ": " << pair.second << "}\n";
    }
}

void printMap(mapa_dim& mapa) {
    /**
     * Prints a map.
     * @param mapa The map to be printed.
     */
    for (auto const& pair : mapa) {
        std::cout << "{" << pair.first.id << ": " << pair.second << "}\n";
    }
}

void printMap(mapa_owm& mapa) {
    /**
     * Prints a map.
     * @param mapa The map to be printed.
     */

    for (auto const& pair : mapa) {
        std::cout << "{Prot: " << pair.first.prot
                  << ", Dim: " << pair.first.dim.id
                  << ", Met: " << pair.first.met << ": " << pair.second
                  << "}\n";
    }
}

std::string fromMapToString(mapa_prot& mapa) {
    /**
     * Converts a map to a string. The prot is an integer.
     * @param mapa The map to be converted.
     * @return The string.
     */

    if (mapa.empty())
        return "";
    std::string out = "{";
    for (auto const& pair : mapa) {
        out += std::to_string(pair.first) + ": " + std::to_string(pair.second) +
               ", ";
    }
    // to remove the last comma and space
    out.pop_back();
    out.pop_back();

    out += "}";
    return out;
}

std::string fromMapToString(mapa_met& mapa) {
    /**
     * Converts a map to a string. In this case the met is already a string.
     * @param mapa The map to be converted.
     * @return The string.
     */

    if (mapa.empty())
        return "";
    std::string out = "{";
    for (auto const& pair : mapa) {
        out += pair.first + ": " + std::to_string(pair.second) + ", ";
    }
    out.pop_back();
    out.pop_back();
    out += "}";
    return out;
}

std::string fromMapToString(mapa_dim& mapa) {
    /**
     * Converts a map to a string. In this case the dim is an object so we take
     * the id.
     * @param mapa The map to be converted.
     * @return The string.
     */

    if (mapa.empty())
        return "";
    std::string out = "{";
    for (auto const& pair : mapa) {
        out += std::to_string(pair.first.id) + ": " +
               std::to_string(pair.second) + ", ";
    }
    out.pop_back();
    out.pop_back();
    out += "}";
    return out;
}

std::string currentTime() {

    /**
     * Returns the current hour, minute and second in a string.
     * @return The string with the current hour, minute and second.
     */
    std::time_t now = time(0);
    std::tm* ltm = localtime(&now);

    // Format the time as HH:MM:SS
    std::stringstream ss;
    ss << std::put_time(ltm, "%H:%M:%S");

    return ss.str();
}


// std::string mutate(std::string binString) {
//     /**
//      * Mutates a binary string.
//      * @param binString The binary string to be mutated.
//      * @return The mutated binary string.
//      */

//     if (MUTATION_PROBABILITY < 0.0f || MUTATION_PROBABILITY > 1.0f) {
//         std::cout << "The probability must be between 0 and 1\n";
//         throw std::invalid_argument("The probability must be between 0 and 1");
//     }

//     std::string newString(binString.size(), ' ');
//     for (int i = 0; i < binString.size(); ++i) {
//         if (std::generate_canonical<double, 5>(GENERATOR) <
//             MUTATION_PROBABILITY)
//             newString[i] = binString[i] == '0' ? '1' : '0';
//         else
//             newString[i] = binString[i];
//     }
//     return newString;
// }


// std::vector<bool> intToBoolVector(int num, int size) {
//     std::vector<bool> result(size, false);
//     for (int i = size - 1; i >= 0 && num > 0; --i) {
//         result[i] = num % 2;
//         num /= 2;
//     }
//     return result; // Ensure the function returns the result
// }
