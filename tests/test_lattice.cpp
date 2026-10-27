#include "kawasaki/lattice.hpp"

#include <stdexcept>
#include <cmath>
#include <random>
#include <vector>

using kawasaki::Lattice;
using kawasaki::Site;

int main() {
    for (const bool periodic_y : {false, true}) {
        for (const int length : {4, 6, 10, 16}) {
            Lattice lattice(length, periodic_y);
            std::mt19937 generator(42 + length);
            for (auto first : lattice.movable_sites())
                lattice.set(first, generator() % 2 ? 1 : -1);

            for (const auto first : lattice.movable_sites()) {
                for (const auto second : lattice.neighbours(first)) {
                    if (second.row < 0 || !lattice.movable(second) ||
                        lattice.at(first) == lattice.at(second)) continue;
                    const double before_energy = lattice.energy();
                    const double before_magnetization = lattice.magnetization();
                    const double local_change = lattice.exchange_energy_delta(first, second);
                    Lattice copy = lattice;
                    copy.exchange(first, second);
                    if (std::abs((copy.energy() - before_energy) - local_change) >= 1e-10)
                        throw std::runtime_error("Local energy delta mismatch");
                    if (copy.magnetization() != before_magnetization)
                        throw std::runtime_error("Magnetization not conserved");
                    if (std::abs(local_change - lattice.exchange_energy_delta(second, first)) >= 1e-10)
                        throw std::runtime_error("Pair selection is asymmetric");
                }
            }
        }
    }
}
