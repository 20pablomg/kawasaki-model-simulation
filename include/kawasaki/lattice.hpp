#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

namespace kawasaki {

struct Site {
    int row;
    int col;
    bool operator==(const Site& rhs) const noexcept {
        return row == rhs.row && col == rhs.col;
    }
};

class Lattice {
public:
    explicit Lattice(int length, bool periodic_y = false)
        : length_(length), periodic_y_(periodic_y), spins_(length * length, -1) {
        if (length < 4 || length % 2 != 0) {
            throw std::invalid_argument("L must be even and at least 4");
        }
        if (!periodic_y_) {
            for (int col = 0; col < length_; ++col) {
                set({0, col}, -1);
                set({length_ - 1, col}, 1);
            }
        }
    }

    int length() const noexcept { return length_; }
    int sites() const noexcept { return length_ * length_; }
    bool periodic_y() const noexcept { return periodic_y_; }
    int at(Site site) const { return spins_.at(index(site)); }
    void set(Site site, int spin) {
        if (spin != -1 && spin != 1) throw std::invalid_argument("Spin must be -1 or +1");
        spins_.at(index(site)) = spin;
    }
    bool movable(Site site) const noexcept {
        return periodic_y_ || (site.row > 0 && site.row < length_ - 1);
    }
    const std::vector<int>& spins() const noexcept { return spins_; }

    std::vector<Site> movable_sites() const {
        std::vector<Site> sites_list;
        sites_list.reserve(sites());
        for (int row = 0; row < length_; ++row)
            for (int col = 0; col < length_; ++col)
                if (movable({row, col})) sites_list.push_back({row, col});
        return sites_list;
    }

    // Returns all nearest neighbours, including the fixed boundary rows.
    std::array<Site, 4> neighbours(Site site) const {
        const int up = site.row == 0 ? (periodic_y_ ? length_ - 1 : -1) : site.row - 1;
        const int down = site.row == length_ - 1 ? (periodic_y_ ? 0 : -1) : site.row + 1;
        return {{{site.row, (site.col + 1) % length_},
                 {site.row, (site.col + length_ - 1) % length_},
                 {down, site.col}, {up, site.col}}};
    }

    double energy() const {
        double energy_value = 0.0;
        for (int row = 0; row < length_; ++row) {
            for (int col = 0; col < length_; ++col) {
                const int spin = at({row, col});
                energy_value -= spin * at({row, (col + 1) % length_});
                if (row + 1 < length_) energy_value -= spin * at({row + 1, col});
                else if (periodic_y_) energy_value -= spin * at({0, col});
            }
        }
        return energy_value;
    }

    // Calculate only the bonds affected by an exchange; count each bond once.
    double exchange_energy_delta(Site first, Site second) const {
        if (at(first) == at(second)) return 0.0;
        std::vector<std::pair<int, int>> bonds;
        bonds.reserve(8);
        for (const auto site : {first, second}) {
            for (const auto neighbour : neighbours(site)) {
                if (neighbour.row < 0) continue;
                const int first_id = index(site);
                const int second_id = index(neighbour);
                const std::pair<int, int> edge = std::minmax(first_id, second_id);
                if (std::find(bonds.begin(), bonds.end(), edge) == bonds.end())
                    bonds.push_back(edge);
            }
        }
        double delta = 0.0;
        const int first_index = index(first);
        const int second_index = index(second);
        for (const auto& [i, j] : bonds) {
            const int before_i = spins_[i];
            const int before_j = spins_[j];
            const int after_i = (i == first_index || i == second_index) ? -before_i : before_i;
            const int after_j = (j == first_index || j == second_index) ? -before_j : before_j;
            delta += -after_i * after_j + before_i * before_j;
        }
        return delta;
    }

    void exchange(Site first, Site second) {
        if (!movable(first) || !movable(second))
            throw std::invalid_argument("Cannot exchange a fixed boundary spin");
        std::swap(spins_[index(first)], spins_[index(second)]);
    }

    double magnetization() const {
        return static_cast<double>(std::accumulate(spins_.begin(), spins_.end(), 0)) / sites();
    }

    std::pair<double, double> domain_magnetizations(double m0 = 0.0) const {
        // For nonzero magnetization, assign (1+m0)/2 of rows to the + domain.
        const int bottom_rows = std::clamp(
            static_cast<int>(std::llround((1.0 + m0) * length_ / 2.0)), 1, length_ - 1);
        const int first_bottom = length_ - bottom_rows;
        int upper = 0, lower = 0;
        for (int row = 0; row < length_; ++row)
            for (int col = 0; col < length_; ++col)
                (row < first_bottom ? upper : lower) += at({row, col});
        return {static_cast<double>(upper) / (first_bottom * length_),
                static_cast<double>(lower) / (bottom_rows * length_)};
    }

    double domain_order(double m0 = 0.0) const {
        const auto [upper, lower] = domain_magnetizations(m0);
        const int bottom_rows = std::clamp(
            static_cast<int>(std::llround((1.0 + m0) * length_ / 2.0)), 1, length_ - 1);
        const double bottom_weight = static_cast<double>(bottom_rows) / length_;
        return (1.0 - bottom_weight) * std::abs(upper) + bottom_weight * std::abs(lower);
    }

    double signed_domain_order(double m0 = 0.0) const {
        const auto [upper, lower] = domain_magnetizations(m0);
        return (lower - upper) / 2.0;
    }

private:
    int index(Site site) const {
        if (site.row < 0 || site.row >= length_ || site.col < 0 || site.col >= length_)
            throw std::out_of_range("Lattice site outside boundaries");
        return site.row * length_ + site.col;
    }
    int length_;
    bool periodic_y_;
    std::vector<int> spins_;
};

}  // namespace kawasaki
