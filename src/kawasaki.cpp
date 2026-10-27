#include "kawasaki/lattice.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;
using kawasaki::Lattice;
using kawasaki::Site;

namespace {
struct Config {
    int length = 24;
    double temperature = 2.3;
    double magnetization = 0.0;
    long sweeps = 3000;
    long burn_in = 1000;
    long stride = 10;
    std::uint64_t seed = 2024;
    bool periodic_y = false;
    std::string output = "outputs/run";
};

struct Moments {
    long count = 0;
    double energy = 0.0, energy_squared = 0.0;
    double order = 0.0, order_squared = 0.0;
    double magnetization = 0.0, signed_order = 0.0;

    void accumulate(double e, double o, double m, double signed_o) {
        ++count;
        energy += e;
        energy_squared += e * e;
        order += o;
        order_squared += o * o;
        magnetization += m;
        signed_order += signed_o;
    }
};

void print_help() {
    std::cout << "2D Ising model with Kawasaki spin-exchange dynamics\n"
              << "Options:\n"
              << "  --L INT             Even lattice size >= 4 (default 24)\n"
              << "  --T FLOAT           Temperature > 0 (default 2.3)\n"
              << "  --m0 FLOAT          Requested magnetization [-1, 1] (default 0)\n"
              << "  --sweeps INT        Monte Carlo sweeps (default 3000)\n"
              << "  --burn INT          Burn-in sweeps (default 1000)\n"
              << "  --stride INT        Sample every INT sweeps after burn-in (default 10)\n"
              << "  --seed INT          Random seed (default 2024)\n"
              << "  --periodic-y 0|1    Periodic y instead of fixed boundary spins\n"
              << "  --output PATH       Output filename prefix (default outputs/run)\n"
              << "  --help              Show help\n";
}

Config parse_args(int argc, char** argv) {
    Config config;
    for (int index = 1; index < argc; ++index) {
        const std::string key = argv[index];
        if (key == "--help" || key == "-h") { print_help(); std::exit(0); }
        if (index + 1 >= argc) throw std::invalid_argument("Missing value for " + key);
        const std::string value = argv[++index];
        if (key == "--L") config.length = std::stoi(value);
        else if (key == "--T") config.temperature = std::stod(value);
        else if (key == "--m0") config.magnetization = std::stod(value);
        else if (key == "--sweeps") config.sweeps = std::stol(value);
        else if (key == "--burn") config.burn_in = std::stol(value);
        else if (key == "--stride") config.stride = std::stol(value);
        else if (key == "--seed") config.seed = std::stoull(value);
        else if (key == "--output") config.output = value;
        else if (key == "--periodic-y") {
            if (value != "0" && value != "1") throw std::invalid_argument("--periodic-y must be 0 or 1");
            config.periodic_y = value == "1";
        } else throw std::invalid_argument("Unknown option: " + key);
    }
    if (config.length < 4 || config.length % 2 != 0 || !std::isfinite(config.temperature) ||
        config.temperature <= 0 || !std::isfinite(config.magnetization) ||
        std::abs(config.magnetization) > 1 || config.sweeps <= 0 || config.burn_in < 0 ||
        config.burn_in >= config.sweeps || config.stride <= 0)
        throw std::invalid_argument("Invalid simulation parameters; run --help");
    return config;
}

void ensure_output_directory(const std::string& prefix) {
    const auto parent = fs::path(prefix).parent_path();
    if (!parent.empty()) fs::create_directories(parent);
}

void validate_stream(const std::ofstream& stream, const std::string& path) {
    if (!stream) throw std::runtime_error("Cannot write to " + path);
}

void simulate(const Config& config) {
    ensure_output_directory(config.output);
    std::mt19937_64 engine(config.seed);
    std::uniform_real_distribution<double> probability(0.0, 1.0);
    std::uniform_int_distribution<int> direction(0, 3);
    Lattice lattice(config.length, config.periodic_y);
    std::vector<Site> mobile = lattice.movable_sites();

    // The requested magnetization is rounded to the nearest feasible spin count.
    const int plus = std::clamp(static_cast<int>(std::llround(
        (1.0 + config.magnetization) * mobile.size() / 2.0)), 0, static_cast<int>(mobile.size()));
    std::shuffle(mobile.begin(), mobile.end(), engine);
    for (std::size_t i = 0; i < mobile.size(); ++i)
        lattice.set(mobile[i], static_cast<int>(i) < plus ? 1 : -1);
    const double realized_magnetization = lattice.magnetization();
    std::uniform_int_distribution<std::size_t> choose_site(0, mobile.size() - 1);

    std::ofstream series(config.output + "_timeseries.csv");
    validate_stream(series, config.output + "_timeseries.csv");
    series << "sweep,energy,magnetization,domain_order,domain_signed\n";
    series << std::setprecision(14);

    Moments moments;
    std::vector<double> density(config.length, 0.0);
    std::uint64_t attempted = 0, accepted = 0;

    for (long sweep = 0; sweep < config.sweeps; ++sweep) {
        for (int attempt = 0; attempt < lattice.sites(); ++attempt) {
            const Site first = mobile[choose_site(engine)];
            Site second = first;
            switch (direction(engine)) {
                case 0: second.col = (first.col + 1) % config.length; break;
                case 1: second.col = (first.col + config.length - 1) % config.length; break;
                case 2:
                    if (first.row + 1 < config.length) second.row = first.row + 1;
                    else if (config.periodic_y) second.row = 0;
                    break;
                case 3:
                    if (first.row > 0) second.row = first.row - 1;
                    else if (config.periodic_y) second.row = config.length - 1;
                    break;
            }
            ++attempted;
            if (!lattice.movable(second) || second == first ||
                lattice.at(first) == lattice.at(second)) continue;
            const double change = lattice.exchange_energy_delta(first, second);
            if (change <= 0.0 || probability(engine) < std::exp(-change / config.temperature)) {
                lattice.exchange(first, second);
                ++accepted;
            }
        }
        if (sweep < config.burn_in || (sweep - config.burn_in) % config.stride != 0)
            continue;
        const double energy = lattice.energy();
        const double order = lattice.domain_order(realized_magnetization);
        const double signed_order = lattice.signed_domain_order(realized_magnetization);
        const double magnetization = lattice.magnetization();
        moments.accumulate(energy, order, magnetization, signed_order);
        for (int row = 0; row < config.length; ++row)
            for (int col = 0; col < config.length; ++col)
                if (lattice.at({row, col}) == 1) density[row] += 1.0 / config.length;
        series << sweep + 1 << ',' << energy << ',' << magnetization << ',' << order
               << ',' << signed_order << '\n';
    }

    std::ofstream snapshot(config.output + "_final.csv");
    validate_stream(snapshot, config.output + "_final.csv");
    for (int row = 0; row < config.length; ++row) {
        for (int col = 0; col < config.length; ++col) {
            if (col) snapshot << ',';
            snapshot << lattice.at({row, col});
        }
        snapshot << '\n';
    }

    std::ofstream profile(config.output + "_density.csv");
    validate_stream(profile, config.output + "_density.csv");
    profile << "y,mean_plus_density\n" << std::setprecision(14);
    for (int row = 0; row < config.length; ++row)
        profile << row << ',' << density[row] / moments.count << '\n';

    const double n = static_cast<double>(moments.count);
    const double number_sites = static_cast<double>(lattice.sites());
    const double energy_variance = std::max(0.0, moments.energy_squared / n -
                                                   std::pow(moments.energy / n, 2));
    const double order_variance = std::max(0.0, moments.order_squared / n -
                                                  std::pow(moments.order / n, 2));
    std::ofstream summary(config.output + "_summary.csv");
    validate_stream(summary, config.output + "_summary.csv");
    summary << "L,T,m0_requested,m0_realized,seed,sweeps,burn,stride,n_samples,"
               "energy_per_site,specific_heat,domain_order,domain_order_variance,"
               "domain_response_proxy,mean_magnetization,mean_domain_signed,"
               "acceptance_rate,periodic_y\n";
    summary << std::setprecision(14) << config.length << ',' << config.temperature << ','
            << config.magnetization << ',' << realized_magnetization << ',' << config.seed << ','
            << config.sweeps << ',' << config.burn_in << ',' << config.stride << ','
            << moments.count << ',' << moments.energy / (n * number_sites) << ','
            << energy_variance / (number_sites * config.temperature * config.temperature) << ','
            << moments.order / n << ',' << order_variance << ','
            << number_sites * order_variance / config.temperature << ','
            << moments.magnetization / n << ',' << moments.signed_order / n << ','
            << static_cast<double>(accepted) / attempted << ',' << config.periodic_y << '\n';

    if (lattice.magnetization() != realized_magnetization)
        throw std::runtime_error("Invariant failed: total magnetization changed");
    std::cout << "Results written to " << config.output
              << "_{summary,timeseries,final,density}.csv\n";
}
} // namespace

int main(int argc, char** argv) {
    try {
        simulate(parse_args(argc, argv));
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
