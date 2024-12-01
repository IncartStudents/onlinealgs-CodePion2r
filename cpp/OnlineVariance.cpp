#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <random>
#include <fstream>

class Welford {
public:
    Welford(int window_size) : count(0), mean(0.0), M2(0.0), window_size(window_size) {}

    void add(double new_value) {
        count++;
        double delta = new_value - mean;
        mean += delta / count;
        double delta2 = new_value - mean;
        M2 += delta * delta2;

        // Сохраняем новое значение
        values.push_back(new_value);
        if (values.size() > window_size) {
            remove(values[0]); // Удаляем самое старое значение
            values.erase(values.begin()); // Удаляем первое значение из вектора
        }
    }

    double get_variance() const {
        if (count < 2) {
            return NAN;
        }
        return M2 / (count - 1);
    }

private:
    void remove(double old_value) {
        if (count > 0) {
            count--;
            double delta = old_value - mean;
            mean -= delta / (count + 1);
            double delta2 = old_value - mean;
            M2 -= delta * delta2;
        }
    }

    int count;      // Количество наблюдений
    double mean;    // Среднее значение
    double M2;      // Сумма квадратов отклонений от среднего
    int window_size; // Размер окна
    std::vector<double> values; // Храним значения для удаления
};

class OnlineVariance {
public:
    OnlineVariance(int window_size) : welford(window_size), index(0), window_size(window_size) {}

    void process(const std::vector<double>& data) {
        for (double value : data) {
            welford.add(value);
            if (index >= window_size - 1) {
                results_welford.push_back(welford.get_variance());
            }
            index++;
        }
    }

    void compute_two_pass_variance(const std::vector<double>& data) {
        for (size_t i = 0; i <= data.size() - window_size; ++i) {
            double mean = 0.0;
            for (size_t j = i; j < i + window_size; ++j) {
                mean += data[j];
            }
            mean /= window_size;

            double variance = 0.0;
            for (size_t j = i; j < i + window_size; ++j) {
                variance += (data[j] - mean) * (data[j] - mean);
            }
            variance /= (window_size - 1);
            results_two_pass.push_back(variance);
        }
    }

    const std::vector<double>& get_welford_results() const {
        return results_welford;
    }

    const std::vector<double>& get_two_pass_results() const {
        return results_two_pass;
    }

private:
    Welford welford;
    int index;
    int window_size;
    std::vector<double> results_welford;
    std::vector<double> results_two_pass;
};

int main() {
    const size_t data_size = 100000; 
    const int num_cycles = 100;       
    double min_value = 0.0;
    double max_value = 100.0;
    int window_size = 500;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(min_value, max_value);

    OnlineVariance online_variance(window_size);

    for (int cycle = 0; cycle < num_cycles; ++cycle) {
        std::vector<double> data(data_size);
        for (size_t i = 0; i < data_size; ++i) {
            data[i] = dis(gen);
        }

        online_variance.process(data);
        online_variance.compute_two_pass_variance(data);
    }

    std::ofstream welford_output_file("Welford_variance.txt");
    if (welford_output_file.is_open()) {
        const auto& welford_results = online_variance.get_welford_results();
        for (double variance : welford_results) {
            welford_output_file << variance << "\n";
            }
        welford_output_file.close();
    }

    std::ofstream two_pass_output_file("TwoPass_variance.txt");
    if (two_pass_output_file.is_open()) {
        const auto& two_pass_results = online_variance.get_two_pass_results();
        for (double variance : two_pass_results) {
            two_pass_output_file << variance << "\n";
        }
        two_pass_output_file.close();
    }

    return 0;
}