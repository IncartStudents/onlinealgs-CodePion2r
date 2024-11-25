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
        return M2 / (count - 1); // Используем n-1 для выборочной дисперсии
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

int main() {
    size_t data_size = 100000;
    double min_value = 0.0;
    double max_value = 100.0;
    int window_size = 500;

    // Генерация случайных данных
    std::vector<double> data(data_size);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(min_value, max_value);

    for (size_t i = 0; i < data_size; ++i) {
        data[i] = dis(gen);
    }

    // Вычисление дисперсии по скользящему окну (алгоритм Уэлфорда)
    Welford welford(window_size);
    std::vector<double> rolling_result;

    for (size_t i = 0; i < data.size(); ++i) {
        welford.add(data[i]);
        if (i >= window_size - 1) {
            rolling_result.push_back(welford.get_variance());
        }
    }

    std::ofstream rolling_output_file("rolling_variance.txt");
    if (rolling_output_file.is_open()) {
        for (double variance : rolling_result) {
            rolling_output_file << variance << "\n";
        }
        rolling_output_file.close();
    }

    // Вычисление дисперсии по двухпроходному алгоритму
    std::vector<double> two_pass_result;
    size_t n = data.size();

    if (n >= window_size) {
        for (size_t i = 0; i <= n - window_size; ++i) {
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

            two_pass_result.push_back(variance);
        }
    }

    std::ofstream two_pass_output_file("two_pass_variance.txt");
    if (two_pass_output_file.is_open()) {
        for (double variance : two_pass_result) {
            two_pass_output_file << variance << "\n";
        }
        two_pass_output_file.close();
    }

    return 0;
}