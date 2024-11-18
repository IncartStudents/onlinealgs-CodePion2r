#include <iostream>
#include <vector>
#include <tuple>
#include <cmath>
#include <limits>
#include <random> 
#include <fstream> 

class Welford {
public:
    Welford() : count(0), mean(0.0), M2(0.0) {}

    void update(double new_value) {
        count++;
        double delta = new_value - mean;
        mean += delta / count;
        double delta2 = new_value - mean;
        M2 += delta * delta2;
    }

    void remove(double old_value) {
        if (count > 0) {
            count--;
            double delta = old_value - mean;
            mean -= delta / (count + 1); 
            double delta2 = old_value - mean;
            M2 -= delta * delta2; 
        }
    }

    std::tuple<double, double, double> finalize() const {
        if (count < 2) {
            return std::make_tuple(NAN, NAN, NAN); // Если недостаточно данных
        } else {
            double variance = M2 / count;
            double sample_variance = M2 / (count - 1);
            return std::make_tuple(mean, variance, sample_variance);
        }
    }

private:
    int count;      // Количество наблюдений
    double mean;    // Среднее значение
    double M2;      // Сумма квадратов отклонений от среднего
};

std::vector<double> rolling_variance(const std::vector<double>& data, int window_size) {
    std::vector<double> variances;
    Welford welford;

    for (size_t i = 0; i < data.size(); i++) {
        welford.update(data[i]);

        if (i >= window_size) {
            welford.remove(data[i - window_size]);
        }

        // Добавляем дисперсию в результат, когда окно заполнено
        if (i >= window_size - 1) {
            double mean, variance, sample_variance;
            std::tie(mean, variance, sample_variance) = welford.finalize();
            variances.push_back(variance); 
        }
    }

    return variances;
}

std::vector<double> generate_random_data(size_t size, double min, double max) {
    std::vector<double> data(size);
    std::random_device rd; 
    std::mt19937 gen(rd()); 
    std::uniform_real_distribution<> dis(min, max);

    for (size_t i = 0; i < size; ++i) {
        data[i] = dis(gen); 
    }

    return data;
}

int main() {
    size_t data_size = 100000; 
    double min_value = 0.0; 
    double max_value = 100.0;
    int window_size = 500;

    std::vector<double> data = generate_random_data(data_size, min_value, max_value);

    std::vector<double> result = rolling_variance(data, window_size);

    std::ofstream output_file("rolling_variance.txt");
    if (output_file.is_open()) {
        for (double variance : result) {
            output_file << variance << "\n"; 
        }
        output_file.close(); 
    } 
    return 0;
}