#include <iostream>
#include <vector>
#include <tuple>
#include <cmath>
#include <limits>

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
            variances.push_back(variance); // Добавляем дисперсию в результат
        }
    }

    return variances;
}

int main() {
    std::vector<double> data = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    int window_size = 5;

    std::vector<double> result = rolling_variance(data, window_size);

    std::cout << "Rolling Variance: ";
    for (double variance : result) {
        std::cout << variance << " ";
    }
    std::cout << std::endl;

    return 0;
}