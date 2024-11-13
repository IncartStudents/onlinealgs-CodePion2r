#include <iostream>
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

int main() {
    Welford welford;

    double data[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

    for (double value : data) {
        welford.update(value);
    }

    double mean, variance, sample_variance;
    std::tie(mean, variance, sample_variance) = welford.finalize();

    std::cout << "Mean: " << mean << std::endl;
    std::cout << "Variance: " << variance << std::endl;
    std::cout << "Sample Variance: " << sample_variance << std::endl;

    return 0;
}