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
        if (count < window_size) {
            count++;
        } else {
            remove(oldest_value);
        }

        // Обновляем новое значение
        double delta = new_value - mean;
        mean += delta / count;
        double delta2 = new_value - mean;
        M2 += delta * delta2;

        // Сохраняем значение для удаления в следующем обновлении
        oldest_value = new_value;
    }

    std::vector<double> get_variances(const std::vector<double>& data) {
        std::vector<double> variances;

        for (size_t i = 0; i < data.size(); i++) {
            add(data[i]);

            if (i >= window_size - 1) {
                variances.push_back(get_variance());
            }
        }

        return variances;
    }

    double get_variance() const {
        if (count < 2) {
            return NAN; 
        }
        return M2 / count; 
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
    double oldest_value; // Сохранение самого старого значения для удаления
    int window_size; // Размер окна
};

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

std::vector<double> rolling_variance(const std::vector<double>& data, int window_size) {
    Welford welford(window_size);
    return welford.get_variances(data); 
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