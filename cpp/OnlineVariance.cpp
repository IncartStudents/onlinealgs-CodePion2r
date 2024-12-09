#include <iostream>
#include <cmath>
#include <limits>
#include <random>
#include <fstream>
#include <vector>
#include <stdexcept>

void write_results_to_file(const std::string& filename, const double* results, size_t result_count);

class CircularBuffer {
public:
    CircularBuffer(size_t size) : size(size), head(0), count(0) {
        buffer = new double[size];
    }

    ~CircularBuffer() {
        delete[] buffer;
    }

    void add(double value) {
        buffer[head] = value;
        head = (head + 1) % size;
        if (count < size) count++;
    }

    void add_batch(const double* values, size_t batch_size) {
        for (size_t i = 0; i < batch_size; ++i) {
            add(values[i]);
        }
    }

    double get(size_t index) const {
        if (index >= count) throw std::out_of_range("Index out of range");
        return buffer[(head + size - count + index) % size];
    }

    size_t get_count() const {
        return count;
    }

private:
    size_t size; // Максимальный размер буфера
    double* buffer; // Массив для хранения данных
    size_t head; // Индекс следующего элемента для записи
    size_t count; // Текущее количество элементов в буфере
};

class Welford {
public:
    Welford() : count(0), mean(0.0), M2(0.0) {}

    void add(double new_value) {
        count++;
        double delta = new_value - mean;
        mean += delta / count;
        double delta2 = new_value - mean;
        M2 += delta * delta2;
    }

    void add_batch(const double* values, size_t batch_size) {
        for (size_t i = 0; i < batch_size; ++i) {
            add(values[i]);
        }
    }

    double get_variance() const {
        if (count < 2) {
            return NAN; 
        }
        return M2 / (count - 1);
    }

private:
    int count;      // Количество наблюдений
    double mean;    // Среднее значение
    double M2;      // Сумма квадратов отклонений от среднего
};

double compute_two_pass_variance(const CircularBuffer& buffer) {
    double mean = 0.0;
    size_t count = buffer.get_count();

    for (size_t j = 0; j < count; ++j) {
        mean += buffer.get(j);
    }
    mean /= count;

    double variance = 0.0;

    for (size_t j = 0; j < count; ++j) {
        variance += (buffer.get(j) - mean) * (buffer.get(j) - mean);
    }
    variance /= (count - 1);

    return variance;
}

class OnlineVariance {
public:
    static constexpr size_t BUFFER_SIZE = 10000; 

    OnlineVariance(size_t max_size, int window_size) 
        : max_size(max_size), window_size(window_size), welford(), buffer(window_size) {
        result_count = 0;
        buffered_results_welford = new double[BUFFER_SIZE]; 
        buffered_results_two_pass = new double[BUFFER_SIZE]; 
        buffered_count = 0; 
    }

    ~OnlineVariance() {
        delete[] buffered_results_welford;
        delete[] buffered_results_two_pass;
    }

    void process(const double* new_values, size_t batch_size) {
        welford.add_batch(new_values, batch_size);
        buffer.add_batch(new_values, batch_size);

        for (size_t i = 0; i < batch_size; ++i) {
            if (result_count < max_size) { 
                buffered_results_welford[buffered_count] = welford.get_variance();
                buffered_results_two_pass[buffered_count] = compute_two_pass_variance(buffer);
                buffered_count++;
                result_count++;
            }
        }
    }

    size_t get_result_count() const {
        return result_count;
    }

    CircularBuffer& get_buffer() {
        return buffer; 
    }

    double* get_buffered_results_welford() {
        return buffered_results_welford; 
    }

    double* get_buffered_results_two_pass() {
        return buffered_results_two_pass; 
    }

    size_t get_buffered_count() const {
        return buffered_count; 
    }

    void reset_buffered_count() {
        buffered_count = 0; 
    }

private:
    Welford welford;
    CircularBuffer buffer; 
    size_t max_size;
    size_t window_size;
    size_t result_count;
    double* buffered_results_welford; 
    double* buffered_results_two_pass; 
    size_t buffered_count; 
};

void write_results_to_file(const std::string& filename, const double* results, size_t result_count) {
    std::ofstream file(filename, std::ios::app);
    if (file.is_open()) {
        for (size_t i = 0; i < result_count; ++i) {
            file << results[i] << "\n";
        }
        file.close();
    } else {
        throw std::runtime_error("Could not open output file");
    }
}

class SignalsFromFile {
public:
    SignalsFromFile(const std::string& filename, size_t channel_count)
        : m_fscan(filename, std::ios::binary), CHANNELS(channel_count) {
        if (!m_fscan.is_open()) {
            throw std::runtime_error("Could not open file");
        }
    }

    bool read_channel(std::vector<double>& channel_data, size_t channel_index, size_t quant) {
        size_t len = CHANNELS * quant * sizeof(double);
        std::vector<double> buffer(CHANNELS * quant);

        if (m_fscan.eof()) return false;

        m_fscan.read(reinterpret_cast<char*>(buffer.data()), len);
        if (m_fscan.gcount() == 0) return false; 

        for (size_t i = 0; i < quant; ++i) {
            channel_data.push_back(buffer[i * CHANNELS + channel_index]);
        }
        return true;
    }

private:
    std::ifstream m_fscan;
    const size_t CHANNELS;
};

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <method: welford|two_pass>" << std::endl;
        return 1;
    }

    bool use_welford = (std::string(argv[1]) == "welford");
    const size_t data_size = 100000; 
    const size_t batch_size = 10000; 
    const size_t window_size = 500; 
    const size_t channel_index = 0; 

    OnlineVariance online_variance(data_size, window_size); 
    SignalsFromFile signals("C:/Users/Timofey/Downloads/8s001456.bin", 12);

    std::vector<double> channel_data;
    size_t quant = 40; 

    while (signals.read_channel(channel_data, channel_index, quant)) {
        online_variance.process(channel_data.data() + (channel_data.size() - quant), quant);
        if (online_variance.get_buffered_count() >= OnlineVariance::BUFFER_SIZE) {
            if (use_welford) {
                write_results_to_file("welford_results.txt", online_variance.get_buffered_results_welford(), online_variance.get_buffered_count());
            } else {
                write_results_to_file("two_pass_results.txt", online_variance.get_buffered_results_two_pass(), online_variance.get_buffered_count());
        }
    online_variance.reset_buffered_count(); 
        }

    }


    return 0;
}