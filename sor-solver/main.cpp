#include <iostream>
#include <cmath>
#include <utility>
#include <vector>
#include <cassert>

class Matrix {

    using arr_type = std::vector<double>;

public:
    [[maybe_unused]]
    Matrix(unsigned n, std::initializer_list<double> list)
            : size_y_(n),
              size_x_(n) {
        data_.reserve(size());
        data_.assign(list);
        while (data_.size() < size()) {
            data_.push_back(0);
        }
    }

    Matrix(unsigned cols, unsigned rows, std::initializer_list<double> list)
            : size_y_(rows),
              size_x_(cols) {
        data_.reserve(size());
        data_.assign(list);
        while (data_.size() < size()) {
            data_.push_back(0);
        }
    }

    double &operator()(unsigned x, unsigned y) {
        assert(x >= 0);
        assert(x < size_x_);
        assert(y >= 0);
        assert(y < size_y_);
        return data_.at(y * size_x_ + x);
    }

    double operator()(unsigned x, unsigned y) const {
        assert(x >= 0);
        assert(x < size_x_);
        assert(y >= 0);
        assert(y < size_y_);
        return data_.at(y * size_x_ + x);
    }

    double &operator()(unsigned i) {
        return data_.at(i);
    }

    double operator()(unsigned i) const {
        return data_.at(i);
    }

    void print() {
        unsigned i = -1;
        for (auto element : data_) {
            if (++i == size_x_) {
                std::cout << "\n";
                i = 0;
            }
            std::cout << element << "\t";
        }
    }

    unsigned x() const { return size_x_; }

    unsigned y() const { return size_y_; }

    unsigned c() const { return size_x_; }

    unsigned r() const { return size_y_; }

    unsigned size() const { return size_y_ * size_x_; }

private:
    unsigned size_x_;
    unsigned size_y_;

    arr_type data_;
};

double norm(const Matrix &input) {
    double sqr_sum = 0;
    for (int i = 0 ; i < input.size(); i++){
        sqr_sum += input(i) * input(i);
    }
    return std::sqrt(sqr_sum);
}

Matrix mul(const Matrix &A, const Matrix &B) {
    Matrix result(B.c(), A.r(), {});
    for (unsigned i = 0; i < A.r(); ++i)
        for (unsigned j = 0; j < B.c(); ++j)
            for (unsigned k = 0; k < A.c(); ++k) {
                result(j, i) += A(k, i) * B(j, k);
            }
    return result;
}

Matrix sub(const Matrix &A, const Matrix &B) {
    Matrix result(A.x(), A.y(), {});
    for (unsigned i = 0; i < A.x() * A.y(); i++) {
        result(i) = A(i) - B(i);
    }
    return result;
}


typedef Matrix Vector;

Vector sor_solver(Matrix A, Vector b,
                  double omega,
                  Vector initial_guess,
                  double max_diff) {

    unsigned step = 0;
    Vector phi = std::move(initial_guess);
    auto residual = norm(sub(mul(A, phi), b));
    while (residual > max_diff) {
        for (int i = 0; i < A.r(); i++) {
            double sigma = 0;
            for (int j = 0; j < A.c(); j++) {
                if (j != i) {
                    sigma += A(i, j) * phi(j);
                }
            }
            phi(i) = (1 - omega) * phi(i) + (omega / A(i, i)) * (b(i) - sigma);
        }
        residual = norm(sub(mul(A, phi), b));
        step++;
        printf("Step %d Residual: %lf\n", step, residual);
    }
    return phi;
}


int main() {
//    Matrix A(3, 2, {1, 2, 3, 4, 5, 6});
    Matrix B(2, 3, {7, 8, 9, 10, 11, 12});
    double residual_convergence = 1e-8;
    double omega = 0.5;  // Relaxation factor

    Matrix A(4, 4, {4, -1, -6, 0,
                    -5, -4, 10, 8,
                    0, 9, 4, -2,
                    1, 0, -7, 5});

    Vector b(1,4,{2, 21, -12, -6});

    A.print();
    std::cout << "\n";
    b.print();


    Vector initial_guess(1,4
                         ,{});

    Vector phi = sor_solver(A, b, omega, initial_guess, residual_convergence);
    std::cout << "Phi: " << std::endl;
    phi.print();


    return 0;
}
