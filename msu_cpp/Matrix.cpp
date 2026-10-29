#include "Matrix.h"
#include <cmath>
#include <stdexcept>
#include <algorithm> // для std::swap

void Matrix::alloc() {
    // выделяем один непрерывный кусок памяти размером n*n
    data = new double[n * n](); 
}

void Matrix::free() {
    delete[] data;
    data = nullptr; // нулевой указатель
}

Matrix::Matrix() : n(0), data(nullptr) {}

Matrix::Matrix(int n) : n(n) { alloc(); }

Matrix::Matrix(int n, double val) : n(n) {
    alloc();
    for (int i = 0; i < n * n; ++i)
        data[i] = val;
}

Matrix::Matrix(const Matrix& other) : n(other.n) {
    alloc();
    for (int i = 0; i < n * n; ++i)
        data[i] = other.data[i];
}
Matrix::~Matrix() { free(); }

Matrix& Matrix::operator=(const Matrix& other) {
    if (this != &other) {
        free();
        n = other.n;
        alloc();
        for (int i = 0; i < n * n; ++i)
            data[i] = other.data[i];
    }
    return *this;
}

int Matrix::size() const {
    return n;
}

double* Matrix::operator[](int i) { 
    return data + i * n; 
}

const double* Matrix::operator[](int i) const { 
    return data + i * n; 
}

double& Matrix::operator()(int i, int j) { 
    return data[i * n + j]; 
}

double Matrix::operator()(int i, int j) const { 
    return data[i * n + j]; 
}

// арифметика
Matrix Matrix::operator+(const Matrix& other) const {
    if (n != other.n) throw std::invalid_argument("size mismatch");
    Matrix res(n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            res.data[i * n + j] = data[i * n + j] + other.data[i * n + j]; 
    return res;
}

Matrix& Matrix::operator+=(const Matrix& other) {
    *this = *this + other;
    return *this;
}

Matrix Matrix::operator-(const Matrix& other) const {
    if (n != other.n) throw std::invalid_argument("size mismatch");
    Matrix res(n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            res.data[i * n + j] = data[i * n + j] - other.data[i * n + j];
    return res;
}

Matrix& Matrix::operator-=(const Matrix& other) {
    *this = *this - other;
    return *this;
}

Matrix Matrix::operator*(const Matrix& other) const {
    if (n != other.n) throw std::invalid_argument("size mismatch");
    Matrix res(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double sum = 0;
            for (int k = 0; k < n; ++k) {
                sum += data[i * n + k] * other.data[k * n + j];
            }
            res.data[i * n + j] = sum;
        }
    }
    return res;
}

Matrix& Matrix::operator*=(const Matrix& other) {
    *this = *this * other;
    return *this;
}

Matrix Matrix::operator*(double scalar) const {
    Matrix res(n);
    for (int i = 0; i < n * n; ++i)
        res.data[i] = data[i] * scalar;
    return res;
}

Matrix& Matrix::operator*=(double scalar) {
    for (int i = 0; i < n * n; ++i)
        data[i] *= scalar;
    return *this;
}

// сравнение
bool Matrix::operator==(const Matrix& other) const {
    if (n != other.n) return false;
    for (int i = 0; i < n * n; ++i)
        if (std::abs(data[i] - other.data[i]) > 1e-9)
            return false;
    return true;
}

bool Matrix::operator!=(const Matrix& other) const { 
    return !(*this == other); 
}

// ввод вывод
std::ostream& operator<<(std::ostream& os, const Matrix& m) {
    for (int i = 0; i < m.n; ++i) {
        for (int j = 0; j < m.n; ++j)
            os << m(i, j) << " "; // Используем оператор (), чтобы не писать m.data[i * m.n + j]
        os << "\n";
    }
    return os;
}

std::istream& operator>>(std::istream& is, Matrix& m) {
    for (int i = 0; i < m.n; ++i)
        for (int j = 0; j < m.n; ++j)
            is >> m(i, j);
    return is;
}

// транспонирование
Matrix Matrix::transpose() const {
    Matrix res(n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            res.data[j * n + i] = data[i * n + j]; // Обрати внимание на индексы!
    return res;
}

// единичная
Matrix Matrix::identity(int n) {
    Matrix I(n);
    for (int i = 0; i < n; ++i)
        I.data[i * n + i] = 1.0;
    return I;
}

// определитель методом Гаусса
double Matrix::determinant() const {
    Matrix tmp(*this);
    double det = 1.0;
    for (int i = 0; i < n; ++i) {
        int pivot = i;
        for (int r = i + 1; r < n; ++r)
            if (std::abs(tmp.data[r * n + i]) > std::abs(tmp.data[pivot * n + i]))
                pivot = r;
        if (std::abs(tmp.data[pivot * n + i]) < 1e-12) return 0.0;
        
        if (pivot != i) {
            // меняем строки местами (т.е. меняем элементы в цикле)
            for (int c = 0; c < n; ++c) {
                std::swap(tmp.data[i * n + c], tmp.data[pivot * n + c]);
            }
            det *= -1;
        }
        
        det *= tmp.data[i * n + i];
        for (int r = i + 1; r < n; ++r) {
            double factor = tmp.data[r * n + i] / tmp.data[i * n + i];
            for (int c = i; c < n; ++c)
                tmp.data[r * n + c] -= factor * tmp.data[i * n + c];
        }
    }
    return det;
}

// обратная (приписываем единичную и т.д.)
Matrix Matrix::inverse() const {
    Matrix A(*this);
    Matrix I = Matrix::identity(n);
    for (int i = 0; i < n; ++i) {
        int pivot = i;
        for (int r = i + 1; r < n; ++r)
            if (std::abs(A.data[r * n + i]) > std::abs(A.data[pivot * n + i]))
                pivot = r;
        if (std::abs(A.data[pivot * n + i]) < 1e-12)
            throw std::runtime_error("matrix is singular");
            
        if (pivot != i) {
            for (int c = 0; c < n; ++c) {
                std::swap(A.data[i * n + c], A.data[pivot * n + c]);
                std::swap(I.data[i * n + c], I.data[pivot * n + c]);
            }
        }
        
        double div = A.data[i * n + i];
        for (int c = 0; c < n; ++c) {
            A.data[i * n + c] /= div;
            I.data[i * n + c] /= div;
        }
        
        for (int r = 0; r < n; ++r) {
            if (r == i) continue;
            double factor = A.data[r * n + i];
            for (int c = 0; c < n; ++c) {
                A.data[r * n + c] -= factor * A.data[i * n + c];
                I.data[r * n + c] -= factor * I.data[i * n + c];
            }
        }
    }
    return I;
}

// ступенчатый вид
void Matrix::toRowEchelon() {
    for (int i = 0; i < n; ++i) {
        int pivot = i;
        for (int r = i + 1; r < n; ++r)
            if (std::abs(data[r * n + i]) > std::abs(data[pivot * n + i]))
                pivot = r;
        if (std::abs(data[pivot * n + i]) < 1e-12) continue;
        
        if (pivot != i) {
            for (int c = 0; c < n; ++c)
                std::swap(data[i * n + c], data[pivot * n + c]);
        }
        
        for (int r = i + 1; r < n; ++r) {
            double factor = data[r * n + i] / data[i * n + i];
            for (int c = i; c < n; ++c)
                data[r * n + c] -= factor * data[i * n + c];
        }
    }
}

// умножение на вектор
Matrix Matrix::multiplyVector(const double* vec) const {
    Matrix out(n);
    for (int i = 0; i < n; ++i) {
        double sum = 0;
        for (int j = 0; j < n; ++j)
            sum += data[i * n + j] * vec[j];
        out.data[i * n + 0] = sum; // результат в первый столбец
    }
    return out;
}
