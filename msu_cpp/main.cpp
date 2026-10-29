#include "Matrix.h"
#include <iostream>

int main() {
    int n = 3;
    Matrix A(n);
    std::cout << "Введите " << n*n << " чисел для матрицы A:\n";
    std::cin >> A;

    std::cout << "A =\n" << A;
    std::cout << "det(A) = " << A.determinant() << "\n";
    std::cout << "A^T =\n" << A.transpose();

    try {
        Matrix inv = A.inverse();
        std::cout << "A^{-1} =\n" << inv;
        std::cout << "A * A^{-1} =\n" << (A * inv);
    } catch (std::exception& e) {
        std::cout << "Ошибка - не удалось найти обратную матрицу \n";
    }

    Matrix B(n, 2.0);
    std::cout << "A + B =\n" << (A + B);
    return 0;
}
