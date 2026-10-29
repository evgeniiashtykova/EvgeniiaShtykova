#ifndef MATRIX_H // чтобы класс не создавался дважды
#define MATRIX_H

#include <iostream>

class Matrix {
private:
    int n;
    double* data;

    void alloc();   // выделяет память под массив
    void free();    // очищает

public:
    // конструкторы
    Matrix();                       
    Matrix(int n);         
    Matrix(int n, double val);      
    Matrix(const Matrix& other);    // копирование
    ~Matrix();                      

    // присваивание
    Matrix& operator=(const Matrix& other);

    Matrix operator+(const Matrix& other) const;
    Matrix& operator+=(const Matrix& other);
    Matrix operator-(const Matrix& other) const;
    Matrix& operator-=(const Matrix& other);
    Matrix operator*(const Matrix& other) const;  
    Matrix& operator*=(const Matrix& other);
    Matrix operator*(double scalar) const;        
    Matrix& operator*=(double scalar);

    // сравниваем
    bool operator==(const Matrix& other) const;
    bool operator!=(const Matrix& other) const;

    double* operator[](int i);    // можем взять строку с индексом i         
    const double* operator[](int i) const; 
    double& operator()(int i, int j);      // A(i,j)
    double operator()(int i, int j) const;

    // ввод и вывод
    friend std::ostream& operator<<(std::ostream& os, const Matrix& m);
    friend std::istream& operator>>(std::istream& is, Matrix& m);

    // операции
    Matrix transpose() const;              
    double determinant() const;            
    Matrix inverse() const;                
    static Matrix identity(int n);  // операция принадлежит всему классу       
    void toRowEchelon();     // к ступенчатому виду              
    Matrix multiplyVector(const double* vec) const; 

    int size() const;      // возвращает размер матрицы
};

#endif
