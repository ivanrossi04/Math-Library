// !V1

#ifndef MATHSS_H
#define MATHSS_H

#include <iostream> // needed for ostream compatibility
#include <cmath> // needed for abs
#include <initializer_list> // needed for vector constructors

// Swaps the content of two variables by performing xor operation
template <typename T> void swap(T *a, T *b);

// Finds the Greatest Common Divisor between two whole numbers
int GCD(unsigned long long a, unsigned long long b);

// Finds the Lowest Common Multiple between two whole numbers
int LCM(unsigned long long a, unsigned long long b);

class Fraction {
private:
    /* data */
    unsigned long numerator;
    unsigned long denominator;
    bool sign; // 0 = positive, 1 = negative

    static void simplify(Fraction *frac);
    
public:
    // Constructors

    // Initialize the fraction at 0
    Fraction();

    // Initialize the fraction at set values
    Fraction(long long numerator, long long denominator);

    // Approximates the given number to a fraction; precision of set value is defined by default at 10e7
    Fraction(double num, int precision = 10e5);
    
    // Destructor
    ~Fraction();

    // Converts given fraction to floating point value; approximation is given by float memory dimension
    explicit operator float() const&;
    explicit operator double() const&;

    int getNumerator();
    int getDenominator();
    int getSign();

    inline bool operator==(const Fraction& rhs);
    inline bool operator!=(const Fraction& rhs);
    inline bool operator< (const Fraction& rhs);
    inline bool operator> (const Fraction& rhs);
    inline bool operator<=(const Fraction& rhs);
    inline bool operator>=(const Fraction& rhs);

    inline void operator=(double num);

    Fraction operator+(const Fraction& frac);
    Fraction& operator+=(const Fraction& frac);
    Fraction operator+(const double& num);
    friend Fraction operator+(const double& num, const Fraction& frac);
    Fraction& operator+=(const double& num);

    Fraction& operator++();
    Fraction operator++(int);


    Fraction operator-(const Fraction& frac);
    Fraction& operator-=(const Fraction& frac);
    Fraction operator-(const double& num);
    friend Fraction operator-(const double& num, const Fraction& frac);
    Fraction& operator-=(const double& num);

    Fraction& operator--();
    Fraction operator--(int);


    Fraction operator*(const Fraction& frac);
    Fraction& operator*=(const Fraction& frac);
    Fraction operator*(const double& num);
    friend Fraction operator*(const double& num, const Fraction& frac);
    Fraction& operator*=(const double& num);


    Fraction operator /(const Fraction& frac);
    Fraction& operator/=(const Fraction& frac);
    Fraction operator /(const double& num);
    friend Fraction operator/(const double& num, const Fraction& frac);
    Fraction& operator/=(const double& num);

    friend std::ostream& operator <<(std::ostream &out, const Fraction& f);
};

class Complex{
    private:
    double a,b;

    public:
    Complex();
    Complex(double a, double b);
    Complex(Complex& c);

    inline bool operator==(const Complex& rhs);
    inline bool operator!=(const Complex& rhs);
    inline bool operator< (const Complex& rhs);
    inline bool operator> (const Complex& rhs);
    inline bool operator<=(const Complex& rhs);
    inline bool operator>=(const Complex& rhs);

    Complex operator+(const Complex& c);
    Complex& operator+=(const Complex& c);

    Complex operator-(const Complex& c);
    Complex& operator-=(const Complex& c);

    Complex operator*(const Complex& c);
    Complex& operator*=(const Complex& c);

    Complex operator/(const Complex& c);
    Complex& operator/=(const Complex& c);

    Complex operator*(const double& n);
    friend Complex operator*(const double& n, const Complex& c);
    Complex& operator*=(const double& n);

    Complex operator/(const double& n);
    friend Complex operator/(const double& n, const Complex& c);
    Complex& operator/=(const double& n);

    static double mod(const Complex& c);
    static double norm(const Complex& c);
    static double Re(const Complex& c);
    static double Im(const Complex& c);
    static double Arg(const Complex& c);

    static Complex conjugate(const Complex& c);
    static Complex inverse(const Complex& c);
    static Complex cpow(const Complex& c, const int exp);
    static Complex* croot(const Complex& c, const unsigned int exp);

    friend std::ostream& operator <<(std::ostream& out, const Complex& c);
};

template <typename T>
class Vec {
    private:
    T *val;
    size_t dim;

    public:
    Vec();
    Vec(const size_t n);
    Vec(const T arr[],const size_t size);
    Vec(const std::initializer_list<T>& arr);
    Vec(const Vec<T>& v);

    ~Vec();

    void operator=(const std::initializer_list<T> & arr);
    void operator=(const Vec<T>& v);
    void operator=(Vec<T>&& v);

    // confrontation operator working with norm
    inline bool operator==(const Vec<T>& rhs);
    inline bool operator!=(const Vec<T>& rhs);
    inline bool operator< (const Vec<T>& rhs);
    inline bool operator> (const Vec<T>& rhs);
    inline bool operator<=(const Vec<T>& rhs);
    inline bool operator>=(const Vec<T>& rhs);

    size_t getDim();
    T getNorm();
    T& operator[](size_t index);

    Vec<T> operator+(const Vec<T>& vec);
    Vec<T>& operator+=(const Vec<T>& vec);
    Vec<T> operator-(const Vec<T>& vec);
    Vec<T>& operator-=(const Vec<T>& vec);
    Vec<T> operator*(const T& val);
    
    template<typename U> friend Vec<U> operator*(const U& val, const Vec<U>& vec);
    Vec<T>& operator*=(const T& val);

    static T dot(const Vec<T>& lhs,const Vec<T>& rhs);
    static Vec<T> cross(const Vec<T>& lhs,const Vec<T>& rhs);

    template<typename U> friend std::ostream& operator <<(std::ostream &out,const Vec<U> &v);
};

// include template vector class implementation
#include "Vec.ipp"

template <typename T>
class Matrix{
    private:
    T** data;
    size_t rows, columns;

    public:
    Matrix();
    Matrix(size_t n, size_t m, T val = 0);
    Matrix(const std::initializer_list<std::initializer_list<T>>& matrix);
    Matrix(const Matrix<T>& m);

    ~Matrix();

    void operator=(const Matrix<T>& m);
    void operator=(const std::initializer_list<std::initializer_list<T>>& matrix);
    void operator=(Matrix<T>&& m);

    static Matrix<T> identity(size_t n);

    T* operator[](size_t index);
    size_t getRows();
    size_t getColumns();

    void swapRows(size_t row1, size_t row2);
    void multiplyRow(size_t row, T val);
    void addRow(size_t row1, size_t row2, T val);

    Matrix<T> operator+(const Matrix<T>& m);
    Matrix<T>& operator+=(const Matrix<T>& m);
    Matrix<T> operator-(const Matrix<T>& m);
    Matrix<T>& operator-=(const Matrix<T>& m);
    Matrix<T> operator*(const Matrix<T>& m);

    template<typename U> friend Matrix<U> operator*(const U& val, const Matrix<U>& m);
    Matrix<T> operator*(const T& val);
    Matrix<T>& operator*=(const T& val);

    static Matrix<T> transpose(const Matrix<T>& m);
    static T det(const Matrix<T>& m);

    // static Matrix<T> inv(const Matrix<T>& m);

    template<typename U> friend std::ostream& operator <<(std::ostream &out, const Matrix<U> &v);
};

// linear algebra
template<typename T>
Matrix<T>* LPU_decomposition(const Matrix<T> &m);

template<typename T>
Matrix<T>* QR_decomposition(const Matrix<T> &m);
// TODO: inverse matrix

// diagonalization of n*n matrices
// linear system resolution

// include template matrix class implementation
#include "Matrix.ipp"

#endif