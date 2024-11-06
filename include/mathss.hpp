// !V1

#ifndef MATHSS_H
#define MATHSS_H

#include <iostream> // needed for ostream compatibility
#include <cmath> // needed for abs, sqrt, trigonometric functions
#include <initializer_list> // needed for vector constructors

// Swaps the content of two variables by performing xor operation
template <typename T> void swap(T* a, T* b);

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

    // simplifies the fraction to the lowest terms (in case of big numbers the fraction is approximated)
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

    // Converts given fraction to floating point value; approximation is given by double memory dimension
    explicit operator double() const&;

    // Returns the numerator of this fraction
    int getNumerator();

    // Returns the denominator of this fraction
    int getDenominator();

    // Returns the sign of this fraction
    int getSign();

    inline bool operator==(const Fraction& rhs);
    inline bool operator!=(const Fraction& rhs);
    inline bool operator< (const Fraction& rhs);
    inline bool operator> (const Fraction& rhs);
    inline bool operator<=(const Fraction& rhs);
    inline bool operator>=(const Fraction& rhs);

    // assign operator for double number
    inline void operator=(double num);

    // mathematical operations between fractions and fraction - double

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

    // print operator overload
    friend std::ostream& operator <<(std::ostream &out, const Fraction& f);
};

class Complex{
    private:
    double a,b;

    public:

    // Initialize the complex number at 0 + 0i
    Complex();

    // Initialize the complex number at the given value a + bi 
    Complex(double a, double b);

    // Copy constructor
    Complex(Complex& c);

    // equality operator checks for match of both parameters a and b
    inline bool operator==(const Complex& rhs);
    inline bool operator!=(const Complex& rhs);

    // compares norm of the complex number to verify the inequality

    inline bool operator< (const Complex& rhs);
    inline bool operator> (const Complex& rhs);
    inline bool operator<=(const Complex& rhs);
    inline bool operator>=(const Complex& rhs);

    // mathematical operations between Complex numbers and Complex 
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

    // get the module (sqrt(norm)) of a complex number
    static double mod(const Complex& c);

    // get the norm of a complex number, calculated as the sum of squares of a and b
    static double norm(const Complex& c);

    // get the real part of the given complex number
    static double Re(const Complex& c);

    // get the imaginary part of the given complex number
    static double Im(const Complex& c);

    // get the trigonometric argument of the given complex number
    static double Arg(const Complex& c);

    // returns the conjugate (a - bi) of the given complex number
    static Complex conjugate(const Complex& c);

    // returns the inverse (1 / (a + bi)) of the given complex number
    static Complex inverse(const Complex& c);

    // calculates the complex power of a given number to an integer exponent
    static Complex cpow(const Complex& c, const int exp);

    // calculates the n roots of a given number using the trigonometric form calculation 
    static Complex* croot(const Complex& c, const unsigned int exp);

    // print operator overload
    friend std::ostream& operator <<(std::ostream& out, const Complex& c);
};

template <typename T>
class Matrix{
    private:
    T* data;
    size_t rows, columns;

    public:

    // initialize a 0 x 0 empty matrix
    Matrix();

    // initialize an n x m matrix to a set value (default is 0)
    Matrix(size_t n, size_t m, T val = 0);

    // initialize a matrix from a nested initializer list; format like ({{1,0,0},{0,1,0},{0,0,1}})
    Matrix(const std::initializer_list<std::initializer_list<T>>& matrix);

    // copy constructor
    Matrix(const Matrix<T>& m);

    ~Matrix();

    // copy assignment operator
    void operator=(const Matrix<T>& m);

    // nested initializer list assignment operator
    void operator=(const std::initializer_list<std::initializer_list<T>>& matrix);

    // move assignment operator
    void operator=(Matrix<T>&& m);

    static Matrix<T> identity(size_t n);

    // returns a row from the matrix
    T* operator[](size_t index);

    // returns the number of rows of this matrix
    size_t getRows();

    // returns the number of columns of this matrix
    size_t getColumns();

    // first gaussian operation: swaps the elements of two rows
    void swapRows(size_t row1, size_t row2);

    // second gaussian operation: multiplies a row by a scalar value
    void multiplyRow(size_t row, T val);

    // third gaussian operation: adds the element of row1 (multiplied by a certain scalar value) to row2
    void addRow(size_t row1, size_t row2, T val = 1);

    // mathematical operations between matrices
    Matrix<T> operator+(const Matrix<T>& m);
    Matrix<T>& operator+=(const Matrix<T>& m);
    Matrix<T> operator-(const Matrix<T>& m);
    Matrix<T>& operator-=(const Matrix<T>& m);
    Matrix<T> operator*(const Matrix<T>& m);

    // matrix * scalar value multiplication
    template<typename U> friend Matrix<U> operator*(const U& val, const Matrix<U>& m);
    Matrix<T> operator*(const T& val);
    Matrix<T>& operator*=(const T& val);

    // returns the transposed matrix of the given one
    static Matrix<T> transpose(const Matrix<T>& m);

    // computes the determinant of the given matrix; uses computational formulas until 3 dimension, gaussian elimination is used for bigger matrices
    static T det(const Matrix<T>& m);

    // print operator overload
    template<typename U> friend std::ostream& operator <<(std::ostream &out, const Matrix<U> &v);
};

// lpu decomposition: given a matrix M, the function returns three matrices L, P, U such that PM = LU
template<typename T>
Matrix<T>* LPU_decomposition(const Matrix<T> &m);

// qr decomposition: given a matrix M, the function returns an orthogonal matrix Q and an upper triangular matrix R such that M = QR 
template<typename T>
Matrix<T>* QR_decomposition(const Matrix<T> &m);

template <typename T>
class Vec {
    private:
    T *val;
    size_t dim;

    public:

    // initialize a 0 dimension vector
    Vec();

    // initialize an n-dimensional vector to a given value(default is 0)
    Vec(const size_t& n, const T& val = 0);

    // initialize an n-dimensional vector from values stored in an array
    Vec(const T arr[],const size_t& n);

    // initialize a vector from an initializer list; format like: ({num-1, num-2, ... , num-n})
    Vec(const std::initializer_list<T>& arr);

    // copy constructor
    Vec(const Vec<T>& v);

    ~Vec();

    // initializer list assignment operator
    void operator=(const std::initializer_list<T> & arr);

    // copy assignment operator
    void operator=(const Vec<T>& v);

    // move operator
    void operator=(Vec<T>&& v);

    // equality operator check for a match of every element of the vectors
    inline bool operator==(const Vec<T>& rhs);
    inline bool operator!=(const Vec<T>& rhs);

    // confrontation operators use norm to solve the inequalities
    inline bool operator< (const Vec<T>& rhs);
    inline bool operator> (const Vec<T>& rhs);
    inline bool operator<=(const Vec<T>& rhs);
    inline bool operator>=(const Vec<T>& rhs);

    // returns the vectors' dimension
    size_t getDim();

    // returns norm of the vector, calculated as the sum of squares of the single elements
    T getNorm();

    // get an element of the vector
    T& operator[](size_t index);

    explicit operator Matrix<T>() const&;

    // mathematical operation between vectors and multiplication with scalar value
    Vec<T> operator+(const Vec<T>& vec);
    Vec<T>& operator+=(const Vec<T>& vec);
    Vec<T> operator-(const Vec<T>& vec);
    Vec<T>& operator-=(const Vec<T>& vec);
    Vec<T> operator*(const T& val);
    
    template<typename U> friend Vec<U> operator*(const U& val, const Vec<U>& vec);
    Vec<T>& operator*=(const T& val);

    // returns the dot product between vectors of the same dimension 
    static T dot(const Vec<T>& lhs,const Vec<T>& rhs);

    // calculates the cross product of 3 dimensional vectors
    static Vec<T> cross(const Vec<T>& lhs,const Vec<T>& rhs);

    // print operator overload
    template<typename U> friend std::ostream& operator <<(std::ostream &out,const Vec<U> &v);
};

// include template vector class implementation
#include "Vec.ipp"

// include template matrix class implementation
#include "Matrix.ipp"

#endif