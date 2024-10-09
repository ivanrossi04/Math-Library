// !V1
#ifndef MATHSS_C
#define MATHSS_C

#include "mathss.hpp"

// ------- Fraction methods implementation -------

Fraction::Fraction(){
    this -> numerator = 0;
    this -> denominator = 1;
    this -> sign = 0;
}

Fraction::Fraction(long long numerator, long long denominator){
    this -> numerator = abs(numerator);
    this -> denominator = abs(denominator);
    this -> sign = (numerator < 0) ^ (denominator < 0);
    
    simplify(this);
}

Fraction::Fraction(double num, int precision){
    if(num){
        this -> sign = num < 0;

        this -> numerator = num * precision * (num >= 0 ?: -1);
        this -> denominator = precision;

        simplify(this);        
    } else Fraction();
}

Fraction::~Fraction(){}

void Fraction::simplify(Fraction *frac){
    if(frac -> numerator == 0) frac -> denominator = 1;
    else {
        while(frac -> denominator >= 10 && frac -> numerator * frac -> denominator > 2E14){
            frac -> numerator /= 10;
            frac -> denominator /= 10;
        }

        int simplify = GCD(frac -> numerator, frac -> denominator);

        frac -> numerator /= simplify;
        frac -> denominator /= simplify;
    }
}

int Fraction::getNumerator(){
    return numerator;
}

int Fraction::getDenominator(){
    return denominator;
}

int Fraction::getSign(){
    return sign ? -1 : 1;
}

Fraction::operator float() const&{
    return ((float) numerator / denominator * (sign ? -1 : 1));
}

Fraction::operator double() const&{
    return ((double) numerator / denominator * (sign ? -1 : 1));
}

void Fraction::operator=(double num){
    *this = Fraction(num);
    return;
}

inline bool Fraction::operator==(const Fraction& rhs){
    return this -> numerator == rhs.numerator && this -> numerator == rhs.denominator;
}

inline bool Fraction::operator!=(const Fraction& rhs){
    return !(*this == rhs);
}

inline bool Fraction::operator<(const Fraction& rhs){
    return (float) (*this) < (float)(rhs);
}

inline bool Fraction::operator>(const Fraction& rhs){
    return (float) (*this) > (float) (rhs);
}

inline bool Fraction::operator<=(const Fraction& rhs){
    return !this -> operator>(rhs);
}

inline bool Fraction::operator>=(const Fraction& rhs){
    return !this -> operator<(rhs);
}

Fraction Fraction::operator+(const Fraction& frac){
    Fraction result;

    result.sign = this -> sign;

    result.denominator = LCM(this -> denominator, frac.denominator);

    unsigned int add1 = this -> numerator * (result.denominator / this -> denominator);
    unsigned int add2 = frac.numerator * (result.denominator /  frac.denominator);

    if(this -> sign == frac.sign){
        result.numerator = add1 + add2;
    } else {
        if(add1 < add2){
            swap(&add1, &add2);
            result.sign = !(this -> sign);
        }
        result.numerator = add1 - add2;
    }

    simplify(&result);
    return result;
}

Fraction& Fraction::operator+=(const Fraction& frac){
    *this = *this + frac;
    return *this;
}

Fraction Fraction::operator+(const double& num){
    return operator+(Fraction(num));
}

Fraction operator+(const double& num, const Fraction& frac){
    return Fraction(num).operator+(frac);
}

Fraction& Fraction::operator+=(const double& num){
    return operator+=(Fraction(num));
}

Fraction& Fraction::operator++(){
    return operator+=(1);
}

Fraction Fraction::operator++(int){
    Fraction old = *this;
    operator++();
    return old;
}

Fraction Fraction::operator-(const Fraction& frac){
    Fraction result;

    result.sign = this -> sign;
    result.denominator = LCM(this -> denominator, frac.denominator);

    unsigned int add1 = this -> numerator * (result.denominator / this -> denominator);
    unsigned int add2 = frac.numerator * (result.denominator /  frac.denominator);

    if(this -> sign == frac.sign){
        if(add1 < add2){
            swap(&add1, &add2);
            result.sign = !(this -> sign);
        }
        result.numerator = add1 - add2;
    } else {
        result.numerator = add1 + add2;
    }

    simplify(&result);
    return result;
}

Fraction& Fraction::operator-=(const Fraction& frac){
    *this = *this - frac;
    return *this;
}

Fraction Fraction::operator-(const double& num){
    return operator-(Fraction(num));
}

Fraction operator-(const double& num, const Fraction& frac){
    return Fraction(num).operator-(frac);
}

Fraction& Fraction::operator-=(const double& num){
    return operator-=(Fraction(num));
}

Fraction& Fraction::operator--(){
    return operator-=(1);
}

Fraction Fraction::operator--(int){
    Fraction old = *this;
    operator--();
    return old;
}

Fraction Fraction::operator*(const Fraction& frac){

    int cross_simplify1 = GCD(this -> numerator, frac.denominator);
    int cross_simplify2 = GCD(this -> denominator, frac.numerator);

    Fraction result;
    result.numerator = this -> numerator / cross_simplify1 * frac.numerator / cross_simplify2;
    result.denominator = this -> denominator / cross_simplify2 * frac.denominator / cross_simplify1;
    result.sign = this -> sign ^ frac.sign;
    return result;
}

Fraction& Fraction::operator*=(const Fraction& frac){
    *this = *this * frac;
    return *this;
}

Fraction Fraction::operator*(const double& num){
    return operator*(Fraction(num));
}

Fraction operator*(const double& num, const Fraction& frac){
    return Fraction(num).operator*(frac);
}

Fraction& Fraction::operator*=(double const& num){
    return operator*=(Fraction(num));
}

Fraction Fraction::operator/(const Fraction& frac){
    if(frac.numerator == 0) throw "Invalid operation: can't divide by zero";

    int cross_simplify1 = GCD(this -> numerator, frac.numerator);
    int cross_simplify2 = GCD(this -> denominator, frac.denominator);

    Fraction result;
    result.numerator = this -> numerator / cross_simplify1 * frac.denominator / cross_simplify2;
    result.denominator = this -> denominator / cross_simplify2 * frac.numerator / cross_simplify1;
    result.sign = this -> sign ^ frac.sign;
    return result;
}

Fraction& Fraction::operator/=(const Fraction& frac){
    try{
        *this = *this / frac;
    } catch(char* err){
        throw err;
    }
    return *this;
}

Fraction Fraction::operator/(const double& num){
    try{
        return operator/(Fraction(num));
    } catch(char* err){
        throw err;
        return *this;
    }
}

Fraction operator/(const double& num, const Fraction& frac){
    try{
        return Fraction(num).operator/(frac);
    } catch(char* err){
        throw err;
        return frac;
    }
}

Fraction& Fraction::operator/=(const double& num){
    try{
        return operator/=(Fraction(num));
    } catch(char* err){
        throw err;
        return *this;
    }
}

std::ostream& operator <<(std::ostream& out, const Fraction& f) {
    out << (f.sign ? "-" : "") << f.numerator << "/" << f.denominator;

    return out;
}

// ------- Complex methods implementation -------
Complex::Complex(){
    this -> a = 0;
    this -> b = 0;
}

Complex::Complex(double a, double b){
    this -> a = a;
    this -> b = b;
}

Complex::Complex(Complex& c){
    this -> a = c.a;
    this -> b = c.b;
}

bool Complex::operator==(const Complex& c){return this -> a == c.a && this -> b == c.b;}
bool Complex::operator!=(const Complex& c){return !operator==(c);}
bool Complex::operator<(const Complex& c){return Complex::norm(*(this)) < Complex::norm(c);}
bool Complex::operator>(const Complex& c){return Complex::norm(*(this)) > Complex::norm(c);}
bool Complex::operator<=(const Complex& c){return !operator>(c);}
bool Complex::operator>=(const Complex& c){return !operator<(c);}

Complex Complex::operator+(const Complex& c){
    Complex result(
        this -> a + c.a,
        this -> b + c.b
    );
    return result;
}

Complex& Complex::operator+=(const Complex& c){
    this -> a += c.a;
    this -> b += c.b;
    return *(this);
}

Complex Complex::operator-(const Complex& c){
    Complex result(
        this -> a - c.a,
        this -> b - c.b
    );
    return result;
}

Complex& Complex::operator-=(const Complex& c){
    this -> a -= c.a;
    this -> b -= c.b;
    return *(this);
}

Complex Complex::operator*(const Complex& c){
    Complex result(
        this -> a * c.a - this -> b * c.b,
        this -> a * c.b + this -> b * c.a
    );
    return result;
}

Complex& Complex::operator*=(const Complex& c){
    double temp = this -> a;
    this -> a = this -> a * c.a - this -> b * c.b;
    this -> b = temp * c.b + this -> b * c.a;
    return *(this);
}

Complex Complex::operator/(const Complex& c){
    Complex result = operator*(conjugate(c)) / norm(c);
    return result;
}

Complex& Complex::operator/=(const Complex& c){
    operator*=(conjugate(c));
    return operator/=(norm(c));
}

Complex Complex::operator*(const double& n){
    Complex result(
        this -> a * n,
        this -> b * n
    );
    return result;
}
Complex operator*(const double& n, const Complex& c){
    Complex result(
        c.a * n,
        c.b * n
    );
    return result;
}
Complex& Complex::operator*=(const double& n){
    this -> a *= n;
    this -> b *= n;
    return *(this);
}
    
Complex Complex::operator/(const double& n){
    Complex result(
        this -> a / n,
        this -> b / n
    );
    return result;
}

Complex operator/(const double& n, const Complex& c){
    return n * Complex::inverse(c);
}

Complex& Complex::operator/=(const double& n){
    this -> a /= n;
    this -> b /= n;
    return *(this);
}

double Complex::mod(const Complex& c){return sqrt(norm(c));}
double Complex::norm(const Complex& c){return c.a * c.a + c.b * c.b;}
double Complex::Re(const Complex& c){return c.a;}
double Complex::Im(const Complex& c){return c.b;}
double Complex::Arg(const Complex& c){return asin(c.b / mod(c)) * (c.a >= 0 ?: -1);}

Complex Complex::conjugate(const Complex& c){
    return Complex(
        c.a,
        -c.b
    );
}

Complex Complex::inverse(const Complex& c){
    return Complex::conjugate(c) / Complex::norm(c);
}

Complex Complex::cpow(const Complex& c, const int exp){
    if(exp){
        double length = pow(mod(c), exp);
        double theta = Complex::Arg(c);
        Complex result(
            cos(abs(exp) * theta),
            sin(abs(exp) * theta)
        );

        if(exp > 0) return length * result;
        return length / result;
    }
    return Complex(1, 0);
}

Complex* Complex::croot(const Complex& c,const unsigned int exp){
    if(exp){
        double length = pow(mod(c), 1/exp);
        double theta = Complex::Arg(c);

        Complex *results = new Complex[exp];
        for(int i = 0; i < exp; i++){
            double root_theta = (theta + 2 * i * M_PI) / exp;
            results[i] = length * Complex(
                cos(root_theta),
                sin(root_theta)
            );
        }

        return results;
    }
    return new Complex(1, 0);
}

std::ostream& operator <<(std::ostream& out,const Complex& c) {
    out << c.a << " + " << c.b << " i";

    return out;
}


// ------- General functions implementation -------

template <typename T>
void swap(T *a, T *b){
    *b = *a ^ *b;
    *a = *a ^ *b;
    *b = *a ^ *b;
}

int GCD(unsigned long long a, unsigned long long b){
    if(b == 0) return a;
    if(a == 0) return b;
    if(a < b) swap(&a, &b);
    int temp;
    do{
        temp = a;
        a = b;
        b = temp % b;
    }while(b != 0);
    return a;
}

int LCM(unsigned long long a, unsigned long long b){
    return a / GCD(a, b) * b;
}

#endif