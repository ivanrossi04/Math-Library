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
    this -> numerator = llabs(numerator);
    this -> denominator = llabs(denominator);
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

        // number precision cutoff: if a fraction is too precise it gets approximated
        while(frac -> denominator >= 10 && frac -> numerator * frac -> denominator > 2E14){
            frac -> numerator /= 10 + (frac -> numerator % 10 >=5);
            frac -> denominator /= 10;
        }

        // find the greatest common divisor between numerator & denominator to simplify the fraction
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
    // returns the division between numerator and denominator
    return ((float) numerator / denominator * (sign ? -1 : 1));
}

Fraction::operator double() const&{
    // returns the division between numerator and denominator
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
    // the match is done between the fraction converted to floats
    return (float) (*this) < (float)(rhs);
}

inline bool Fraction::operator>(const Fraction& rhs){
    // the match is done between the fraction converted to floats
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

    // result sign initialization
    result.sign = this -> sign;

    // denominator is set to the Lowest Common Denominator
    result.denominator = LCM(this -> denominator, frac.denominator);

    // the addends are calculated multiplying by the fraction with the LCD
    unsigned int add1 = this -> numerator * (result.denominator / this -> denominator);
    unsigned int add2 = frac.numerator * (result.denominator /  frac.denominator);

    // if the fraction have the same sign simply add the numerators
    if(this -> sign == frac.sign){
        result.numerator = add1 + add2;
    } else {
        // a subtraction between unsigned integers is performed

        // need to check which number is bigger to avoid underflow
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

    // result sign initialization
    result.sign = this -> sign;

    // denominator is set to the Lowest Common Denominator
    result.denominator = LCM(this -> denominator, frac.denominator);

    // the operands are calculated multiplying by the fraction with the LCD
    unsigned int add1 = this -> numerator * (result.denominator / this -> denominator);
    unsigned int add2 = frac.numerator * (result.denominator /  frac.denominator);

    if(this -> sign == frac.sign){
        // the subtraction is performed

        if(add1 < add2){
            swap(&add1, &add2);
            result.sign = !(this -> sign);
        }
        result.numerator = add1 - add2;
    } else {
        // if signs are different the operation becomes a sum
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

    // cross simplifications are calculated to keep the multiplication at lowest terms possible
    int cross_simplify1 = GCD(this -> numerator, frac.denominator);
    int cross_simplify2 = GCD(this -> denominator, frac.numerator);

    Fraction result;

    // the multiplication is performed, simplifying before the actual operation
    result.numerator = (this -> numerator / cross_simplify1) * (frac.numerator / cross_simplify2);
    result.denominator = (this -> denominator / cross_simplify2) * (frac.denominator / cross_simplify1);
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
    // throwing possible exception
    if(frac.numerator == 0) throw "Invalid operation: can't divide by zero";

    // division is performed like a multiplication with the inverse fraction of the second operand

    int cross_simplify1 = GCD(this -> numerator, frac.numerator);
    int cross_simplify2 = GCD(this -> denominator, frac.denominator);

    Fraction result;
    result.numerator = this -> numerator / cross_simplify1 * frac.denominator / cross_simplify2;
    result.denominator = this -> denominator / cross_simplify2 * frac.numerator / cross_simplify1;
    result.sign = this -> sign ^ frac.sign;
    return result;
}

Fraction& Fraction::operator/=(const Fraction& frac){
    *this = *this / frac;
    return *this;
}

Fraction Fraction::operator/(const double& num){
    return operator/(Fraction(num));
}

Fraction operator/(const double& num, const Fraction& frac){
    return Fraction(num).operator/(frac);
}

Fraction& Fraction::operator/=(const double& num){
    return operator/=(Fraction(num));
}

std::ostream& operator <<(std::ostream& out, const Fraction& f) {

    // print example: std::cout << Fraction(-2, 1); -> -2/1
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
    Complex result = operator*(inverse(c));
    return result;
}

Complex& Complex::operator/=(const Complex& c){
    operator*=(inverse(c));
    return *(this);
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
    // the inverse is calculated as 1/z = 1/(a+bi) * (a-bi)/(a-bi) = (a-bi)/(a^2+b^2)
    return Complex::conjugate(c) / Complex::norm(c);
}

Complex Complex::cpow(const Complex& c, const int exp){
    // based on deMoivre's formula (cos(x) + i sin(x))^n = (cos(nx) + i sin(nx))
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
    // based on deMoivre's formula
    if(exp){
        double length = pow(mod(c), 1/exp);
        double theta = Complex::Arg(c);

        Complex *results = new Complex[exp];
        for(int i = 0; i < exp; i++){
            // arg(z_k) = (arg(Z) + 2kPI) / n
            double root_theta = (theta + 2 * i * M_PI) / exp;

            // z_k = mod(Z)(cos(arg(z_k)) + i sin(arg(z_k)))
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

    // print example: std::cout << Complex(-2, -1); -> -2-1i
    out << c.a << (c.b >= 0 ? "+" : "") << c.b << "i";

    return out;
}


// ------- General functions implementation -------

int GCD(unsigned long long a, unsigned long long b){
    // implementation of euclidean algorithm
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