// !V1
#include <iostream>
#include "mathss.hpp"

unsigned int malloc_count = 0;
void print_allocations(){
    std::cout << "Allocations: " << malloc_count << "\n";
}

void* operator new(size_t size){
    malloc_count++;
    void *p = malloc(size);
    return p;
}

int main(int argc, char *argv[]){
    /*
    int a = 15, b = 42;

    std::cout << a << " " << b << std::endl;
    swap(&a, &b);
    std::cout << a << " " << b << std::endl; 

    std::cout << GCD(a, b) << " " << LCM(a, b) << std::endl;


    std::cout << "\nFraction class - Test\n";
    Fraction f1 = Fraction();
    Fraction f2 = Fraction(2, 10);
    Fraction f3 = Fraction(-3.0);
    Fraction f4 = 4.5;
    Fraction f5 = f1;

    std::cout << f1 << " | " << f2 << " | " << f3 << " | " << f4 << " | " << f5 << "\n";

    f1 = f2 / f3;
    std::cout << f1 << "\n";

    f1 /= f2;
    std::cout << f1 << "\n";

    f1 /= 0.5;
    std::cout << f1 << "\n";

    std::cout << f3 << "\n";

    std::cout<< "\nComplex class - Test\n";

    Complex c1 = Complex();
    Complex c2 = Complex(3.4, 2);
    Complex c3 = Complex(c2);
    
    std::cout << c1 << " | " << c2 << " | " << c3 << "\n"; 

    c1 = c2 / c3;
    std::cout << c1 << "\n";

    c1 /= c2;
    std::cout << c1 << "\n";

    std::cout << Complex::mod(c1) << " | " << Complex::Arg(c1) << "\n";

    c1 = Complex::conjugate(c1);
    std::cout << c1 << "\n";

    c1 = Complex::inverse(c1);
    std::cout << c1 << "\n";
    */
    /*
    std::cout << "\nVector class - Test\n";

    float nums[3] = {4,5,1};
    Vec<float> v1 = Vec<float>();
    Vec<float> v2 = Vec<float>(nums, 3);
    Vec<float> v3 = Vec<float>({2,3,4});

    std::cout << v1 << " | " << v2 << " | " << v3 << "\n";
    print_allocations();

    v1 = Vec<float>({2,1});
    std::cout << v1 << "\n";
    print_allocations();
    
    v1 = v2 + v3;
    std::cout << v1 << "\n";
    print_allocations();
    
    v1 += v2;
    std::cout << v1 << "\n";
    print_allocations();
    
    v1 = v1 - v2;
    std::cout << v1 << "\n";
    print_allocations();
    
    v1 -= v2;
    std::cout << v1 << "\n";
    print_allocations();
    
    std::cout << Vec<float>::cross(v1, v2) << "\n";
    print_allocations();
    */

    std::cout << "Matrix class - Test\n";
    
    Matrix<double> m1 = Matrix<double>();
    Matrix<double> m2({
        {2,3,4},
        {2,3,6},
        {-4,3,3}
    });
    
    std::cout << m1 << "\n\n" << m2 << "\n\n";
    
    m1={{2,1,2},
        {4,-2,0},
        {1,0,4}};
    std::cout << m1 << "\n";

    m1 = m1 + m2;
    std::cout << m1 << "\n";

    m1 += m2;
    std::cout << m1 << "\n";

    m1 = m1 - m2;
    std::cout << m1 << "\n";

    m1 -= m2;
    std::cout << m1 << "\n";

    Matrix<double> m3({
        {2,3,1},
        {-1,0,0}
    });

    m1 = m3 * m1;
    std::cout << m1 << "\n";

    m3 *= 2;
    std::cout << Matrix<double>::transpose(m2) << "\n";  
    std::cout << "det(A) = " << Matrix<double>::det(m2) << "\n";

    Matrix<double> *result = LPU_decomposition(m2);
    std::cout << result[0] << "\n" << result[2];
    Matrix<double> LU = result[0] * result[2];
    std::cout << "L * U:\n" << LU << "\ndet(LU) = " << Matrix<double>::det(LU) << "\n";
    std::cout << "P * A:\n" << result[1] * m2;
    
    Matrix<double> *QR = QR_decomposition(m2);
    std::cout << QR[0] << "\n" << QR[1];
    std::cout << "\n" << QR[0] * QR[1];

    print_allocations();

    return 0;
}