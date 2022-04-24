#pragma once

#include <ostream>

namespace linalg {

////////////////////////////////////////////////////////////////////////
//                            DESCRIPTION                             //
////////////////////////////////////////////////////////////////////////

// An implementation of complex numbers.

////////////////////////////////////////////////////////////////////////
//                            DECLARATIONS                            //
////////////////////////////////////////////////////////////////////////

template <typename Real>
class Complex {
   public:
    Real real, imag;

    Complex() : real(0), imag(0){};
    Complex(Real real) : real(real), imag(0){};
    Complex(Real real, Real imag) : real(real), imag(imag){};
    Complex(const Complex &c) : real(c.real), imag(c.imag){};

    bool operator==(Complex c) const;
    bool operator!=(Complex c) const;

    Complex operator+(Complex c) const;
    Complex operator-(Complex c) const;
    Complex operator-() const;

    Complex operator*(Complex c) const;
    Complex operator/(Complex c) const;

    Complex conj() const;
    Real mag2() const;
};

template <typename Real>
std::ostream &operator<<(std::ostream &os, Complex<Real> c);

////////////////////////////////////////////////////////////////////////
//                            DEFINITIONS                             //
////////////////////////////////////////////////////////////////////////

template <typename Real>
bool Complex<Real>::operator==(Complex<Real> c) const {
    return real == c.real && imag == c.imag;
}

template <typename Real>
bool Complex<Real>::operator!=(Complex<Real> c) const {
    return real != c.real || imag != c.imag;
}

template <typename Real>
Complex<Real> Complex<Real>::operator+(Complex<Real> c) const {
    return Complex(real + c.real, imag + c.imag);
}

template <typename Real>
Complex<Real> Complex<Real>::operator-(Complex<Real> c) const {
    return Complex(real - c.real, imag - c.imag);
}

template <typename Real>
Complex<Real> Complex<Real>::operator-() const {
    return Complex(-real, -imag);
}

template <typename Real>
Complex<Real> Complex<Real>::operator*(Complex<Real> c) const {
    return Complex(real * c.real - imag * c.imag, real * c.imag + imag * c.real);
}

template <typename Real>
Complex<Real> Complex<Real>::operator/(Complex<Real> c) const {
    Real mag2 = c.mag2();
    return Complex((real * c.real + imag * c.imag) / mag2, (imag * c.real - real * c.imag) / mag2);
}

template <typename Real>
Complex<Real> Complex<Real>::conj() const {
    return Complex(real, -imag);
}

template <typename Real>
Real Complex<Real>::mag2() const {
    return real * real + imag * imag;
}

template <typename Real>
std::ostream &operator<<(std::ostream &os, Complex<Real> c) {
    os << c.real << "+i" << c.imag;
    return os;
}

}  // namespace linalg