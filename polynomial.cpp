#include <iostream>
#include <cassert>
#include <vector>
<<<<<<< HEAD


class Polynomial{
	public:

=======
#include <algorithm>

class Polynomial
{
    // Coefficient vector for monom basis
    std::vector<int> a_;
>>>>>>> 5f527ff16714fca7f00bf949c8eb8b2e19c7e4c5

    // Coefficient vector for newton basis
    std::vector<int> n_;



  public:

    // 1-D. Conversion to 2-D is necessary
    Polynomial()
    : a_(1, int())
    {}

    Polynomial(int n)
    : a_(std::max(1, n+1), int())
    {}

    Polynomial(std::initializer_list<int> coeffs)
    : a_{coeffs}
    {}

    int degree() const
    {
        return a_.size()-1;
    }

    bool operator==(Polynomial const & other) const
    {
        return a_ == other.a_;
    }

    int operator[](int i) const
    {
        return a_[i];
    }

    int & operator[](int i)
    {
        return a_[i];
    }

    int operator()(int x) const
    {
        int res = a_[degree()];
        for(int k=degree()-1; k>=0; --k)
        {
            res = res*x + a_[k];
        }
        return res;
    }

    Polynomial derivative() const
    {
        Polynomial res(degree()-1);
        for(int k=1; k<=degree(); ++k)
            res[k-1] = k*a_[k];
        return res;
    }

    Polynomial derivative(int d) const
    {
        Polynomial res(*this);
        for(int k=0; k<d; ++k)
            res = res.derivative();
        return res;
    }

    Polynomial operator+(Polynomial const & other) const
    {
        if(degree() < other.degree())
        {
            Polynomial res(other);
            for(int k=0; k<=degree(); ++k)
                res[k] += a_[k];
            return res;
        }
        else
        {
            Polynomial res(*this);
            for(int k=0; k<=other.degree(); ++k)
                res[k] += other.a_[k];
            return res;
        }
    }

    Polynomial operator*(Polynomial const & other) const
    {
        Polynomial res(degree()+other.degree());
        for(int k=0; k<=degree(); ++k)
            for(int l=0; l<=other.degree(); ++l)
                res[k+l] += a_[k]*other.a_[l];
        return res;
    }
};

int main()
{
    // Testen des Standard-Konstruktors
    Polynomial p0;
    assert(p0.degree() == 0);
    assert(p0[0] == 0);

    // Testen des Konstruktors mit Gradangabe
    Polynomial pm(-1);
    assert(pm == p0);

    Polynomial p2(2);
    assert(p2.degree() == 2);

    // Testen des Konstruktors mit Koeffizientenarray
    Polynomial p{1, 2, 3};
    assert(p.degree() == 2);

    // Testen der Werte der Koeffizienten
    assert(p2[0] == 0);
    assert(p2[1] == 0);
    assert(p2[2] == 0);

    assert(p[0] == 1);
    assert(p[1] == 2);
    assert(p[2] == 3);

    // Testen, dass die Polynom fuer verschiedene Argumente x
    // korrekt ausgerechnet werden
    assert(p0(1) == 0);
    assert(p0(2) == 0);
    assert(p0(3) == 0);

    assert(p2(1) == 0);
    assert(p2(2) == 0);
    assert(p2(3) == 0);

    assert(p(1) == 6);
    assert(p(2) == 17);
    assert(p(3) == 34);

    // Testen der Polynom-Addition
    Polynomial pa1 = p + p;
    Polynomial pa1_expected{2, 4, 6};
    assert(pa1 == pa1_expected);

    Polynomial pa2 = p + Polynomial{3, 6};
    Polynomial pa2_expected{4, 8, 3};
    assert(pa2 == pa2_expected);

    // Testen der Polynom-Multiplikation
    Polynomial pm1 = p * p;
    assert(pm1.degree() == 4);
    assert(pm1(2) == p(2)*p(2));
    Polynomial pm1_expected{1, 4, 10, 12, 9};
    assert(pm1 == pm1_expected);

    Polynomial pm2 = p * Polynomial{2, 6};
    assert(pm2.degree() == 3);
    Polynomial pm2_expected{2, 10, 18, 18};
    assert(pm2 == pm2_expected);

    // Testen der ersten Ableitung
    Polynomial derivative1_expected{2, 6};
    assert(p.derivative()  == derivative1_expected);
    assert(p.derivative(1) == derivative1_expected);

    // Testen der zweiten bis vierten Ableitung
    Polynomial derivative2_expected{6};
    assert(p.derivative(2)  == derivative2_expected);
    Polynomial derivative3_expected{0};
    assert(p.derivative(3)  == derivative3_expected);
    assert(p.derivative(4)  == derivative3_expected);

    // Testen der 

<<<<<<< HEAD


int main(){


	return 0;
}
=======
    std::cout << "Alle Tests erfolgreich!\n";
}

>>>>>>> 5f527ff16714fca7f00bf949c8eb8b2e19c7e4c5
