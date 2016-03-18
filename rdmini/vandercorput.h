#ifndef VANDERCORPUT_H_
#define VANDERCORPUT_H_

#include <iostream>
#include <random>
#include <limits> 

//-----------------------------------------------------------------------------
class Linear_RNG {
    public:
        typedef unsigned result_type;

    private:
        result_type _s  = 0;
        static constexpr result_type Nmax = std::numeric_limits<unsigned int>::max();

    public:
        result_type operator()(){ return (_s++)%Nmax; }

        static result_type min() { return 0; }
        static result_type max() { return Nmax-1; }
};


//-----------------------------------------------------------------------------
template <class RealType = double > 
class VC_uniform_real_distribution {
    public:
        typedef RealType result_type;
        typedef RealType param_type;

    private:
        param_type _p;

    public:
        explicit VC_uniform_real_distribution(param_type p=1.0) : _p(p) { }

        void reset() {}

        param_type param() const { return _p; }

        void param(param_type p) { _p= p; }  

        template <class Generator>
            result_type operator()(Generator &g) {
                return (*this)(g,_p);
            }

        template <class Generator>
            result_type operator()(Generator &g, const param_type& p){
                result_type r;
                double k=0.1;
                unsigned n=g();
                for (r=0.0; n!=0; n/=10) {
                    r = r + (n%10)*k;
                    k/=10.;
                }
                return p*r;
            }

        result_type min() const { return 0; }
        result_type max() const { return _p; }
};



//-----------------------------------------------------------------------------
int main(int argc, const char ** argv){
    Linear_RNG rng;
    std::uniform_real_distribution<double> ud(0,1);
    VC_uniform_real_distribution<> vc(1);

    for (unsigned  i=1; i<100; ++i)
        std::cout << vc(rng) << std::endl;

    return 0;   
}

#endif //ndef VANDERCORPUT_H_ 
