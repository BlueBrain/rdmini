#ifndef VANDERCORPUT_H_
#define VANDERCORPUT_H_

#include <iostream>
#include <random>
#include <limits> 

//-----------------------------------------------------------------------------
namespace rdmini{
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
class vc_uniform_real_distribution {
    public:
        typedef RealType result_type;

        // Nested class: param_type
        class param_type {
        private:
            RealType _M_a;
            RealType _M_b;
        public: 
            explicit param_type(RealType a=0, RealType b=1) : _M_a(a), _M_b(b) {}
            result_type a() const { return _M_a;}
            result_type b() const { return _M_b;}
        };

        // Constructor (1)
        explicit vc_uniform_real_distribution(RealType a=0, RealType b=1) : _M_param(a,b) {}
        
        // Constructor (2)
        explicit vc_uniform_real_distribution(const param_type& p) : _M_param(p) {}

        // Reset the distribution rate (nothing to be done)
        void reset() {}

        // Getter of the parameter
        param_type param() const { return _M_param; }

        // Setter of the parameter
        void param(const param_type& p) { _M_param = p; }  

        // Overload call operator (1)
        template <class Generator>
            result_type operator()(Generator &g) {
                return (*this)(g,_M_param);
            }
        
        // Overload call operator (2)
        template <class Generator>
            result_type operator()(Generator &g, const param_type& p){
                result_type r;
                double k=0.1;
                unsigned n=g();
                for (r=0.0; n!=0; n/=10) {
                    r = r + (n%10)*k;
                    k/=10.;
                }
                return (p.b()-p.a())*r + p.a();
            }

        // Getters of min and max possible values
        result_type min() const { return _M_param.a(); }
        result_type max() const { return _M_param.b(); }

    private:
        param_type _M_param;
};
}

#endif //ndef VANDERCORPUT_H_ 
