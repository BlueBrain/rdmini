#ifndef VANDERCORPUT_H_
#define VANDERCORPUT_H_

#include <iostream>
#include <random>
#include <limits> 

//-----------------------------------------------------------------------------
namespace rdmini{
class counting_generator {
    public:
        typedef unsigned result_type;

    private:
        result_type s  = 0;

    public:
        result_type operator()(){ return s++; }

        static result_type min() { return 0; }
        static result_type max() { return std::numeric_limits<result_type>::max(); }
};


//-----------------------------------------------------------------------------
template <class RealType = double > 
class vdc_uniform_real_distribution {
    public:
        typedef RealType result_type;

        // Nested class: param_type
        class param_type {
        private:
            RealType M_a;
            RealType M_b;
        public: 
            explicit param_type(RealType a=0, RealType b=1) : M_a(a), M_b(b) {}
            result_type a() const { return M_a;}
            result_type b() const { return M_b;}
        };

        // Constructor (1)
        explicit vdc_uniform_real_distribution(RealType a=0, RealType b=1) : M_param(a,b) {}
        
        // Constructor (2)
        explicit vdc_uniform_real_distribution(const param_type& p) : M_param(p) {}

        // Reset the distribution rate (nothing to be done)
        void reset() {}

        // Getter of the parameter
        param_type param() const { return M_param; }

        // Setter of the parameter
        void param(const param_type& p) { M_param = p; }  

        // Overload call operator (1)
        template <class Generator>
            result_type operator()(Generator &g) {
                return (*this)(g,M_param);
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
        result_type min() const { return M_param.a(); }
        result_type max() const { return M_param.b(); }

    private:
        param_type M_param;
};
}

#endif //ndef VANDERCORPUT_H_ 
