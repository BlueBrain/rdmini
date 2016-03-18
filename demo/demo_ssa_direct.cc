#include <string>
#include <cstring>
#include <cstddef>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <vector>
#include <random>

#include "rdmini/ssa_direct.h"

// throw to clean-up and exit
struct fatal_error: std::exception {
    fatal_error(const std::string &what_str_): what_str(what_str_) {}
    const char *what() const throw() { return what_str.c_str(); }

private:
    std::string what_str;
};
    
// throw on cmdline argument error
struct usage_error: fatal_error {
    usage_error(const std::string &what_str_): fatal_error(what_str_) {}
};

// usage info
const char *usage_text=
    "[OPTION]\n"
    "  -k K    Number of processes\n"
    "  -n N    Number of samples\n";

struct cl_args {
    int k=10; // number of processes
    int n=10000; // number of samples
};

cl_args parse_cl_args(int argc,char **argv) {
    cl_args A;

    enum parse_state_enum { no_opt, opt_n, opt_k } parse_state = no_opt;

    int i=0;
    while (++i<argc) {
        const char *arg=argv[i];
        switch (parse_state) {
        case no_opt:
            if (arg[0]!='-') throw usage_error("unexpected argument");
            switch (arg[1]) {
            case 'n': parse_state=opt_n; break;
            case 'k': parse_state=opt_k; break;
            default:
                throw usage_error("unrecognized option "+std::string(arg));
            }
            break;
        case opt_n:
            A.n=std::stoi(arg);
            parse_state=no_opt;
            break;
        case opt_k:
            A.k=std::stoi(arg);
            parse_state=no_opt;
            break;
        }
    }

    if (parse_state!=no_opt)
        throw usage_error("missing option argument");

    return A;
}

void run_demo_static(std::ostream &O,const cl_args &A) {
    std::minstd_rand R;
    std::uniform_real_distribution<double> U(0.5,1);

    // generate propensities
    size_t n_proc=A.k;
    size_t n_samp=A.n;

    std::vector<double> prop(n_proc);
    for (size_t i=0; i<n_proc; ++i) {
        prop[i]=std::ldexp(U(R),-i);
    }
    std::shuffle(prop.begin(),prop.end(),R);

    double total=std::accumulate(prop.begin(),prop.end(),0.0);
    
    O << "#propensities\nk,p\n";
    for (size_t i=0; i<prop.size(); ++i)
        O << i << ',' << prop[i] << "\n";

    O << "\n";

    typedef rdmini::ssa_direct<size_t,double> S;
    S ssa(prop.size());
    for (size_t i=0; i<prop.size(); ++i)
        ssa.update(i,prop[i]);

    std::vector<S::event_type> trace;
    for (size_t j=0; j<n_samp; ++j)
        trace.emplace_back(ssa.next(R));

    O << "#samples\nt,k\n";
    double t=0;
    for (auto ev: trace)
        O << (t+=ev.dt()) << ',' << ev.key() << "\n";
}


int main(int argc, char **argv) {
    const char *basename=strrchr(argv[0],'/');
    basename=basename?basename+1:argv[0];
    int rc=0;

    try {
        cl_args A=parse_cl_args(argc,argv);

        run_demo_static(std::cout,A);
    }
    catch (usage_error &E) {
        std::cerr << basename << ": " << E.what() << "\n";
        std::cerr << "Usage: " << basename << " " << usage_text;
        rc=2;
    }
    catch (std::exception &E) {
        std::cerr << basename << ": " << E.what() << "\n";
        rc=1;
    }
    
    return rc;
}

