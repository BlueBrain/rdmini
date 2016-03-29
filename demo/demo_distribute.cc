/** Test various distribution procedures */

#include <stdexcept>
#include <cmath>
#include <random>
#include <vector>
#include <cstring>
#include <cassert>

#include "rdmini/categorical.h"
#include "rdmini/sampler.h"
#include "rdmini/timer.h"
#include "rdmini/util/functor_iterator.h"

using std::size_t;

struct fatal_error: std::exception {
    fatal_error(const std::string &what_str_): what_str(what_str_) {}
    const char *what() const throw() { return what_str.c_str(); }

private:
    std::string what_str;
};
    
struct usage_error: fatal_error {
    usage_error(const std::string &what_str_): fatal_error(what_str_) {}
};

const char *usage_text=
    "[OPTION]\n"
    "  -m METHOD Method to use: steps, multinomial, oss\n"
    "  -c N      Count to distribute\n"
    "  -c N-M    Select counts uniformly in interval [N,M]\n"
    "  -b N      Distribute among N bins\n"
    "  -g RATIO  Distribute weights geometrically by RATIO\n"
    "  -n N      Run N trials (default 1)\n"
    "  -C        Report raw counts, not normalised values\n"
    "  -S        Print just summary statistics\n\n"
    "  -T        Print timing statistics\n\n"
    "The multinomial method assigns rounded-down values and then\n"
    "distributes the remainder multinomially.\n"
    "The oss method assigns rounded-down and then distributes the\n"
    "remainder by ordered systematic sampling.\n\n"
    "Normalised results are scaled by inverse bin weight; weights\n"
    "are scaled so that the total weight is the number of bins.\n";

enum method_enum { STEPS, MULTINOMIAL, OSS, UNKNOWN_METHOD };

struct cl_args {
    int N=1;    // number of trials
    int b=1;    // number of bins
    std::pair<int,int> c={1,1};  // (inclusive) range of counts to distribute

    bool summary=false; // print just summary statistics
    bool raw_counts=false; // don't scale by weight in report
    bool emit_timing=false;
    double weight_ratio=1; // ratio in weights between successive bins

    method_enum method=STEPS;
};

std::pair<const char *,method_enum> method_tbl[]={
    {"steps", STEPS},
    {"multinomial", MULTINOMIAL},
    {"oss", OSS},
    {0, UNKNOWN_METHOD}
};

template <typename T>
T keyword_lookup(const std::pair<const char *,T> *tbl,const char *kw) {
    for (;;) {
        if (!tbl->first || !std::strcmp(tbl->first,kw)) return tbl->second;
        ++tbl;
    }
}

std::pair<int,int> parse_range(const char *s) {
    char *s2=0;
    int cmin=0,cmax=0;

    cmin=(int)std::strtol(s,&s2,0);
    if (*s2 && *s2!='-') throw usage_error(std::string("failed to parse range: ")+s);

    if (!*s2) return std::pair<int,int>(cmin,cmin);
    ++s2;

    cmax=(int)std::strtol(s2,&s2,0);
    if (*s2 || cmax<cmin) throw usage_error(std::string("failed to parse range: ")+s);
    return std::pair<int,int>(cmin,cmax);
}

cl_args parse_cl_args(int argc,char **argv) {
    cl_args A;

    enum parse_state_enum { no_opt, opt_c, opt_b, opt_g, opt_n, opt_m } parse_state = no_opt;

    int i=0;
    while (++i<argc) {
        const char *arg=argv[i];
        switch (parse_state) {
        case no_opt:
            if (arg[0]!='-') throw usage_error("unexpected argument");
            switch (arg[1]) {
            case 'c': parse_state=opt_c; break;
            case 'b': parse_state=opt_b; break;
            case 'g': parse_state=opt_g; break;
            case 'n': parse_state=opt_n; break;
            case 'm': parse_state=opt_m; break;
            case 'C': // no argument
                A.raw_counts=true;
                break;
            case 'T': // no argument
                A.emit_timing=true;
                break;
            case 'S': // no argument
                A.summary=true;
                break;
            default:
                throw usage_error("unrecognized option "+std::string(arg));
            }
            break;
        case opt_c:
            A.c=parse_range(arg);
            parse_state=no_opt;
            break;
        case opt_b:
            A.b=std::stoi(arg);
            parse_state=no_opt;
            break;
        case opt_g:
            A.weight_ratio=std::stod(arg);
            parse_state=no_opt;
            break;
        case opt_n:
            A.N=std::stoi(arg);
            parse_state=no_opt;
            break;
        case opt_m:
            A.method=keyword_lookup(method_tbl,arg);
            if (A.method==UNKNOWN_METHOD)
                throw usage_error("unrecognized method "+std::string(arg));
            parse_state=no_opt;
            break;
        }
    }

    if (parse_state!=no_opt)
        throw usage_error("missing option argument");

    return A;
}

void run_test(const cl_args &);

int main(int argc,char **argv) {
    char *basename=strrchr(argv[0],'/');
    basename=basename?basename+1:argv[0];
    int rc=0;

    try {
        cl_args A=parse_cl_args(argc,argv);

        run_test(A);
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

// Adapted from STEPS implementation:

template <typename Rng>
void distribute_steps(unsigned c, Rng &R,
                      std::vector<unsigned> &bin,
                      const std::vector<double> &weight)
{
    std::uniform_real_distribution<double> U;

    assert(bin.size()==weight.size());
    if (bin.empty()) return;

    std::fill(bin.begin(), bin.end(), 0);

    double total_weight=std::accumulate(weight.begin(),weight.end(),0.0);
    assert(total_weight > 0);

    unsigned nremoved=0;
    size_t nelem=bin.size();
    for (size_t i=0; i<nelem; ++i) {
        if (c==0 || nremoved == c) {
            bin[i]=0;
            continue;
        }

        double fract=c*(weight[i]/total_weight);
        unsigned n3=static_cast<unsigned>(fract);

        double n3_frac=fract-n3;
        if (n3_frac>0.0) {
            if (U(R)<n3_frac) ++n3;
        }

        nremoved+=n3;
        if (nremoved>=c) {
            n3-=nremoved-c;
            nremoved = c;
        }

        bin[i]=n3;
    }

    assert(nremoved <= c);
    c -= nremoved;
    for (; c != 0; --c) {
        // pick item by weighted selection
        double selector = U(R) * total_weight;
        
        for (size_t i = 0; i < nelem; ++i) {
            selector -= weight[i];
            if (selector < 0) {
                ++bin[i];
                break;
            }
        }
    }
}

// Distribute rounded-down values, return deficit and write residuals back
// into weights.
size_t distribute_common(unsigned c,std::vector<unsigned> &bin,
                         std::vector<double> &weight) {
    size_t N=bin.size();
    if (!N) return 0;
    if (N!=weight.size()) throw fatal_error("bin and weight sizes differ");

    std::fill(bin.begin(), bin.end(), 0);
    double total_weight = std::accumulate(weight.begin(), weight.end(), 0.0);

    double oo_total=1.0/total_weight;
    size_t asum=0;
    for (size_t i=0; i<N; ++i) {
        double q=weight[i]*oo_total*c;
        size_t a=(size_t)q;

        bin[i]=a;
        weight[i]=q-(double)a;
        asum+=a;
    }
    assert(asum<=c);

    return c-asum;
}

template <typename RNG>
void distribute_multinomial(unsigned c, RNG &R,
                      std::vector<unsigned> &bin,
                      std::vector<double> weight)
{
    static std::uniform_real_distribution<double> U(0,1);

    size_t r=distribute_common(c,bin,weight);
    if (r==0) return;

    rdmini::multinomial_draw_sampler S(r,weight.begin(),weight.end());
    S.sample(bin.begin(),bin.end(),rdmini::functor_iterator([](unsigned &b) { ++b; }),R);
}

template <typename RNG>
void distribute_oss(unsigned c, RNG &R,
                      std::vector<unsigned> &bin,
                      std::vector<double> weight,
                      double total_weight = 0)
{
    static std::uniform_real_distribution<double> U(0,1);

    size_t r=distribute_common(c,bin,weight);
    if (r==0) return;

    rdmini::ordered_systematic_sampler S(weight.begin(),weight.end());
    S.sample(bin.begin(),bin.end(),rdmini::functor_iterator([](unsigned &b) { ++b; }),R);
}

// running stats

struct running_stats {
    running_stats() { clear(); }

    double mean() const { return m; }
    double variance() const { return n<2?0:m2/(n-1); }
    double min() const { return xmin; }
    double max() const { return xmax; }

    void clear() { 
        n=0;
        m=0;
        m2=0;
        xmin=0;
        xmax=0;
    }

    void insert(double x) {
        double s=x-m;

        ++n;
        m+=s/n;
        m2+=s*(x-m);

        if (!n || xmin>x) xmin=x;
        if (!n || xmax<x) xmax=x;
    }

    int n;
    double m,m2;
    double xmin,xmax;
};

// harness

void run_test(const cl_args &A) {
    std::mt19937_64 R;

    // print header
    if (!A.summary) {
        std::cout << "trial";
        for (int bin=0; bin<A.b; ++bin) std::cout << ",B" << (1+bin);
        std::cout << '\n';
    }

    std::vector<double> weight(A.b);
    if (A.weight_ratio==1) {
        std::fill(weight.begin(),weight.end(),1.0);
    }
    else {
        weight[0]=A.b*(A.weight_ratio-1)/(std::pow(A.weight_ratio,A.b)-1);
        for (size_t i=1; i<A.b; ++i) weight[i]=A.weight_ratio*weight[i-1];
    }

    std::vector<unsigned> bin(A.b,0);
    std::vector<running_stats> stats(A.b);

    std::uniform_int_distribution<int> U(A.c.first,A.c.second);
    rdmini::timer::hr_timer timer;

    for (int trial=0; trial<A.N; ++trial) {
        int count=U(R);
        timer.resume();
        switch (A.method) {
        case STEPS:
            distribute_steps(count,R,bin,weight);
            break;
        case MULTINOMIAL:
            distribute_multinomial(count,R,bin,weight);
            break;
        case OSS:
            distribute_oss(count,R,bin,weight);
            break;
        default:
            throw fatal_error("unrecognized method");
        }
        timer.stop();

        if (A.summary) {
            for (size_t i=0; i<A.b; ++i) {
                double x=bin[i];
                if (!A.raw_counts && weight[i]!=0) x/=weight[i];
                stats[i].insert(x);
            }
        }
        else {
            std::cout << (1+trial);
            for (size_t i=0; i<A.b; ++i) {
                std::cout << ',';
                if (A.raw_counts)
                    std::cout << bin[i];
                else
                    std::cout << (weight[i]!=0?bin[i]/weight[i]:bin[i]);
            }
            std::cout << '\n';
        }
    }

    if (A.summary) {
        std::cout << "bin,mean,variance";
        if (A.raw_counts) std::cout << ",min,max";
        std::cout << "\n";

        for (size_t i=0; i<A.b; ++i) {
            const auto &S=stats[i];

            std::cout << (i+1) << ',' << S.mean() << ',' << S.variance();
            if (A.raw_counts) std::cout << ',' << S.min() << ',' << S.max();
            std::cout << '\n';
        }
    }

    if (A.emit_timing) {
        double us=1.e6*timer.time()/A.N;
        std::cout << "mean execution time [µs]: " << std::setprecision(3) << std::fixed << us << "\n";
    }
}
