/** Demo sampler implementations */

#include <stdexcept>
#include <cmath>
#include <random>
#include <vector>
#include <cstring>
#include <sstream>
#include <cassert>

#include "rdmini/sampler.h"
#include "rdmini/util/iterator.h"
#include "rdmini/timer.h"

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
    "[OPTION] SAMPLER\n"
    "Options:\n"
    "  -c N       Sample size\n"
    "  -N N       Population size\n"
    "  -g RATIO   Distribute expectation values geometrically\n"
    "  -l RATIO   Distribute expectation values linearly\n"
    "  -m MU1,... Set expectation values explicitly to MU1,...\n"
    "             (values will be scaled to sum to sample size)\n"
    "  -n N       Run N trials (default 1)\n"
    "  -d SEED    Initialise RNG seed to SEED\n"
    "  -o STAT    Emit statistics according to STAT (see below)\n"
    "  -T         Print timing data\n\n"
    "RATIO parameters describe the ratio between the first and last expectations.\n"
    "Mote that for without-replacement samplers, the expectations will equal the\n"
    "inclusion probabilities.\n\n"
    "SAMPLER is one of: multinomial, adjpareto, efraimidis, oss, cpsrej\n\n"
    "STAT is one of:\n"
    "    raw:    output result of each sample\n"
    "    mu:     ouptut mean across samples\n"
    "    pi:     ouptut empirical inclusion probability\n"
    "    pi2     output second order empirical inclusion probabilities\n";

enum sampler_enum { MULTINOMIAL, OSS, ADJPARETO, EFRAIMIDIS, CPSREJ, UNKNOWN_SAMPLER };
enum mu_enum { CONSTANT, GEOMETRIC, LINEAR, EXPLICIT };
enum stat_enum { RAW, MU, PI, PI2, UNKNOWN_STAT };

inline bool sampler_is_wr(enum sampler_enum s) {
    return s==MULTINOMIAL;
}

struct cl_args {
    int trials=1;    // number of trials
    int N=1;         // population size
    int c=1;         // sample size

    double ratio=1;          // inlusion probability ratio
    unsigned long seed=0;
    enum stat_enum stats=MU;

    mu_enum mu_spec=CONSTANT;
    sampler_enum sampler=UNKNOWN_SAMPLER;
    std::vector<double> mu;

    bool emit_timing=false; 
};

std::pair<const char *,sampler_enum> sampler_tbl[]={
    {"multinomial", MULTINOMIAL},
    {"oss", OSS},
    {"adjpareto", ADJPARETO},
    {"efraimidis", EFRAIMIDIS},
    {"cpsrej", CPSREJ},
    {0, UNKNOWN_SAMPLER}
};

std::pair<const char *,stat_enum> stat_tbl[]={
    {"raw", RAW},
    {"mu", MU},
    {"pi", PI},
    {"pi2", PI2},
    {0, UNKNOWN_STAT}
};

template <typename T>
T keyword_lookup(const std::pair<const char *,T> *tbl,const char *kw) {
    for (;;) {
        if (!tbl->first || !std::strcmp(tbl->first,kw)) return tbl->second;
        ++tbl;
    }
}

template <typename T>
std::vector<T> parse_csv_option(const std::string &s) {
    std::string error_string="failed to parse option argument: "+s;
    std::vector<T> values;
    std::stringstream ss(s);

    while (ss) {
        char comma;
        T x;

        ss >> x;
        if (!ss) throw usage_error(error_string);
  
        ss >> comma;
        if (comma!=',' && !ss.eof()) throw usage_error(error_string);

        values.push_back(std::stod(s));
    }
    if (!ss.eof()) throw usage_error(error_string);

    return values;
}

cl_args parse_cl_args(int argc,char **argv) {
    cl_args A;

    enum parse_state_enum { no_opt, opt_c, opt_N, opt_g, opt_l, opt_p, opt_n, opt_d, opt_o }
        parse_state = no_opt;

    int i=0;
    while (++i<argc) {
        const char *arg=argv[i];
        switch (parse_state) {
        case no_opt:
            if (arg[0]!='-') {
                // must be sampler name
                if (A.sampler!=UNKNOWN_SAMPLER)
                    throw usage_error("unexpected argument");

                A.sampler=keyword_lookup(sampler_tbl,arg);
                if (A.sampler==UNKNOWN_SAMPLER)
                    throw usage_error("unrecognized sampler "+std::string(arg));
            }
            else {
                switch (arg[1]) {
                case 'c': parse_state=opt_c; break;
                case 'N': parse_state=opt_N; break;
                case 'g': parse_state=opt_g; break;
                case 'l': parse_state=opt_l; break;
                case 'p': parse_state=opt_n; break;
                case 'n': parse_state=opt_n; break;
                case 'd': parse_state=opt_d; break;
                case 'o': parse_state=opt_o; break;
                case 'T':
                    A.emit_timing=true;
                    break;
                default:
                    throw usage_error("unrecognized option "+std::string(arg));
                }
            }
            break;
        case opt_c:
            A.c=std::stoi(arg);
            parse_state=no_opt;
            break;
        case opt_N:
            A.N=std::stoi(arg);
            parse_state=no_opt;
            break;
        case opt_g:
            A.ratio=std::stod(arg);
            A.mu_spec=GEOMETRIC;
            parse_state=no_opt;
            break;
        case opt_l:
            A.ratio=std::stod(arg);
            A.mu_spec=LINEAR;
            parse_state=no_opt;
            break;
        case opt_n:
            A.trials=std::stoi(arg);
            parse_state=no_opt;
            break;
        case opt_p:
            A.mu=parse_csv_option<double>(arg);
            A.mu_spec=EXPLICIT;
            parse_state=no_opt;
            break;
        case opt_d:
            A.seed=std::stoul(arg);
            parse_state=no_opt;
            break;
        case opt_o:
            A.stats=keyword_lookup(stat_tbl,arg);
            if (A.stats==UNKNOWN_STAT)
                throw usage_error("unrecognized output statistic "+std::string(arg));
            parse_state=no_opt;
            break;
        }
    }

    if (parse_state!=no_opt)
        throw usage_error("missing option argument");

    if (A.sampler==UNKNOWN_SAMPLER)
        throw usage_error("missing sampler argument");

    // compute expectations mu according to options.
    switch (A.mu_spec) {
    case CONSTANT:
        A.mu.assign(A.N,(double)A.c/A.N);
        break;
    case EXPLICIT:
        A.mu.resize(A.N);
        { 
            double scale=A.c/std::accumulate(A.mu.begin(),A.mu.end(),0.0);
            for (auto &mu: A.mu) mu*=scale;
        }
        break;
    case LINEAR:
        A.mu.resize(A.N);
        if (A.N>1) {
            double a=2.0/(A.N-1)*(A.ratio-1)/(A.ratio+1);
            for (size_t i=0; i<A.N; ++i) 
                A.mu[i]=(1+a*(i-(A.N-1)*0.5))*A.c/(double)A.N;
        }
        else A.mu.assign(A.N,1.0);
        break;
    case GEOMETRIC:
        A.mu.resize(A.N);
        if (A.N>1) {
            A.mu.resize(A.N);
            double a=std::pow(A.ratio,1.0/(A.N-1));
            A.mu[0]=A.c*(a-1)/(std::pow(a,A.N)-1);
            for (size_t i=1; i<A.N; ++i) A.mu[i]=a*A.mu[i-1];
        }
        else A.mu.assign(A.N,1.0);
        break;
    default:
        throw usage_error("error in parsing mu specification");
    }

    // check mu for out-of-range values
    
    for (auto mu: A.mu) 
        if (mu<0) throw usage_error("negative expectation specified");

    if (!sampler_is_wr(A.sampler)) {
        double mu_max=*std::max_element(A.mu.begin(),A.mu.end());
        if (mu_max>1) {
            if (A.mu_spec==EXPLICIT)
                throw usage_error("maximum expectation for a without-replacement sampler is 1.");

            double min_ratio=1,max_ratio=1;
            switch (A.mu_spec) {
            case LINEAR:
                min_ratio=2.0*A.c/(A.N-1.0);
                max_ratio=1/min_ratio;
                break;
            case GEOMETRIC:
                // Okay, this is a bit hairy. We do Newton a few times
                // starting from twice the stationary point.
                {
                    double n=A.c;
                    double N=A.N;
                    double x=2*std::pow(n/(n-1)*(N-1)/N,N-1);
                    double a=std::pow(x,1/(N-1));
                    std::cout << "x=" << x << "; a=" << a << "\n";
                    for (int i=0; i<5; ++i) {
                        x-=((n-1)*x*a-n*x+1)/(N*(n-1)/(N-1)*a-n);
                        a=std::pow(x,1/(N-1));
                        std::cout << "x=" << x << "; a=" << a << "\n";
                    }
                    max_ratio=x;
                    min_ratio=1/x;
                }
                break;
            default:
                throw std::logic_error("internal error");
            }
            std::stringstream msg;
            msg << "expectation over one for a without-replacement sampler; "
                << "for sampling " << A.c << " from " << A.N << ", valid ratio range is "
                << std::setprecision(4) << min_ratio << " to " << max_ratio << ".";
    
            throw usage_error(msg.str());
        }
    }

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

// Sample with a draw-sampler, return vector of counts across population.

template <typename DSampler,typename RNG>
std::vector<unsigned> sample_draw(unsigned N,DSampler &S,RNG &R)
{
    std::vector<unsigned> sample(N);
    S.sample(sample.begin(),sample.end(),rdmini::functor_iterator([](unsigned &b) { ++b; }),R);
    return sample;
}

// Sample with a rejection or reservoir sampler, return vector of counts across population.

template <typename RSampler,typename RNG>
std::vector<unsigned> sample_rr(unsigned N,RSampler &S,RNG &R)
{
    using count=rdmini::counting_iterator<unsigned>;

    std::vector<unsigned> sample(N);
    std::vector<unsigned> items(S.max());

    typename RSampler::size_type n=S.sample(count(0),count(N),items.begin(),R);

    items.resize(n);
    for (auto i: items) ++sample[i];
    return sample;
}

// running stats

struct running_mean {
    running_mean() { clear(); }

    double mean() const { return m; }

    void clear() { 
        n=0;
        m=0;
    }

    void insert(double x) { m+=(x-m)/++n; }

    int n;
    double m;
};


// harness

void run_test(const cl_args &A) {
    std::mt19937_64 R(A.seed);
    unsigned N=A.N;

    // print header and initialise stats vectors
    std::vector<running_mean> means,means2;

    switch (A.stats) {
    case RAW:
        std::cout << "trial";
        for (int i=0; i<N; ++i) std::cout << ",i" << (1+i);
        std::cout << '\n';
        break;
    case MU:
        std::cout << "item,model_mu,mu\n";
        means.resize(N);
        break;
    case PI:
        std::cout << "item,model_pi,pi\n";
        means.resize(N);
        break;
    case PI2:
        std::cout << "item,model_pi";
        for (int i=0; i<N; ++i) std::cout << ",pi" << (1+i);
        std::cout << '\n';
        means.resize(N);
        means2.resize((N*(N-1))/2);
        break;
    }

    // Compute model pi for with-replacement samplers
    // (currently only multinomial!)
    std::vector<double> model_pi;
    if (!sampler_is_wr(A.sampler))
        model_pi=A.mu;
    else {
        switch (A.sampler) {
        case MULTINOMIAL:
            model_pi=A.mu;
            for (auto &pi: model_pi) pi=1-std::pow((1-pi/A.c),A.c);
            break;
        default:
            throw std::runtime_error("unable to compute inclusion probabilities for this sampler");
        }
    }

    rdmini::timer::hr_timer timer;

    for (int trial=0; trial<A.trials; ++trial) {
        std::vector<unsigned> sample;
        timer.resume();
        switch (A.sampler) {
        case MULTINOMIAL:
            {
                rdmini::multinomial_draw_sampler S(A.c,A.mu.begin(),A.mu.end());
                sample=sample_draw(N,S,R);
            }
            break;
        case OSS:
            {
                rdmini::ordered_systematic_sampler S(A.mu.begin(),A.mu.end());
                sample=sample_draw(N,S,R);
            }
            break;
        case ADJPARETO:
            {
                rdmini::adjusted_pareto_sampler S(A.c,A.mu.begin(),A.mu.end());
                sample=sample_rr(N,S,R);
            }
            break;
        case EFRAIMIDIS:
            {
                rdmini::efraimidis_spirakis_sampler S(A.c,A.mu.begin(),A.mu.end());
                sample=sample_rr(N,S,R);
            }
            break;
        case CPSREJ:
            {
                rdmini::cps_multinomial_rejective S(A.c,A.mu.begin(),A.mu.end());
                sample=sample_rr(N,S,R);
            }
            break;
        default:
            throw fatal_error("unrecognized sampler");
        }
        timer.stop();

        switch (A.stats) {
        case RAW:
            std::cout << (1+trial);
            for (size_t i=0; i<A.N; ++i)
                std::cout << ',' << sample[i];
            std::cout << "\n";
            break;
        case MU:
            for (unsigned i=0; i<N; ++i) means[i].insert(sample[i]);
            break;
        case PI:
            for (unsigned i=0; i<N; ++i) means[i].insert(sample[i]>0);
            break;
        case PI2:
            for (unsigned i=0; i<N; ++i) means[i].insert(sample[i]>0);
            for (unsigned i=1; i<N; ++i) {
                for (unsigned j=0; j<i; ++j) {
                    unsigned index=(i*(i-1))/2+j;
                    means2[index].insert(sample[i]>0 && sample[j]>0);
                }
            }
            break;
        }
    }

    switch (A.stats) {
    case RAW:
        break;
    case MU:
        for (unsigned i=0; i<N; ++i)
            std::cout << (i+1) << "," << A.mu[i] << "," << means[i].mean() << "\n";
        break;
    case PI:
        for (unsigned i=0; i<N; ++i) 
            std::cout << (i+1) << "," << model_pi[i] << "," << means[i].mean() << "\n";
        break;
    case PI2:
        for (unsigned i=0; i<N; ++i) {
            std::cout << (i+1) << "," << model_pi[i];
            for (unsigned j=0; j<N; ++j) {
                double x;
                if (i==j)
                    x=means[i].mean();
                else if (j<i)
                    x=means2[(i*(i-1))/2+j].mean();
                else
                    x=means2[(j*(j-1))/2+i].mean();

                std::cout << ',' << x;;
            }
            std::cout << "\n";
        }
        break;
    }

    if (A.emit_timing) {
        double us=1.e6*timer.time()/A.trials;
        std::cout << "mean execution time [µs]: " << std::setprecision(3) << std::fixed << us << "\n";
    }
}
