#include <string>
#include <cstring>
#include <cstddef>
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>

#include "rdmini/rdmodel.h"
#include "rdmini/serial_ssa.h"
#include "rdmini/rdmini_version.h"

const char *demo_sim_version="0.0.1";

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
    "[OPTION] [model-file]\n"
    "  -m MODEL    Load the model named MODEL\n"
    "  -n N        Run simulation N steps\n"
    "  -t TIME     Run simulation for TIME simulated seconds\n"
    "  -d N/TIME   Sample simulation every N steps or TIME seconds\n"
    "  -v          Verbose output\n"
    "  -B          Batch output\n"
    "\n"
    "  -h          Print usage information\n"
    "  -V          Print version information\n"
    "\nOne of -n or -t must be specified.\n";

struct cl_args {
    std::string model_file;
    std::string model_name;
    double sample_delta=0;
    double t_end=0;
    size_t n_events=0;
    int verbosity=0;
    bool batch=false;

    bool help=false;
    bool version=false;
};

cl_args parse_cl_args(int argc,char **argv) {
    cl_args A;

    enum parse_state_enum { no_opt, opt_m, opt_n, opt_t, opt_d } parse_state = no_opt;
    bool has_opt_m=false;
    bool has_opt_n=false;
    bool has_opt_t=false;
    bool has_opt_d=false;
    bool has_file=false;

    int i=0;
    while (++i<argc) {
        const char *arg=argv[i];
        switch (parse_state) {
        case no_opt:
            if (arg[0]=='-') {
                switch (arg[1]) {
                case 'm':
                    parse_state=opt_m;
                    break;
                case 'n':
                    parse_state=opt_n;
                    break;
                case 't':
                    parse_state=opt_t;
                    break;
                case 'd':
                    parse_state=opt_d;
                    break;
                case 'v':
                    ++A.verbosity;
                    break;
                case 'B':
                    A.batch=true;
                    break;
                case 'h':
                    A.help=true; // and return!
                    return A;
                case 'V':
                    A.version=true; // and return!
                    return A;
                default:
                    throw usage_error("unrecognized option "+std::string(arg));
                }
            }
            else {
                if (has_file) throw usage_error("unexpected argument");
                A.model_file=arg;
                has_file=true;
            }
            break;
        case opt_m:
            if (has_opt_m)
                throw usage_error("-m specified multiple times");
            A.model_name=arg;
            has_opt_m=true;
            parse_state=no_opt;
            break;
        case opt_n:
            if (has_opt_n)
                throw usage_error("-n specified multiple times");
            A.n_events=std::stoull(arg);
            has_opt_n=true;
            parse_state=no_opt;
            break;
        case opt_t:
            if (has_opt_t)
                throw usage_error("-t specified multiple times");
            A.t_end=std::stod(arg);
            has_opt_t=true;
            parse_state=no_opt;
            break;
        case opt_d:
            if (has_opt_d)
                throw usage_error("-d specified multiple times");
            A.sample_delta=std::stod(arg);
            has_opt_d=true;
            parse_state=no_opt;
            break;
        }
    }

    if (parse_state!=no_opt)
        throw usage_error("missing option argument");

    return A;
}


struct emit_sim {
    explicit emit_sim(const rd_model &M, bool batch_=false, size_t expected_samples=0): n_species(M.n_species()), n_cells(M.n_cells()), batch(batch_) {
        // prepare csv-style header
        std::stringstream s;
        s << "time,cell";

        for (size_t i=0; i<n_species; ++i) 
            s << ',' << M.species[i].name;
        
        s << '\n';
        header=s.str();

        if (batch) {
            batch_n=0;
            batch_stride=n_species*n_cells;
            batch_t.reserve(expected_samples);
            batch_counts.reserve(batch_stride*expected_samples);
        }
    }
    
    std::ostream &emit_header(std::ostream &O) {
        return batch?O:O << header;
    }

    template <typename Sim>
    std::ostream &emit_state(std::ostream &O, double t, const Sim &sim) {
        if (!batch) {
            for (size_t i=0; i<n_cells; ++i) {
                O << t << ',' << i;
                for (size_t j=0; j<n_species; ++j)
                    O << ',' << sim.count(j,i);
                O << '\n';
            }
        }
        else {
            // consider exposing access to ssa population counts by population index directly
            batch_t.push_back(t);
            size_t offset=batch_stride*batch_n;
            batch_counts.resize(offset+batch_stride);
            ++batch_n;

            for (size_t i=0; i<n_cells; ++i)
                for (size_t j=0; j<n_species; ++j)
                    batch_counts[offset++]=sim.count(j,i);
        }
	return O;
    }

    std::ostream &flush(std::ostream &O) {
        if (!batch) return O;

        O << header;
        for (size_t k=0; k<batch_n; ++k) {
            size_t offset=batch_stride*k;
            double t=batch_t[k];
            for (size_t i=0; i<n_cells; ++i) {
                O << t << ',' << i;
                for (size_t j=0; j<n_species; ++j)
                    O << ',' << batch_counts[offset++];
                O << '\n';
            }
        }
        return O;
    }

    bool batch;
    size_t n_species,n_cells;
    std::string header;

    std::vector<double> batch_t;
    std::vector<size_t> batch_counts;
    size_t batch_stride;
    size_t batch_n;
};


int main(int argc, char **argv) {
    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
    std::chrono::duration<double,std::nano> elapsed_time;

    const char *basename=strrchr(argv[0],'/');
    basename=basename?basename+1:argv[0];
    int rc=0;

    try {
        cl_args A=parse_cl_args(argc,argv);

        if (A.help) {
            std::cout << "Usage: " << basename << " " << usage_text;
            return 0;
        }

        if (A.version) {
            std::cout << basename << " version " << demo_sim_version << "\n";
            std::cout << "rdmini library version " << rdmini_version << "\n";
            return 0;
        }

        // read in model specification

        rd_model M;

        if (A.model_file.empty() || A.model_file=="-")
            M=rd_model_read(std::cin,A.model_name);
        else {
            std::ifstream file(A.model_file);
            if (!file) throw fatal_error("unable to open file for reading");

            M=rd_model_read(file,A.model_name);
        }

        // set up simulator
            
        std::minstd_rand g;
        serial_ssa S(M,0);

        // run simulation

        size_t expected_samples=0;
        if (A.n_events>0) {
            if (A.sample_delta<1) A.sample_delta=1;
            expected_samples=1+A.n_events/(size_t)(A.sample_delta);
        }
        else {
            if (A.sample_delta==0) A.sample_delta=A.t_end;
            expected_samples=1+(size_t)(A.t_end/A.sample_delta);
        }

        emit_sim emitter(M,A.batch,expected_samples);
        emitter.emit_header(std::cout);

        emitter.emit_state(std::cout,0,S);
        if (A.verbosity) std::cout << S;


        if (A.n_events>0) { // event-by-event simulation
            double t;
            size_t dn=(size_t)A.sample_delta;

            start = std::chrono::high_resolution_clock::now();
            for (size_t n=0; n<A.n_events; n+=dn) {
                for (size_t i=0; i<dn; ++i)
                    t=S.advance(g);

                emitter.emit_state(std::cout,t,S);
                if (A.verbosity) std::cout << S;
            }
            end = std::chrono::high_resolution_clock::now();
        }
        else {
            double t=0;

            start = std::chrono::high_resolution_clock::now();
            while (t<A.t_end) {
                t=S.advance(t+A.sample_delta,g);

                emitter.emit_state(std::cout,t,S);
                if (A.verbosity) std::cout << S;
            }
            end = std::chrono::high_resolution_clock::now();
        }
        elapsed_time = end-start;
        emitter.flush(std::cout);
        std::cout << "---------------- \n";
        std::cout << "elapsed time: " << elapsed_time.count() << " [nano s] \n";
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

