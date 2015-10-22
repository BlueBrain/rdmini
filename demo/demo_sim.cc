#include <string>
#include <cstring>
#include <cstddef>
#include <iostream>
#include <fstream>
#include <sstream>

#include "rdmini/timer.h"
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
    "  -P N        Run N independent instances\n"
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
    int n_instances=1;

    bool help=false;
    bool version=false;
};

cl_args parse_cl_args(int argc,char **argv) {
    cl_args A;

    enum parse_state_enum { no_opt, opt_m, opt_n, opt_t, opt_d, opt_P } parse_state = no_opt;
    bool has_opt_m=false;
    bool has_opt_n=false;
    bool has_opt_t=false;
    bool has_opt_d=false;
    bool has_opt_P=false;
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
                case 'P':
                    parse_state=opt_P;
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
        case opt_P:
            if (has_opt_P)
                throw usage_error("-P specified multiple times");
            A.n_instances=std::stoi(arg);
            has_opt_P=true;
            parse_state=no_opt;
            break;
        }
    }

    if (parse_state!=no_opt)
        throw usage_error("missing option argument");

    return A;
}


struct emit_sim {
    explicit emit_sim(const rd_model &M, size_t ni, bool batch_=false, size_t expected_samples=0): n_species(M.n_species()), n_cells(M.n_cells()), n_instances(ni), batch(batch_) {
        // prepare csv-style header
        std::stringstream s;
        s << "instance,time,cell";

        for (size_t i=0; i<n_species; ++i) 
            s << ',' << M.species[i].name;
        
        s << '\n';
        header=s.str();

        if (batch) {
            batch_sample_width=n_species*n_cells;
            batch_samples.reserve(n_instances*expected_samples);
            batch_count_data.reserve(n_instances*batch_sample_width*expected_samples);
        }
    }
    
    std::ostream &emit_header(std::ostream &O) {
        return batch?O:O << header;
    }

    template <typename Sim>
    std::ostream &emit_state(std::ostream &O, size_t instance, double t, const Sim &sim) {
        if (!batch) {
            for (size_t i=0; i<n_cells; ++i) {
                O << instance << ',' << t << ',' << i;
                for (size_t j=0; j<n_species; ++j)
                    O << ',' << sim.count(j,i);
                O << '\n';
            }
        }
        else {
            // consider exposing access to ssa population counts by population index directly
            size_t offset=batch_count_data.size();
            batch_sample b={instance, t, offset}; 
            batch_samples.push_back(b);

            batch_count_data.resize(offset+batch_sample_width);

            for (size_t i=0; i<n_cells; ++i)
                for (size_t j=0; j<n_species; ++j)
                    batch_count_data[offset++]=sim.count(j,i);
        }
	return O;
    }

    std::ostream &flush(std::ostream &O) {
        if (!batch) return O;

        O << header;
        for (const auto &sample: batch_samples) {
            size_t offset=sample.count_data_offset;
            for (size_t i=0; i<n_cells; ++i) {
                O << sample.instance << ',' << sample.t << ',' << i;
                for (size_t j=0; j<n_species; ++j)
                    O << ',' << batch_count_data[offset++];
                O << '\n';
            }
        }
        return O;
    }

    bool batch;
    size_t n_species,n_cells,n_instances;
    std::string header;

    struct batch_sample {
        size_t instance;
        double t;
        size_t count_data_offset;
    };
    size_t batch_sample_width;

    std::vector<batch_sample> batch_samples;
    std::vector<size_t> batch_count_data;
};

void run_sim_by_steps(serial_ssa &S,emit_sim &emitter,size_t n,size_t dn,bool verbose) {
    std::minstd_rand g;

    double t;
    for (size_t i=0; i<n; i+=dn) {
        for (size_t j=0; j<dn; ++j)
            t=S.advance(g);

        emitter.emit_state(std::cout,1,t,S);
        if (verbose) std::cout << S;
    }
}

void run_sim_by_time(serial_ssa &S,emit_sim &emitter,double t_end,double dt,bool verbose) {
    std::minstd_rand g;

    double t=0;
    while (t<t_end) {
        t=S.advance(t+dt,g);

        emitter.emit_state(std::cout,1,t,S);
        if (verbose) std::cout << S;
    }
}


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

        // set up data emitter and timer

        timer::hr_timer T;

        size_t expected_samples=0;
        if (A.n_events>0) {
            if (A.sample_delta<1) A.sample_delta=1;
            expected_samples=1+A.n_events/(size_t)(A.sample_delta);
        }
        else {
            if (A.sample_delta==0) A.sample_delta=A.t_end;
            expected_samples=1+(size_t)(A.t_end/A.sample_delta);
        }

        emit_sim emitter(M,A.n_instances,A.batch,expected_samples);
        emitter.emit_header(std::cout);

        // set up simulator
            
        serial_ssa S(M,0);

        // emit initial state

        emitter.emit_state(std::cout,1,0,S);
        if (A.verbosity) std::cout << S;

        // run simulation

        if (A.n_events>0) {
            auto _(timer::guard(T));
            run_sim_by_steps(S,emitter,A.n_events,(size_t)A.sample_delta,A.verbosity>0);
        }
        else {
            auto _(timer::guard(T));
            run_sim_by_time(S,emitter,A.t_end,A.sample_delta,A.verbosity>0);
        }
        emitter.flush(std::cout);

        std::cout << "#---------------- \n";
        std::cout << "#elapsed time: " << T.time()*1.0e9 << " [nano s] \n";
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

