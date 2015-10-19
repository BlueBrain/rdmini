#include <string>
#include <cstring>
#include <cstddef>
#include <iostream>
#include <fstream>

#include "rdmini/rdmodel.h"
#include "rdmini/serial_ssa.h"

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
    "  -m MODEL    Load the model named MODEL\n";
    "  -n N        Run simulation N steps\n";
    "  -t TIME     Run simulation for TIME simulated seconds\n";
    "  -d TIME     Sample simulation every TIME seconds\n";
    "\nOne of -n or -t must be specified.\n";

struct cl_args {
    std::string model_file;
    std::string model_name;
    double dt=0;
    double t_end=0;
    size_t n_events=0;
};

cl_args parse_cl_args(int argc,char **argv) {
    cl_args A;

    enum parse_state_enum { no_opt, opt_m, opt_n, opt_t, opt_dt } parse_state = no_opt;
    bool has_opt_m=false;
    bool has_opt_n=false;
    bool has_opt_t=false;
    bool has_opt_dt=false;
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
                default:
                    throw usage_error("unrecognized option "+std::string(arg));
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
            A.dt=std::stod(arg);
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
    explicit emit_sim(const rd_model &M): n_species(M.n_species()), n_cells(M.n_cells) {
        // prepare csv-style header
        std::stringstream s;
        s << "time,cell";

        for (size_t i=0; i<n_species; ++i) 
            s << ',' << M.species[i].name;
        
        s << '\n';
        header=s.str();
    }
    
    std::ostream &emit_header(std::ostream &O) { return O << header; }

    template <typename Sim>
    std::ostream &emit_state(std::ostream &O, double t, const Sim &sim) {
        O << t;
        for (size_t i=0; i<n_cell; ++i) {
            O << ',' << i;
            for (size_t j=0; j<n_species; ++j)
                O << ',' << sim.count(j,i) ;
            O << '\n';
        }
        return O;
    }

    size_t n_species,n_cells;
    std::string header;
};


int main(int argc, char **argv) {
    const char *basename=strrchr(argv[0],'/');
    basename=basename?basename+1:argv[0];
    int rc=0;

    try {
        cl_args A=parse_cl_args(argc,argv);

        rd_model M;

        if (A.model_file.empty() || A.model_file=="-")
            M=rd_model_read(std::cin,A.model_name);
        else {
            std::ifstream file(A.model_file);
            if (!file) throw fatal_error("unable to open file for reading");

            M=rd_model_read(file,A.model_name);
        }


        emit_sim emitter(M);
        emitter.emit_header(std::cout,M);

        std::minstd_rand g;
        serial_ssa S(M,0);
        emit_state(0,S);
        if (A.n_events>0) { // event-by-event simulation
            for (size_t n=0; n<A.n_events; ++n) {
                double t=S.advance(g);
                emitter.emit_state(std::cout,t,S);
            }
        }
        else {
            double t=0;
            while (t<A.t_end) {
                t=S.advance(t+A.dt,g);
                emitter.emit_state(std::cout,t,S);
            }
        }
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

