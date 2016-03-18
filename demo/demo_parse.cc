#include <string>
#include <cstring>
#include <cstddef>
#include <iostream>
#include <fstream>

#include "rdmini/rdmodel.h"

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

struct cl_args {
    std::string model_file;
    std::string model_name;
};

cl_args parse_cl_args(int argc,char **argv) {
    cl_args A;

    enum parse_state_enum { no_opt, opt_m } parse_state = no_opt;
    bool has_opt_m=false;
    bool has_file=false;

    int i=0;
    while (++i<argc) {
        const char *arg=argv[i];
        switch (parse_state) {
        case no_opt:
            if (!strcmp(arg,"-m"))
                parse_state=opt_m;
            else if (arg[0]=='-' && arg[1])
                throw usage_error("unrecognized option "+std::string(arg));
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
        }
    }

    if (parse_state!=no_opt)
        throw usage_error("missing option argument");

    return A;
}

int main(int argc, char **argv) {
    const char *basename=strrchr(argv[0],'/');
    basename=basename?basename+1:argv[0];
    int rc=0;

    try {
        cl_args A=parse_cl_args(argc,argv);

        rdmini::rd_model M;

        if (A.model_file.empty() || A.model_file=="-")
            M=rdmini::rd_model_read(std::cin,A.model_name);
        else {
            std::ifstream file(A.model_file);
            if (!file) throw fatal_error("unable to open file for reading");

            M=rdmini::rd_model_read(file,A.model_name);
        }

        std::cout << M << "\n";
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

