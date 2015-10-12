#ifndef RDMODEL_H_
#define RDMODEL_H_

#include <iosfwd>
#include <string>
#include <set>
#include <stdexcept>

#include "rdmini/named_collection.h"

struct model_io_error: std::runtime_error {
    model_io_error(const std::string &what_arg): std::runtime_error(what_arg) {}
};

struct species_info {
    std::string name;
    double diffusivity;
};

struct reaction_info {
    std::string name;
    std::multiset<int> left,right;
    double rate;
};

struct rd_model {
    std::string name;
    named_collection<species_info> species;
    named_collection<reaction_info> reactions;

    void clear() {
        species.clear();
        reactions.clear();
    }

    friend std::ostream &operator<<(std::ostream &O,const rd_model &M);
};

rd_model rd_model_read(std::istream &,const std::string &model_name="");
rd_model rd_model_read(const std::string &,const std::string &model_name="");

#endif // ndef RDMODEL_H_

