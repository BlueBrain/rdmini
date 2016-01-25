#ifndef RDMODEL_H_
#define RDMODEL_H_

#include <iosfwd>
#include <string>
#include <set>
#include <stdexcept>

#include "rdmini/named_collection.h"
#include "rdmini/check_valid.h"

struct model_io_error: std::runtime_error {
    model_io_error(const std::string &what_arg): std::runtime_error(what_arg) {}
};

struct model_invalid_error: std::runtime_error {
    model_invalid_error(const std::string &what_arg): std::runtime_error(what_arg) {}
};

struct species_info: rdmini::check_valid_api<species_info> { 
    std::string name;
    double diffusivity=0;
    double concentration=0;

    species_info() =default;
    species_info(const std::string &name_,double diff_,double conc_):
        name(name_),diffusivity(diff_),concentration(conc_) {}

    rdmini::valid_info is_valid() const {
        if (diffusivity<0) return "negative diffusivity";
        if (concentration<0) return "negative concentration";
        return true;
    }
};

struct reaction_info: rdmini::check_valid_api<reaction_info> {
    std::string name;
    std::multiset<int> left,right;
    double rate=0;

    reaction_info() =default;
    reaction_info(const std::string &name_,const std::multiset<int> &left_,const std::multiset<int> &right_,double rate_):
        name(name_),left(left_),right(right_),rate(rate_) {}

    rdmini::valid_info is_valid() const {
        if (rate<0) return "negative reaction rate constant";
        return true;
    }
};

struct cell_info: rdmini::check_valid_api<cell_info> {
    struct neighbour_data {
        size_t cell_id;
        double diff_coef;

        explicit neighbour_data(size_t cell_id_=0, double diff_coef_=0): cell_id(cell_id_), diff_coef(diff_coef_) {}
    };

    double volume=0;
    std::vector<neighbour_data> neighbours;

    rdmini::valid_info is_valid() const {
        if (volume<=0) return "non-positive cell volume";
        return true;
    }
};

struct cell_set {
    std::string name;
    std::vector<size_t> cells;
};

struct rd_model {
    std::string name;
    named_collection<species_info> species;
    named_collection<reaction_info> reactions;
    named_collection<cell_set> cell_sets;
    std::vector<cell_info> cells;

    void clear() {
        species.clear();
        reactions.clear();
    }

    friend std::ostream &operator<<(std::ostream &O,const rd_model &M);

    size_t n_species() const { return species.size(); }
    size_t n_reactions() const { return reactions.size(); }
    size_t n_cells() const { return cells.size(); }
};

rd_model rd_model_read(std::istream &,const std::string &model_name="");
rd_model rd_model_read(const std::string &,const std::string &model_name="");

#endif // ndef RDMODEL_H_

