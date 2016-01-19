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

struct model_incorrectBiologicalValue_error: std::runtime_error {
    model_incorrectBiologicalValue_error(const std::string &what_arg): std::runtime_error(what_arg) {}
};


/// Use this as functor to throw on invalid models
/// This implementations is used as a work around to avoid adding getters/setters to defined struct 
/// Could be made more generic to throw any type of error from boolean and type of test
struct throwInvalidModl {
    void operator()(bool b) { throw model_incorrectBiologicalValue_error("Negative biological value"); } 

};

struct species_info { 
    std::string name;
    double diffusivity;
    double concentration;

    bool isModlValid() {
        if ((diffusivity < 0.0f) || (concentration < 0.0f))
            return false;
        return true;
   }
};

struct reaction_info {
    std::string name;
    std::multiset<int> left,right;
    double rate;

    bool isModlValid() {
        if (rate < 0.0f) 
            return false;
        return true;
    }
};

struct cell_info {
    struct neighbour_data {
        size_t cell_id;
        double diff_coef;

        explicit neighbour_data(size_t cell_id_=0, double diff_coef_=0): cell_id(cell_id_), diff_coef(diff_coef_) {}
    };

    double volume;
    std::vector<neighbour_data> neighbours;

    bool isModlValid() {
        if (volume<0.0f)
            return false;
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

