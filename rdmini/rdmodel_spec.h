#ifndef RDMODEL_SPEC_H
#define RDMODEL_SPEC_H

/** Code representation of YAML model specification.
 *
 * RX schema: rdmodel_spec.rx
*/

struct species_spec {
    std::string name;
    double diffusivity;

    // initial concentrations or counts specified by region
    std::vector<std::pair<std::string,double>> concentration;
    std::vector<std::pair<std::string,double>> counts;
};

struct reaction_spec {
    std::string name;
    std::vector<std::string> left;
    std::vector<std::string> right;

    // rates specified by region
    std::vector<std::pair<std::string,double>> rate;
};

struct region_box_spec {
    std::vector<double> min_point;
    std::vector<double> max_point;
};

struct region_interval_spec {
    std::string axis;
    double min,max;
};

struct region_spec {
    std::string name;
    enum region_type_enum { box, interval } region_type; // add more as required

    union {
        region_box_spec box_spec;
        region_interval_spec interval_spec;
    };
};

struct cells_spec {
    

};


struct rdmodel_spec {
    std::string model;

    std::vector<species_spec> species;
    std::vector<reaction_spec> reaction;
    std::vector<region_spec> regions;
    std::vector<cells_spec> cells;
};

#endif // ndef RDMODEL_SPEC_H
