/* Representation of a simple reaction-diffusion system over a mesh of well-mixed volumes. */

#include "point.h"

/* wmvols */

typedef unsigned int wmvol_id;

struct wmvol_info {
    bbox3d extent;
    double volume=0;
};

struct wmvol_adj_info {
    int idx=-1;
    double diff_coeff=0.0;
};

struct wmvol_info {
    wmvol_info(size_t n_wmvol_,size_t max_adj)) {
        n_wmvol=n_wmvol_;
        max_adj=max_adj_;

        extents.resize(n_wmvol);
        adjancencies.resize(n_wmvol*max_adj);
    }

    wmvol_extent &extent(int i) { return extents[i]; }
    const wmvol_extent &extent(int i) const { return extents[i]; }

    size_t size() { return extents.size(); }

    unsigned int n_adj(size_t i) {
        int n=0;
        int offset=i*max_adj;
        for (int j=0;j<max_adj;++j) n+=(adjacencies[offset+j].idx>=0);
        return n;
    }

    wmvol_adj_info *operator[](size_t i) { return adjacencies[i*max_adj]; }
    const wmvol_adj_info *operator[](size_t i) const { return adjacencies[i*max_adj]; }

    struct wmvol_extent {
        bbox3d bounds;
        double volume=0;
    };
    std::vector<wmvol_extent> extents;
    std::vector<wmvol_adj_info> adjacencies;
};

wmvol_info wmvols;

/* populations:
 *
 * A population represents the (count of) a species in a well-mixed volume.
 *
 * Two maps:
 *     pop_idx -> (wmvol_idx, species_idx)
 *     (wmvol_idx, species_idx) -> pop_idx or -1
 *
 * Currently, every species can be in every wmvol, so the maps are trivial.
 */

struct pop_info {
    pop_info(size n_wmvol_,size_t n_species_):
        n_wmvol(n_wmvol_),n_species(n_species_),count(n_wmvol_*n_species_) {}

    size_t n_wmvol,n_species;
    std::vector<size_t> initial_count;

    int species_idx(size_t p) const {
        return p%n_wmvol;
    }

    int wmvol_idx(size_t p) const {
        return p/n_wmvol;
    }

 // ...?
};

pop_info pops;
    

/* event interface: sim consumer receives a sequence of sim events
 * and processes them for analysis, display, etc. */

struct sim_event {
    enum ev_type_enum { SIMEV_REAC, SIMEV_POP_DELTA

};


/* externally visable total system state */

struct {
    double t; // (simulation) time
    
    std::vector<size_t> pop_count;

};

