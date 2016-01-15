#include <iostream>
#include <iomanip>
#include <string>
//#include <assert.h>
#include <limits>

#include "rdmini/rdmodel.h"
#include "rdmini/yamlview.h"
#include "rdmini/point.h"
#include "rdmini/range_seq.h"

static std::ostream &emit_reaction_expr(std::ostream &O,const rd_model &M,const std::multiset<int> &x) {
    bool first=true;
    if (x.empty()) return O << u8"\u00d8"; // crayCC barfs on utf8 in source

    for (auto i=x.begin();i!=x.end();i=x.upper_bound(*i)) {
        if (!first) O << " + ";
        first=false;

        int c=x.count(*i);
        if (c>1) O << c;
        O << M.species[*i].name;
    }
    return O;
}

std::ostream &operator<<(std::ostream &O,const rd_model &M) {
    O << "cells:\n";
    for (const auto &c: M.cell_sets) {
        O << " " << std::setw(10) << std::right << (c.name+":") << " "
          << range_seq<size_t>(c.cells.begin(),c.cells.end()) << "\n";
    }

    O << "species:\n";
    for (const auto &s: M.species) {
        O << " " << std::setw(10) << std::right << (s.name+":") << " ";
        O << "diffusivity=" << std::setw(10) << std::left << s.diffusivity << "\n";
    }
    O << "reactions:\n";
    for (const auto &r: M.reactions) {
        O << " " << std::setw(10) << std::right << (r.name+":") << " ";
        O << "rate=" << std::setw(10) << std::left << r.rate << "\t";
        emit_reaction_expr(O,M,r.left);
        O << " -> ";
        emit_reaction_expr(O,M,r.right);
        O << "\n";
    }
    return O;
}

static rd_model rd_model_read_yaml(yaml_parser,const std::string &);

rd_model rd_model_read(std::istream &I,const std::string &model_name) {
    return rd_model_read_yaml(yaml_parser(I),model_name);
}

rd_model rd_model_read(const std::string &s,const std::string &model_name) {
    return rd_model_read_yaml(yaml_parser(s),model_name);
}

template <typename NamedCollection>
static std::string check_or_make_unique_name(const NamedCollection &C,const yaml_node_view &name_node,const std::string &fallback) {
    if (name_node) {
        std::string name=name_node.str();
        if (C.index(name)>=0) 
            throw model_io_error(name+" already specified before use at: "+name_node.where());
        return name;
    }
    else
        return C.unique_key(fallback);
}

// Parse species info

static void parse_species(rd_model &M,const yaml_node_view &S) {
    std::string name=check_or_make_unique_name(M.species,S["name"],"_s");

    try {
        double diff_value=0;
        yaml_node_view diff=S["diffusivity"];
        if (diff) diff_value=std::stod(diff.str());

        double conc_value=0;
        yaml_node_view conc=S["concentration"];
        if (conc) conc_value=std::stod(conc.str());

        species_info species={name,diff_value,conc_value};
        if (species.diffusivity<std::numeric_limits<double>::epsilon()) 
          throw model_incorrectBiologicalValue_error("Value of diffusivity is negative");
        if (species.concentration<std::numeric_limits<double>::epsilon()) 
          throw model_incorrectBiologicalValue_error("Value of concentration is negative");
        //assert(species.diffusivity<std::numeric_limits<double>::epsilon());
        //assert(species.concentration<std::numeric_limits<double>::epsilon());

        M.species.insert(species);
    }
    catch (yaml_error &error) {
        throw model_io_error("parsing species failue: "+error.where());
    }
}

static std::multiset<int> parse_species_list(const rd_model &M,const yaml_node_view &species) {
    std::multiset<int> s;

    if (species.is_map()) goto error;

    for (int i=0;i<species.size();++i) {
        yaml_node_view sitem=species[i];
        if (!sitem.is_scalar()) goto error;

        int j=M.species.index(species[i].str());
        if (j<0) goto error;

        s.insert(j);
    }
    return s;
    
error:
    throw model_io_error("improper species list in reaction specification: "+species.where());
}

// Parse reaction info

static void parse_reaction(rd_model &M,const yaml_node_view &R) {
    std::string name=check_or_make_unique_name(M.reactions,R["name"],"_r");

    yaml_node_view rate_node=R["rate"];
    double rate=0,rate_rev=0;
    if (!rate_node || rate_node.is_map() || rate_node.size()<1 || rate_node.size()>2)
        throw model_io_error("unknown reaction rate specification: "+R.where());
    
    rate=std::stod(rate_node[0].str());
    std::multiset<int> left=parse_species_list(M,R["left"]);
    std::multiset<int> right=parse_species_list(M,R["right"]);

    reaction_info reaction={name,left,right,rate};
    //assert(reaction.rate<std::numeric_limits<double>::epsilon());
    if (reaction.rate<std::numeric_limits<double>::epsilon()) 
          throw model_incorrectBiologicalValue_error("Value of reaction rate is negative");

    M.reactions.insert(reaction);

    if (rate_node.size()>1) {
        std::string name_rev=M.reactions.unique_key(name+"_rev");
        rate_rev=std::stod(rate_node[1].str());

        reaction_info reaction={name_rev,right,left,rate_rev};
        M.reactions.insert(reaction);
    }
}

// Parse cell info

static void parse_cells_selection(rd_model &M,const yaml_node_view &e) {
    throw model_io_error("cell selections not supported yet!");
}

static point3d parse_point(const yaml_node_view &e) {
    return point3d{std::stod(e[0].str()), std::stod(e[1].str()), std::stod(e[2].str())};
}


static void parse_cells_wmvol(rd_model &M,const yaml_node_view &e) {
    try {
        std::string name=check_or_make_unique_name(M.cell_sets,e["name"],"_wmvol");
        size_t c0=M.cells.size();

        cell_info ci={std::stod(e["volume"].str())};
        if (ci.volume<std::numeric_limits<double>::epsilon()) 
          throw model_incorrectBiologicalValue_error("Value of volume is negative");
//        assert(ci.volume<std::numeric_limits<double>::epsilon());
        M.cells.push_back(ci);

        cell_set cs={name,{c0}};
        M.cell_sets.insert(cs);
    }
    catch (yaml_error &) {
        throw model_io_error("error parsing cells wmvol entry: "+e.where());
    }
}

static void parse_cells_grid(rd_model &M,const yaml_node_view &e) {
    try {
        std::string name=check_or_make_unique_name(M.cell_sets,e["name"],"_grid");

        double scale=1;
        if (auto scale_node=e["scale"]) scale=std::stod(scale_node.str());

        auto extent_lb=parse_point(e["extent"][0])*scale;
        auto extent_ub=parse_point(e["extent"][1])*scale;
        point3d d=extent_ub-extent_lb;
        
        size_t n[3];
        for (unsigned i=0; i<3; ++i) n[i]=std::stoull(e["counts"][i].str());

        d={d[0]/n[0], d[1]/n[1], d[2]/n[2] };
        double vol=d[0]*d[1]*d[2];

        size_t c0=M.cells.size();
        size_t nc=n[0]*n[1]*n[2];

        auto cidx=[=](size_t i, size_t j, size_t k) { return c0+i+j*n[0]+k*n[0]*n[1]; };

        point3d dc=1.0/(d*d);
        
        for (size_t k=0; k<n[2]; ++k) {
            for (size_t j=0; j<n[1]; ++j) {
                for (size_t i=0; i<n[0]; ++i) {
                    cell_info ci={vol};
                    if (i>0)      ci.neighbours.emplace_back(cidx(i-1,j,k),dc[0]);
                    if (i<n[0]-1) ci.neighbours.emplace_back(cidx(i+1,j,k),dc[0]);
                    if (j>0)      ci.neighbours.emplace_back(cidx(i,j-1,k),dc[1]);
                    if (j<n[1]-1) ci.neighbours.emplace_back(cidx(i,j+1,k),dc[1]);
                    if (k>0)      ci.neighbours.emplace_back(cidx(i,j,k-1),dc[2]);
                    if (k<n[2]-1) ci.neighbours.emplace_back(cidx(i,j,k+1),dc[2]);

                    M.cells.push_back(ci);
                }
            }
        }

        cell_set cs={name};
        for (size_t c=c0; c<c0+nc; ++c) cs.cells.push_back(c);
        M.cell_sets.insert(cs);
    }
    catch (yaml_error &E) {
        throw model_io_error("error parsing cells grid entry: "+e.where()+": "+E.what());
    }
}

static void parse_cells(rd_model &M,const yaml_node_view &R) {
    std::vector<yaml_node_view> selection_nodes;

    for (int i=0; i<R.size(); ++i) {
        // node is either a (named) collection of cells, or a a view (subset)
        // for now, let's not allow views on views to simplify parsing.
        yaml_node_view e=R[i];
        if (e=="select")
            selection_nodes.push_back(e.value());
        else if (e=="wmvol")
            parse_cells_wmvol(M,e.value());
        else if (e=="grid")
            parse_cells_grid(M,e.value());
        else
            throw model_io_error("unrecognised entry in cells specification: "+e.where());
    }
        
    for (auto e: selection_nodes) parse_cells_selection(M,e);
}

rd_model rd_model_read_yaml(yaml_parser Y,const std::string &model_name) {
    yaml_node_view root,model_node;
    do {
        yaml_document D=Y.next_document();
        if (!D) throw model_io_error("model specification not found");

        root=D.root();
        model_node=root["model"];
    } while (!model_node || !model_node.is_scalar() || !model_name.empty() && model_name!=model_node.str());

    rd_model M;
    M.name=model_node.str();

    yaml_node_view cells=root["cells"];
    if (!cells) throw model_io_error("missing cells specification");
    parse_cells(M,cells);

    // add all species
    for (int i=0;i<root.size();++i) {
        yaml_node_view e=root[i];
        if (e!="species") continue;

        parse_species(M,e.value());
    }

    // add all reactions
    for (int i=0;i<root.size();++i) {
        yaml_node_view e=root[i];
        if (e!="reaction") continue;

        parse_reaction(M,e.value());
    }
    
    return M;
}

