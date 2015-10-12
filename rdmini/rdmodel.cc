#include <iostream>
#include <iomanip>
#include <string>

#include "rdmini/rdmodel.h"
#include "rdmini/yamlview.h"

static std::ostream &emit_reaction_expr(std::ostream &O,const rd_model &M,const std::multiset<int> &x) {
    bool first=true;
    if (x.empty()) return O << "Ø";

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
    O << "species:\n";
    for (auto &s: M.species) {
        O << " " << std::setw(10) << std::right << (s.name+":") << " ";
        O << "diffusivity=" << std::setw(10) << std::left << s.diffusivity << "\n";
    }
    O << "reactions:\n";
    for (auto &r: M.reactions) {
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


static void parse_geometry(rd_model &M,const yaml_node_view &G) {
    // two cases: simple wmvol, or box.
    yaml_node_view params=G["wmvol"];
    if (params) {
        // wmvol

        return;
    }

    params=G["box"];
    if (params) {
        // box

        return;
    }

    throw model_io_error("unknown geometry specification: "+G.where());
}

static void parse_species(rd_model &M,const yaml_node_view &S) {
    try {
        double diff_value=0;
        yaml_node_view diff=S["diffusivity"];
        if (diff) diff_value=std::stod(diff.str());

        species_info species={S["name"].str(),diff_value};
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

static void parse_reaction(rd_model &M,const yaml_node_view &R) {
    std::string name;

    yaml_node_view name_node=R["name"];
    if (name_node) {
        name=name_node.str();
        if (M.reactions.index(name)>=0)
            throw model_io_error("reaction "+name+" already specified: "+R.where());
    }
    else {
        name=M.reactions.unique_key("_r");
    }

    yaml_node_view rate_node=R["rate"];
    double rate=0,rate_rev=0;
    if (!rate_node || rate_node.is_map() || rate_node.size()<1 || rate_node.size()>2)
        throw model_io_error("unknown reaction rate specification: "+R.where());
    
    rate=std::stod(rate_node[0].str());
    std::multiset<int> left=parse_species_list(M,R["left"]);
    std::multiset<int> right=parse_species_list(M,R["right"]);

    reaction_info reaction={name,left,right,rate};
    M.reactions.insert(reaction);

    if (rate_node.size()>1) {
        std::string name_rev=M.reactions.unique_key(name+"_rev");
        rate_rev=std::stod(rate_node[1].str());

        reaction_info reaction={name_rev,right,left,rate_rev};
        M.reactions.insert(reaction);
    }
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

    yaml_node_view geom=root["geometry"];
    if (!geom) throw model_io_error("missing geometry specification");
    parse_geometry(M,geom);

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

