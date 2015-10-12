#ifndef YAMLVIEW_H_
#define YAMLVIEW_H_

/** Simple partial C++ wrapper for libyaml parser for safer resource management */

#include <iostream>
#include <memory>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>

extern "C" {
#include <yaml.h>
}

struct yaml_error: std::runtime_error {
    yaml_error(const std::string &what_str,const std::string &where_str_=""):
        std::runtime_error(what_str),where_str(where_str_) {}

    const std::string &where() const { return where_str; }
private:
    std::string where_str;
};

struct yaml_parser;
struct yaml_doucment;

struct yaml_node_view {
    yaml_node_view(): n(null_node), n_bis(null_node) {
        tag_std=TAG_NULL;
        sz=0;
    }

    yaml_node_view(const std::shared_ptr<yaml_document_t> &doc_,yaml_node_t n_,yaml_node_t n_bis_=null_node):
        doc(doc_),n(n_),n_bis(n_bis_)
    {
        tag_std=tag_str_to_enum(n.tag);
        switch (n.type) {
        case YAML_SCALAR_NODE:
            sz=1;
            break;
        case YAML_SEQUENCE_NODE:
            sz=(size_t)(n.data.sequence.items.top-n.data.sequence.items.start);
            break;
        case YAML_MAPPING_NODE:
            sz=(size_t)(n.data.mapping.pairs.top-n.data.mapping.pairs.start);
            break;
        default:
            sz=0;
        }
    }

    enum std_tag_enum {
        TAG_NULL, TAG_BOOL, TAG_STR, TAG_INT, TAG_FLOAT, TAG_TIMESTAMP, TAG_SEQ, TAG_MAP, TAG_OTHER
    };

    bool has_standard_tag() const { return tag_std!=TAG_OTHER; }
    std_tag_enum standard_tag() const { return tag_std; }

    std::string tag() const { return reinterpret_cast<char *>(n.tag); }

    bool is_scalar() const { return n.type==YAML_SCALAR_NODE; }
    bool is_seq() const { return n.type==YAML_SEQUENCE_NODE; }
    bool is_map() const { return n.type==YAML_MAPPING_NODE; }

    operator bool() const { return n.type!=YAML_NO_NODE; }
    
    size_t size() const { return sz; }

    /* Returns ith node or node pair in sequence or mapping */
    yaml_node_view operator[](int i) const {
        if (i<0 || i>=size()) throw yaml_error("index out of range",where());

        switch (n.type) {
        case YAML_SCALAR_NODE:
            return *this;
        case YAML_SEQUENCE_NODE:
            return yaml_node_view(doc,*yaml_document_get_node(doc.get(),n.data.sequence.items.start[i]));
        case YAML_MAPPING_NODE:
            return yaml_node_view(doc,
                    *yaml_document_get_node(doc.get(),n.data.mapping.pairs.start[i].key),
                    *yaml_document_get_node(doc.get(),n.data.mapping.pairs.start[i].value));
        default:
            return yaml_node_view();
        }
    }

    /* Returns first value node in mapping with key matching k */
    yaml_node_view operator[](const std::string &k) const {
        if (n.type!=YAML_MAPPING_NODE) throw yaml_error("not a mapping node",where());

        int n_entries=(int)size();
        for (int i=0;i<n_entries;++i) {
            yaml_node_view key=(*this)[i];
            if (key.is_scalar() && key.str()==k)
                return yaml_node_view(doc,*yaml_document_get_node(doc.get(),n.data.mapping.pairs.start[i].value));
        }

        return yaml_node_view();
    }

    yaml_node_view operator[](const char *k) const { return (*this)[std::string(k)]; }
        
    std::string str() const {
        if (n.type!=YAML_SCALAR_NODE) throw yaml_error("not a scalar node",where());
        return reinterpret_cast<char *>(n.data.scalar.value);
    }

    bool operator==(const std::string &text) const {
        return n.type==YAML_SCALAR_NODE && text==reinterpret_cast<char *>(n.data.scalar.value);
    }

    bool operator!=(const std::string &text) const { return !(*this==text); }

    // return value of (key,value) pair from a mapping entry
    yaml_node_view value() const { return yaml_node_view(doc,n_bis); }

    std::string where() const {
        if (!*this) return "";
        else return "line "+std::to_string(n.start_mark.line)+" column "+std::to_string(n.start_mark.column);
    }

private:
    static std_tag_enum tag_str_to_enum(yaml_char_t *tag_) {
        const char *tag=reinterpret_cast<char *>(tag_);
        if (!strcmp(tag,YAML_NULL_TAG)) return TAG_NULL;
        if (!strcmp(tag,YAML_BOOL_TAG)) return TAG_BOOL;
        if (!strcmp(tag,YAML_STR_TAG)) return TAG_STR;
        if (!strcmp(tag,YAML_INT_TAG)) return TAG_INT;
        if (!strcmp(tag,YAML_FLOAT_TAG)) return TAG_FLOAT;
        if (!strcmp(tag,YAML_TIMESTAMP_TAG)) return TAG_TIMESTAMP;
        if (!strcmp(tag,YAML_SEQ_TAG)) return TAG_SEQ;
        if (!strcmp(tag,YAML_MAP_TAG)) return TAG_MAP;
        return TAG_OTHER;
    }

    static yaml_node_t null_node;

    std::shared_ptr<yaml_document_t> doc;
    yaml_node_t n,n_bis;

    size_t sz;
    std_tag_enum tag_std;
};

yaml_node_t yaml_node_view::null_node = {YAML_NO_NODE, (yaml_char_t *)YAML_NULL_TAG}; 

struct yaml_document {
    struct yaml_document_deleter {
        void operator()(yaml_document_t *p) {
            yaml_document_delete(p);
            delete p;
        }
    };

    std::shared_ptr<yaml_document_t> doc;

    yaml_node_view root() {
        return yaml_node_view(doc,*yaml_document_get_root_node(doc.get()));
    }

    yaml_node_view operator[](int index) {
        return yaml_node_view(doc,*yaml_document_get_node(doc.get(),index));
    }

    yaml_document(): doc(new yaml_document_t,yaml_document_deleter()) {
        yaml_document_initialize(doc.get(),nullptr,nullptr,nullptr,1,1);
    }
    
    yaml_document(yaml_parser &p);

    operator bool() const {
        return yaml_document_get_root_node(doc.get())!=nullptr;
    }

};

inline int yaml_istream_reader(void *data,unsigned char *buffer,size_t size,size_t *size_read) {
    std::istream &I=*static_cast<std::istream *>(data);

    I.read(reinterpret_cast<char *>(buffer),size);
    *size_read=I.gcount();
    return !I.bad();
}

struct yaml_parser {
    struct yaml_parser_deleter {
        void operator()(yaml_parser_t *p) {
            yaml_parser_delete(p);
            delete p;
        }
    };

    std::shared_ptr<yaml_parser_t> P;

    explicit yaml_parser(FILE *source): P(new yaml_parser_t,yaml_parser_deleter()) {
        yaml_parser_initialize(P.get());
        yaml_parser_set_input_file(P.get(),source);
    }
        
    explicit yaml_parser(const std::string &source): P(new yaml_parser_t,yaml_parser_deleter()) {
        yaml_parser_initialize(P.get());
        yaml_parser_set_input_string(P.get(),reinterpret_cast<const unsigned char *>(source.c_str()),source.size());
    }

    explicit yaml_parser(std::istream &I): P(new yaml_parser_t,yaml_parser_deleter()) {
        yaml_parser_initialize(P.get());
        yaml_parser_set_input(P.get(),yaml_istream_reader,(void *)&I);
    }

    yaml_document next_document() { return yaml_document(*this); }
};

yaml_document::yaml_document(yaml_parser &p): doc(new yaml_document_t,yaml_document_deleter()) {
    if (!yaml_parser_load(p.P.get(),doc.get())) {
        // create empty document
        yaml_document_initialize(doc.get(),nullptr,nullptr,nullptr,1,1);
    }
}

#endif // ndef YAMLVIEW_H_

