#include "yamlview.h"

namespace rdmini {

yaml_node_t yaml_node_view::null_node = {YAML_NO_NODE, (yaml_char_t *)YAML_NULL_TAG}; 

yaml_document::yaml_document(yaml_parser &p): doc(new yaml_document_t,yaml_document_deleter()) {
    if (!yaml_parser_load(p.P.get(),doc.get())) {
        // create empty document
        yaml_document_initialize(doc.get(),nullptr,nullptr,nullptr,1,1);
    }
}

int yaml_istream_reader(void *data,unsigned char *buffer,size_t size,size_t *size_read) {
    std::istream &I=*static_cast<std::istream *>(data);

    I.read(reinterpret_cast<char *>(buffer),size);
    *size_read=I.gcount();
    return !I.bad();
}

} // namespace rdmini
