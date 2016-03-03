#ifndef RDMINI_EXCEPTIONS_H_
#define RDMINI_EXCEPTIONS_H_

#include <stdexcept>

namespace rdmini {

/** Thrown when an implementation does not support an optional interface. */

struct operation_not_supported: std::logic_error {
    operation_not_supported(const std::string &m): std::logic_error(m) {}
    operation_not_supported(const char *m): std::logic_error(m) {}
};


/** Thrown when a domain or range check on a parameter fails. */

struct invalid_value: std::runtime_error {
    invalid_value(const std::string &m): std::runtime_error(m) {}
    invalid_value(const char *m): std::runtime_error(m) {}
};

/** Represents internal error in an SSA implementation */

struct ssa_error: std::runtime_error {
    ssa_error(const std::string &m): std::runtime_error(m) {}
    ssa_error(const char *m): std::runtime_error(m) {}
};

} // namespace rdmini

#endif // ndef RDMINI_EXCEPTIONS_H_
