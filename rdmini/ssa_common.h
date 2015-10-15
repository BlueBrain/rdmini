#ifndef SSA_COMMON_H_
#define SSA_COMMON_H_

#include <stdexcept>

struct ssa_error: std::runtime_error {
    ssa_error(const std::string &m): std::runtime_error(m) {}
    ssa_error(const char *m): std::runtime_error(m) {}
};

#endif // ndef SSA_COMMON_H_
