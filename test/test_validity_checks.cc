#include <cassert>
#include <gtest/gtest.h>

#include "rdmini/check_valid.h"

using namespace rdmini;

struct dummy_class: check_valid<dummy_class> {
    dummy_class(): valid_flag(true) {}

    void bad_method() { valid_flag=false; }

    bool is_valid() const { return valid_flag; }
    void assert_valid() const { assert(is_valid()); }

private:
    bool valid_flag;
};
    
TEST(check_valid,check_valid) {


}
