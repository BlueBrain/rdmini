#include <cassert>
#include <utility>
#include <gtest/gtest.h>

#include "rdmini/check_valid.h"

struct dummy_class: rdmini::check_valid_api<dummy_class> {
    void bad_method() { valid_flag=false; }
    bool is_valid() const { return valid_flag; }

private:
    bool valid_flag=true;
};
    
TEST(check_valid,check_valid) {
    dummy_class x;

    ASSERT_TRUE(x.is_valid());
    ASSERT_NO_THROW(x.check_valid());

    x.bad_method();

    ASSERT_FALSE(x.is_valid());
    ASSERT_THROW(x.check_valid(),rdmini::validation_failure);

    const auto &const_x=x;
    ASSERT_FALSE(const_x.is_valid());
    ASSERT_THROW(const_x.check_valid(),rdmini::validation_failure);
}

TEST(check_valid,check_valid_user_message) {
    dummy_class x;
    x.bad_method();
    EXPECT_TRUE(false);

    try {
        const auto &const_x=x;
        const_x.check_valid("foobar");
        FAIL() << "did not throw exception";
    }
    catch (rdmini::validation_failure &ex) {
        ASSERT_EQ(std::string("foobar"),ex.what());
    }
    catch (...) {
        FAIL() << "did not throw rdmini::validation_failure";
    }
}

TEST(check_valid,check_valid_ex) {
    dummy_class x;
    x.bad_method();

    using funny_exception=std::tuple<const char *,double>;

    try {
        const auto &const_x=x;
        const_x.check_valid_ex<funny_exception>("quux",17.0);
        FAIL() << "did not throw exception";
    }
    catch (funny_exception &ex) {
        ASSERT_STREQ("quux",std::get<0>(ex));
        ASSERT_EQ(17.0,std::get<1>(ex));
    }
    catch (...) {
        FAIL() << "did not throw custom exception";
    }
}

// use a perfect forwarding wrapper to ensure conversions are
// implicit...

template <typename X>
rdmini::valid_info return_as_valid_info_type(X &&x) {
    return rdmini::valid_info(std::forward<X>(x));
}

TEST(check_valid,valid_info) {
    // default ctor -> false
    ASSERT_FALSE(rdmini::valid_info());

    // construction from bool
    ASSERT_FALSE(return_as_valid_info_type(false));
    ASSERT_TRUE(return_as_valid_info_type(true));

    // construction from C string
    ASSERT_FALSE(return_as_valid_info_type("foobar"));
    ASSERT_EQ(std::string("foobar"),return_as_valid_info_type("foobar").what());
}

struct zero_n: rdmini::check_valid_api<zero_n> {
    int n=0;

    rdmini::valid_info is_valid() const {
        if (n<0) return "n is negative";
        if (n>0) return "n is positive";
        return true;
    }
};

TEST(check_valid,check_valid_what) {
    zero_n z;

    ASSERT_TRUE(z.is_valid());
    ASSERT_NO_THROW(z.check_valid());

    z.n=3;
    ASSERT_FALSE(z.is_valid());

    try {
        z.check_valid();
        FAIL() << "did not throw exception";
    }
    catch (rdmini::validation_failure &ex) {
        ASSERT_EQ(std::string("n is positive"),ex.what());
    }
    catch (...) {
        FAIL() << "did not throw rdmini::validation_failure";
    }
}

TEST(check_valid,assert_valid_ok) {
    dummy_class x;

    ASSERT_EXIT({x.assert_valid(); ::exit(0);},::testing::ExitedWithCode(0),"");

    const auto &const_x=x;
    ASSERT_EXIT({const_x.assert_valid(); ::exit(0);},::testing::ExitedWithCode(0),"");
}

TEST(check_valid,assert_valid_not_ok) {
    dummy_class x;
    x.bad_method();
    ASSERT_DEATH(x.assert_valid(),"validation failure");
}

TEST(check_valid,assert_valid_message) {
    zero_n z;
    z.n=-10;

    ASSERT_DEATH(z.assert_valid(),"validation failure: n is negative");
}

TEST(check_valid,check_valid_guard_ref) {
    auto set_zero_n=[](zero_n &z,int value) {
        auto _(check_valid_guard(z));
        z.n=value;
    };

    zero_n z;
    ASSERT_NO_THROW(set_zero_n(z,0));

    // fail postcondition
    ASSERT_THROW(set_zero_n(z,4),rdmini::validation_failure);

    z.n=-3;
    // fail precondition
    ASSERT_THROW(set_zero_n(z,0),rdmini::validation_failure);
}

struct zero_n_safer: zero_n {
    int set_n(int value) {
        auto _(check_valid_guard(this));
        n=value;
    }

    int set_n_assert(int value) {
        auto _(assert_valid_guard(this));
        n=value;
    }
};

TEST(check_valid,check_valid_guard_this) {
    zero_n_safer z;

    ASSERT_NO_THROW(z.set_n(0));

    // fail postcondition
    ASSERT_THROW(z.set_n(4),rdmini::validation_failure);

    z.n=-3;
    // fail precondition
    ASSERT_THROW(z.set_n(0),rdmini::validation_failure);
}

TEST(check_valid,assert_valid_guard_this) {
    zero_n_safer z;

    EXPECT_TRUE(false);
    ASSERT_EXIT({z.set_n_assert(0); ::exit(0);},::testing::ExitedWithCode(0),"");

    // fail postcondition
    ASSERT_DEATH(z.set_n_assert(4),"validation failure.*postcondition.*n is positive");

    z.n=-3;
    // fail precondition
    ASSERT_DEATH(z.set_n_assert(0),"validation failure.*precondition.*n is FOO negative");
}
