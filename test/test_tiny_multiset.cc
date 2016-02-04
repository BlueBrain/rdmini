#include <utility>
#include <cmath>
#include <gtest/gtest.h>

#include "rdmini/tiny_multiset.h"

template <typename T>
class xmultiset: public ::testing::Test {
public:
    typedef typename T::value_type value_type;
};

// use this to check for correct ctor, dtor behaviour in tiny_multisets

int g_dtor_count=0;
int g_ctor_count=0;

void reset_counts() {
    g_dtor_count=g_ctor_count=0;
}

struct int_nontrivial {
    int_nontrivial(int n_): n(n_) { ++g_ctor_count; }
    int_nontrivial(const int_nontrivial &x): n(x.n) { ++g_ctor_count; }
    int_nontrivial(int_nontrivial &&x): n(x.n) { ++g_ctor_count; }

    int_nontrivial &operator=(const int_nontrivial &x) { n=x.n; return *this; }
    int_nontrivial &operator=(int_nontrivial &&x) { n=x.n; return *this; }

    ~int_nontrivial() { ++g_dtor_count; }

    operator int() const { return n; }
    int n;
};

using multiset_types=::testing::Types<tiny_multiset<int,20>,tiny_multiset<int_nontrivial,20>,small_multiset<int>>;
TYPED_TEST_CASE(xmultiset,multiset_types);


TYPED_TEST(xmultiset,small_ctor) {
    using value_type=typename TestFixture::value_type;
    using mset=TypeParam;

    reset_counts();

    {
        mset m_ilist({1,2,3,2,3,4,3,4,5});

        ASSERT_EQ(9,m_ilist.size());
    }

    ASSERT_EQ(g_dtor_count,g_ctor_count);

    {
        std::vector<value_type> ns={3,3,4,4,4,5,5,5,5,3};
        mset m_ipair(ns.begin(),ns.end());

        ASSERT_EQ(10,m_ipair.size());

        mset m_copy(m_ipair);
        
        ASSERT_EQ(m_ipair.size(),m_copy.size());

        mset m_move(std::move(m_ipair));
        
        ASSERT_EQ(m_copy.size(),m_move.size());
    }

    ASSERT_EQ(g_dtor_count,g_ctor_count);
}

TYPED_TEST(xmultiset,empty) {
    using value_type=typename TestFixture::value_type;
    using mset=TypeParam;

    mset m;
    ASSERT_TRUE(m.empty());
    ASSERT_EQ(m.end(),m.begin());
    ASSERT_EQ(0,m.size());
}

TYPED_TEST(xmultiset,clear) {
    using value_type=typename TestFixture::value_type;
    using mset=TypeParam;

    reset_counts();

    {
        mset m_ilist({1,2,3,2,3,4,3,4,5});
        m_ilist.clear();
        ASSERT_TRUE(m_ilist.empty());
    }

    ASSERT_EQ(g_dtor_count,g_ctor_count);
}

TYPED_TEST(xmultiset,equality) {
    using value_type=typename TestFixture::value_type;
    using mset=TypeParam;
    
    mset m1({1,2,3,2,3,4,3,4,5});
    mset m2({5,4,4,2,3,2,3,3,1});
    mset m3({5,4,4,2,3,2,3,3});
    mset m4({5,4,4,2,3,2,3,3,2});

    ASSERT_EQ(m1,m2);
    ASSERT_NE(m1,m3);
    ASSERT_NE(m1,m4);
}

TYPED_TEST(xmultiset,insert) {
    using value_type=typename TestFixture::value_type;
    using mset=TypeParam;

    reset_counts();

    {
        mset m;
        value_type v(3);

        m.insert(v);
        m.insert(value_type(4));

        mset m_bis({4,3});
        ASSERT_EQ(m,m_bis);
    }

    ASSERT_EQ(g_dtor_count,g_ctor_count);
}

TYPED_TEST(xmultiset,swap) {
    using value_type=typename TestFixture::value_type;
    using mset=TypeParam;

    reset_counts();

    {
        mset m1({1,2,3,2,3,4,3,4,5});
        mset m1_copy=m1;

        mset m2({7,6,6,5});
        mset m2_copy=m2;

        m1.swap(m2);

        ASSERT_EQ(4,m1.size());
        ASSERT_EQ(m2_copy,m1);

        ASSERT_EQ(9,m2.size());
        ASSERT_EQ(m1_copy,m2);
    }

    ASSERT_EQ(g_dtor_count,g_ctor_count);
}

TYPED_TEST(xmultiset,count) {
    using value_type=typename TestFixture::value_type;
    using mset=TypeParam;

    mset m1({1,2,3,2,3,4,3,4,5});
    ASSERT_EQ(1,m1.count(1));
    ASSERT_EQ(2,m1.count(2));
    ASSERT_EQ(3,m1.count(3));
    ASSERT_EQ(2,m1.count(4));
    ASSERT_EQ(1,m1.count(5));
}

TYPED_TEST(xmultiset,erase) {
    using value_type=typename TestFixture::value_type;
    using mset=TypeParam;

    reset_counts();

    {
        mset m1({1,2,3,4,4,5});
        ASSERT_EQ(2,m1.erase(4));
        ASSERT_EQ(0,m1.erase(4));
        ASSERT_EQ(1,m1.erase(3));
        ASSERT_EQ(3,m1.size());
    }

    ASSERT_EQ(g_dtor_count,g_ctor_count);
}

TYPED_TEST(xmultiset,iter_erase) {
    using value_type=typename TestFixture::value_type;
    using mset=TypeParam;

    reset_counts();

    {
        mset m1({1,2,3,2,3,4,3,4,5});
        size_t initial_size=m1.size();
        int leap=3;

        auto i=m1.begin();
        std::advance(i,leap);

        int n_erase=0;
        while (i!=m1.end()) {
            i=m1.erase(i);
            ++n_erase;
            ASSERT_EQ(initial_size,n_erase+m1.size());
        }

        ASSERT_EQ(leap,m1.size());
    }

    ASSERT_EQ(g_dtor_count,g_ctor_count);
}

template <typename T>
class xmultiset_nonstd_eq: public ::testing::Test {
public:
    typedef typename T::value_type value_type;
    typedef typename T::key_equal key_equal;
};

// non-standard stateful equality functor
struct eq_mod_k {
    eq_mod_k(): k(2) {}
    explicit eq_mod_k(int k_): k(k_) {}
    bool operator()(int a,int b) const { return std::abs(a-b)%k==0; }

    int k;
};

using multiset_nonstd_eq_types=::testing::Types<tiny_multiset<int,20,eq_mod_k>,tiny_multiset<int_nontrivial,20,eq_mod_k>,small_multiset<int,eq_mod_k>>;
TYPED_TEST_CASE(xmultiset_nonstd_eq,multiset_nonstd_eq_types);

TYPED_TEST(xmultiset_nonstd_eq,count) {
    using value_type=typename TestFixture::value_type;
    using mset=TypeParam;

    // default eq_mod_k is mod 2 (i.e. equal if same parity)
    mset m1({1,2,3,4,5});
    ASSERT_EQ(3,m1.count(1));
    ASSERT_EQ(2,m1.count(2));
    
    // test with stateful eq_mod_k, k==3
    mset m2({1,2,3,4,5},eq_mod_k(3));
    ASSERT_EQ(2,m2.count(1));
    ASSERT_EQ(2,m2.count(2));
    ASSERT_EQ(1,m2.count(3));
}

TYPED_TEST(xmultiset_nonstd_eq,erase) {
    using value_type=typename TestFixture::value_type;
    using mset=TypeParam;

    mset m2({1,2,3,4,5},eq_mod_k(3));
    ASSERT_EQ(5,m2.size());

    size_t k=m2.erase(1);
    ASSERT_EQ(2,k);
    ASSERT_EQ(3,m2.size());

    k=m2.erase(2);
    ASSERT_EQ(2,k);
    ASSERT_EQ(1,m2.size());
}

TYPED_TEST(xmultiset_nonstd_eq,key_eq) {
    using value_type=typename TestFixture::value_type;
    using key_equal=typename TestFixture::key_equal;
    using mset=TypeParam;

    mset m({1,2,3,4,5},eq_mod_k(3));
    key_equal eq=m.key_eq();
    ASSERT_EQ(3,eq.k);
}

