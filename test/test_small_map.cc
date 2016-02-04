#include <utility>
#include <cmath>
#include <gtest/gtest.h>

#include "tiny_map.h"

template <typename T>
class small_map: public ::testing::Test {
public:
    typedef typename T::value_type value_type;
};

// use this to check for correct ctor, dtor behaviour in tiny_map

int g_dtor_count=0;
int g_ctor_count=0;

void reset_counts() {
    g_dtor_count=g_ctor_count=0;
}

struct int_nontrivial {
    int_nontrivial() { ++g_ctor_count; }
    int_nontrivial(int n_): n(n_) { ++g_ctor_count; }
    int_nontrivial(const int_nontrivial &x): n(x.n) { ++g_ctor_count; }
    int_nontrivial(int_nontrivial &&x): n(x.n) { ++g_ctor_count; }

    int_nontrivial &operator=(const int_nontrivial &x) { n=x.n; return *this; }
    int_nontrivial &operator=(int_nontrivial &&x) { n=x.n; return *this; }

    ~int_nontrivial() { ++g_dtor_count; }

    operator int() const { return n; }
    int n;
};


using map_types=::testing::Types<small_map<int,int>,small_map<int_nontrivial,int_nontrivial>>;
TYPED_TEST_CASE(small_map,map_types);


TYPED_TEST(small_map,ctor) {
    using value_type=typename TestFixture::value_type;
    using map=TypeParam;

    reset_counts();

    {
        map m_ilist({{1,1},{3,2},{4,1},{3,7}});

        ASSERT_EQ(3,m_ilist.size());
    }

    ASSERT_EQ(g_dtor_count,g_ctor_count);

    {
        std::vector<value_type> ns={{1,1},{3,2},{4,1},{3,7}};
        map m_ipair(ns.begin(),ns.end());

        ASSERT_EQ(3,m_ipair.size());

        map m_copy(m_ipair);
        
        ASSERT_EQ(m_ipair.size(),m_copy.size());

        map m_move(std::move(m_ipair));
        
        ASSERT_EQ(m_copy.size(),m_move.size());
    }

    ASSERT_EQ(g_dtor_count,g_ctor_count);
}

TYPED_TEST(small_map,empty) {
    using map=TypeParam;

    map m;
    ASSERT_TRUE(m.empty());
    ASSERT_EQ(m.end(),m.begin());
    ASSERT_EQ(0,m.size());
}

TYPED_TEST(small_map,clear) {
    using map=TypeParam;

    reset_counts();

    {
        map m_ilist({{1,1},{3,2},{4,1},{3,7}});
        m_ilist.clear();
        ASSERT_TRUE(m_ilist.empty());
    }

    ASSERT_EQ(g_dtor_count,g_ctor_count);
}

TYPED_TEST(small_map,equality) {
    using map=TypeParam;
    
    map m1({{1,1},{3,2},{4,1},{3,7}});
    map m2({{3,7},{1,1},{4,1}});
    map m3({{3,7},{1,1},{4,1},{1,2}});
    map m4({{3,7},{1,1},{4,1},{5,6}});

    ASSERT_EQ(m1,m2);
    ASSERT_NE(m1,m3);
    ASSERT_NE(m1,m4);
}

TYPED_TEST(small_map,insert) {
    using value_type=typename TestFixture::value_type;
    using map=TypeParam;

    reset_counts();

    {
        map m;
        value_type v(3,8);

        m.insert(v);
        m.insert(value_type(4,9));

        map m_bis({{4,9},{3,8}});
        ASSERT_EQ(m,m_bis);
    }

    ASSERT_EQ(g_dtor_count,g_ctor_count);
}

TYPED_TEST(small_map,swap) {
    using map=TypeParam;

    reset_counts();

    {
        map m1({{3,9},{4,1}});
        map m2({{3,7},{1,1},{4,1},{5,6}});

        map m1_copy=m1;
        map m2_copy=m2;

        m1.swap(m2);

        ASSERT_EQ(4,m1.size());
        ASSERT_EQ(m2_copy,m1);

        ASSERT_EQ(2,m2.size());
        ASSERT_EQ(m1_copy,m2);
    }

    ASSERT_EQ(g_dtor_count,g_ctor_count);
}

TYPED_TEST(small_map,count) {
    using map=TypeParam;

    map m1({{1,2},{3,2},{3,5},{4,5}});
    ASSERT_EQ(1,m1.count(1));
    ASSERT_EQ(0,m1.count(2));
    ASSERT_EQ(1,m1.count(3));
    ASSERT_EQ(1,m1.count(4));
    ASSERT_EQ(0,m1.count(5));
}

TYPED_TEST(small_map,erase) {
    using map=TypeParam;

    reset_counts();

    {
        map m1({{1,2},{3,4},{4,5},{3,7}});
        ASSERT_EQ(1,m1.erase(4));
        ASSERT_EQ(0,m1.erase(4));
        ASSERT_EQ(1,m1.erase(3));
        ASSERT_EQ(1,m1.size());
    }

    ASSERT_EQ(g_dtor_count,g_ctor_count);
}

TYPED_TEST(small_map,iter_erase) {
    using value_type=typename TestFixture::value_type;
    using map=TypeParam;

    reset_counts();

    {
        map m1({{1,2},{3,2},{3,4},{4,5},{5,6},{7,8}});
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

TYPED_TEST(small_map,bracket) {
    using map=TypeParam;

    reset_counts();

    {
        map m1;
        m1[3]=5;

        ASSERT_EQ(1,m1.size());
        ASSERT_EQ(5,m1[3]);

        m1[3]=4;
        ASSERT_EQ(1,m1.size());
        ASSERT_EQ(4,m1[3]);
    }

    ASSERT_EQ(g_dtor_count,g_ctor_count);
}

TYPED_TEST(small_map,at) {
    using map=TypeParam;

    reset_counts();

    {
        map m1({{1,2},{2,3}});
        ASSERT_NO_THROW(m1.at(1));
        ASSERT_NO_THROW(m1.at(2));
        ASSERT_EQ(2,m1.at(1));
        ASSERT_EQ(3,m1.at(2));

        m1.at(2)=5;
        ASSERT_EQ(5,m1.at(2));

        ASSERT_THROW(m1.at(9),std::out_of_range);

        // check const version too...
        const map &m2(m1);
        ASSERT_NO_THROW(m2.at(1));
        ASSERT_NO_THROW(m2.at(2));
        ASSERT_EQ(2,m2.at(1));
        ASSERT_EQ(5,m2.at(2));

        ASSERT_THROW(m2.at(9),std::out_of_range);
    }

    ASSERT_EQ(g_dtor_count,g_ctor_count);
}

