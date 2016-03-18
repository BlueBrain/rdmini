#include <algorithm>
#include <vector>
#include <random>

#include <gtest/gtest.h>

#include "rdmini/ssa_direct.h"

class ssa_selector : public ::testing::Test
{
    protected:

    virtual void SetUp()
    {
        /// Define vector of propensities
        prop_size = 100;
        propensities.resize(prop_size);

        /// Randomize values of propensities
        std::uniform_real_distribution<double> UniformDistrib(0.5,1.0);
        for (size_t i=0; i<prop_size; ++i) {
            propensities[i]=std::ldexp(UniformDistrib(R),-i);
        }
        std::shuffle(propensities.begin(),propensities.end(),R);

        /// Create ssa specification with defined propensities
        selector.reset(prop_size);
    }

    virtual void TearDown() { }

    std::minstd_rand R;
    std::vector<double> propensities;
    size_t prop_size;

    rdmini::ssa_direct<size_t,double> selector;
};


TEST_F(ssa_selector,initialSpecification) {
    /// This model is partially defined in the test header 
    ASSERT_EQ(prop_size,selector.size()); 

    /// Update propensities to ssa data structure
    /// and test value of propensities from data structure
    double total=0.0;
    for (size_t i=0; i<prop_size; ++i)
    {
        selector.update(i,propensities[i]);
        total += propensities[i];    
    }
    
    /// Verifying value of propensities stored
    for (size_t i=0; i<prop_size; ++i)
        ASSERT_EQ(propensities[i],selector.propensity(i));

    /// Verifying update of value
    ASSERT_DOUBLE_EQ(total, selector.total_propensity());
}


