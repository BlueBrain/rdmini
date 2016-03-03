/*
 * test_modelspec.cc: Unit testing of RDmini model specification
 * description: Test model specification and associated data structures 
 */

#include <algorithm>
#include <vector>
#include <random>

#include <gtest/gtest.h>
#include "rdmini/rdmodel.h"
#include "rdmini/ssa_direct.h"

class ssa : public ::testing::Test
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
        ssa_solver.reset(prop_size);
    }

    virtual void TearDown() { }

    std::minstd_rand R;
    typedef ssa_direct<size_t,double> S;// ssa_solver;
    std::vector<double> propensities;
    size_t prop_size;
    S ssa_solver;
};


TEST_F(ssa,initialSpecification) {
    /// This model is partially defined in the test header 
    ASSERT_TRUE(ssa_solver.size()==prop_size); 

    /// Update propensities to ssa data structure
    /// and test value of propensities from data structure
    double total=0.0;
    for (size_t i=0; i<prop_size; ++i)
    {
        ssa_solver.update(i,propensities[i]);
        total += propensities[i];    
    }
    
    /// Verifying value of propensities stored
    for (size_t i=0; i<prop_size; ++i)
        ASSERT_TRUE(ssa_solver.propensity(i)==propensities[i]);

    /// Verifying update of value
    ASSERT_DOUBLE_EQ(ssa_solver.total_propensity(), total);
}


