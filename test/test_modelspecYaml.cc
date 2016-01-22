/*
 * test_modelspec.cc: Unit testing of RDmini model specification
 * description: Test model specification defined in the equivalent 
 *              of a yaml file.  
 */

#include <gtest/gtest.h>
#include "rdmini/rdmodel.h"
#include <iostream>

class yamlSpec : public ::testing::Test
{
    protected:

    virtual void SetUp()
    {
        YamlSpecification= "---\nmodel: modelTest\n\n\
cells:\n    wmvol:\n        volume: 1\n\n\
species:\n    name: A\n    concentration: 10\n    diffusivity: 1.0e-9\n\n\
species:\n    name: B\n    concentration: 10\n    diffusivity: 50.0 \n\n\
reaction:\n    left: [ A, A, B ]\n    right: [ A, A, A ]\n    rate: 4e-5\n\n...";
        M = rd_model_read(YamlSpecification,"modelTest");
    }

    virtual void TearDown()
    {
        M.species.clear();
        M.reactions.clear();
    }

    std::string YamlSpecification; 
    rd_model M; 
};


TEST_F(yamlSpec,initialTest) {
    ASSERT_TRUE(true);
}

TEST_F(yamlSpec,initialSpecification) {
    /// This model is partially defined in the Setup function of the fixture 
    ASSERT_TRUE(!M.species.empty());
    ASSERT_TRUE(M.n_species()==2);
    ASSERT_TRUE(M.n_reactions()==1);
    ASSERT_TRUE(M.n_cells()==1);

    /// Test model specy specification in more details
    ASSERT_TRUE(M.species.index("A")==0);
    ASSERT_TRUE(M.species.index("B")==1);
    ASSERT_TRUE(M.species.index("C")==-1);
}

TEST_F(yamlSpec, missingModlSpecification) {
    ASSERT_THROW(rd_model_read(YamlSpecification,"wrongTest"), model_io_error);
}

TEST_F(yamlSpec, readSpecTwice)
{
    /// Read model specification again
    M = rd_model_read(YamlSpecification,"modelTest");
    ASSERT_TRUE(M.n_species()==2);
    ASSERT_TRUE(M.n_reactions()==1);
    ASSERT_TRUE(M.n_cells()==1);
} 

TEST_F(yamlSpec,negativeConcentration1) {
    std::string wrongSpecConcentration= "---\nmodel: modelTest1\n\n\
cells:\n    wmvol:\n        volume: 1\n\n\
species:\n    name: A\n    concentration: -10\n\n\
species:\n    name: B\n    concentration: 10\n\n\
reaction:\n    left: [ A, A, B ]\n    right: [ A, A, A ]\n    rate: 4e-5\n\n...";
 
    // Adding new specy with negative values: Should throw 
    ASSERT_THROW(rd_model_read(wrongSpecConcentration,"modelTest1"), model_incorrectBiologicalValue_error);
}

TEST_F(yamlSpec,negativeConcentration2) {
    std::string wrongSpecConcentration= "---\nmodel: modelTest2\n\n\
cells:\n    wmvol:\n        volume: 1\n\n\
species:\n    name: A\n    concentration: 10\n\n\
species:\n    name: B\n    concentration: -10\n\n\
reaction:\n    left: [ A, A, B ]\n    right: [ A, A, A ]\n    rate: 4e-5\n\n...";
 
    // Adding new specy with negative values: Should throw 
    ASSERT_THROW(rd_model_read(wrongSpecConcentration,"modelTest2"), model_incorrectBiologicalValue_error);
}

TEST_F(yamlSpec,negativeRate) {
    std::string WrongSpecRate= "---\nmodel: modelTest3\n\n\
cells:\n    wmvol:\n        volume: 1\n\n\
species:\n    name: A\n    concentration: 10\n\n\
species:\n    name: B\n    concentration: 10\n\n\
reaction:\n    left: [ A, A, B ]\n    right: [ A, A, A ]\n    rate: -4e-5\n\n...";

    // Adding new specy with negative values: Should throw 
    ASSERT_THROW(rd_model_read(WrongSpecRate,"modelTest3"), model_incorrectBiologicalValue_error);
}

TEST_F(yamlSpec,negativeVolume) {
    std::string wrongSpecVol= "---\nmodel: modelTest4\n\n\
cells:\n    wmvol:\n        volume: -1\n\n\
species:\n    name: A\n    concentration: 10\n\n\
species:\n    name: B\n    concentration: 10\n\n\
reaction:\n    left: [ A, A, B ]\n    right: [ A, A, A ]\n    rate: 4e-5\n\n...";

    // Adding new specy with negative values: Should throw 
    ASSERT_THROW(rd_model_read(wrongSpecVol,"modelTest4"), model_incorrectBiologicalValue_error);
}

