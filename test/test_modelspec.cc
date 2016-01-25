/*
 * test_modelspec.cc: Unit testing of RDmini model specification
 * description: Test model specification and associated data structures 
 */

#include <gtest/gtest.h>
#include "rdmini/rdmodel.h"

class modelspec : public ::testing::Test
{
    protected:

    virtual void SetUp()
    {
        schnakenbergModl.name ="schnakenberg";
        specyA={"A",0.01,10.0};
        specyB={"B",0.02,20.0};
        schnakenbergModl.species.insert(specyA);
        schnakenbergModl.species.insert(specyB);

        std::multiset<int> left = {1,1,2};
        std::multiset<int> right = {1,1,1};
        reactionA={"reactionA",left,right,4e-5};
        schnakenbergModl.reactions.insert(reactionA);
    }

    virtual void TearDown()
    {
        schnakenbergModl.species.clear();
        schnakenbergModl.reactions.clear();
    }

    rd_model schnakenbergModl;
    species_info specyA, specyB;
    reaction_info reactionA;
};


TEST_F(modelspec,initialTest) {
    ASSERT_TRUE(true);
}

TEST_F(modelspec,initialSpecification) {
   /// This model is partially defined in the test header 
    ASSERT_TRUE(!schnakenbergModl.species.empty());
    ASSERT_TRUE(schnakenbergModl.n_species()==2);
    ASSERT_TRUE(schnakenbergModl.n_reactions()==1);
    ASSERT_TRUE(schnakenbergModl.n_cells()==0);
}


TEST_F(modelspec,addSpecies) {
    /// Completing model specification
    // Adding specy which does not exist
    species_info speciesC={"C",0.05,15.0}; 
    // Verifying specy specification 
    ASSERT_STREQ("C",speciesC.name.c_str());
    ASSERT_DOUBLE_EQ(0.05,speciesC.diffusivity);
    ASSERT_DOUBLE_EQ(15.0,speciesC.concentration);

    // Verifying that model is correctly updated  
    schnakenbergModl.species.insert(speciesC);
    ASSERT_TRUE(schnakenbergModl.n_species()==3);
    size_t index = schnakenbergModl.species.index("C"); 
    ASSERT_TRUE(index==2);
    ASSERT_STREQ("C",schnakenbergModl.species[index].name.c_str());
    ASSERT_DOUBLE_EQ(0.05,schnakenbergModl.species[index].diffusivity);
    ASSERT_DOUBLE_EQ(15.0,schnakenbergModl.species[index].concentration); 

    species_info speciesD={"D",-0.05,15.0};
    ASSERT_FALSE(speciesD.is_valid());
    species_info speciesE={"E",0.05,-15.0};
    ASSERT_FALSE(speciesE.is_valid());

    ASSERT_THROW(speciesD.check_valid(),rdmini::validation_failure);
}

TEST_F(modelspec,addReaction) {
    /// Adding regular reaction
    std::multiset<int> left = {1,1,2};
    std::multiset<int> right = {1,1,5,-5};
    reaction_info reactionB={"reactionB",left,right,10};   
    schnakenbergModl.reactions.insert(reactionB);

    ASSERT_TRUE(schnakenbergModl.n_reactions()==2);
    size_t index = schnakenbergModl.reactions.index("reactionA");
    ASSERT_TRUE(index==0);

    const auto &model_reactionA=schnakenbergModl.reactions[index];
    ASSERT_STREQ("reactionA",model_reactionA.name.c_str());
    ASSERT_DOUBLE_EQ(4e-5,model_reactionA.rate);
    ASSERT_TRUE(model_reactionA.is_valid());

    index = schnakenbergModl.reactions.index("reactionB");
    ASSERT_TRUE(index==1);

    const auto &model_reactionB=schnakenbergModl.reactions[index];
    ASSERT_STREQ("reactionB",model_reactionB.name.c_str());
    ASSERT_DOUBLE_EQ(10,model_reactionB.rate);
    ASSERT_TRUE(model_reactionB.left.size()==3);
    ASSERT_TRUE(model_reactionB.right.size()==4);  
    ASSERT_TRUE(model_reactionB.is_valid());
    ASSERT_NO_THROW(model_reactionB.check_valid());

    /// Add improper reaction
    reaction_info reactionC={"reactionC",left,right,-10};   
    schnakenbergModl.reactions.insert(reactionC);

    ASSERT_THROW(schnakenbergModl.reactions["reactionC"].check_valid(),rdmini::validation_failure);
}

/*
TEST(modelspec,reactions) {


}

TEST(modelspec,cell) {

}


TEST(modelspec,model) {


}
*/



