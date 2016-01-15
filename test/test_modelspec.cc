/*
 * test_modelspec.cc: Unit testing of RDmini model specification
 * description: Test model specification and associated data structures 
 */

#include <gtest/gtest.h>
#include "rdmini/rdmodel.h"

class rdmini : public ::testing::Test
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


TEST_F(rdmini,initialTest) {
    ASSERT_TRUE(true);
}

TEST_F(rdmini,initialSpecification) {
   /// This model is partially defined in the test header 
    ASSERT_TRUE(!schnakenbergModl.species.empty());
    ASSERT_TRUE(schnakenbergModl.n_species()==2);
    ASSERT_TRUE(schnakenbergModl.n_reactions()==1);
    ASSERT_TRUE(schnakenbergModl.n_cells()==0);
}


TEST_F(rdmini,addSpecy) {
    /// Completing model specification
    // Adding specy which does not exist
    species_info specyC={"C",0.05,15.0}; 
    // Verifying specy specification 
    ASSERT_STREQ(specyC.name.c_str(),"C");
    ASSERT_DOUBLE_EQ(0.05,specyC.diffusivity);
    ASSERT_DOUBLE_EQ(15.0,specyC.concentration);

    // Verifying that model is correctly updated  
    schnakenbergModl.species.insert(specyC);
    ASSERT_TRUE(schnakenbergModl.n_species()==3);
    size_t index = schnakenbergModl.species.index("C"); 
    ASSERT_TRUE(index==2);
    ASSERT_STREQ(schnakenbergModl.species[index].name.c_str(),"C");
    ASSERT_DOUBLE_EQ(0.05,schnakenbergModl.species[index].diffusivity);
    ASSERT_DOUBLE_EQ(15.0,schnakenbergModl.species[index].concentration); 
}

TEST_F(rdmini,addReaction) {
    /// Adding regular reaction
    std::multiset<int> left = {1,1,2};
    std::multiset<int> right = {1,1,5,-5};
    reaction_info reactionB={"reactionB",left,right,10};   
    schnakenbergModl.reactions.insert(reactionB);

    ASSERT_TRUE(schnakenbergModl.n_reactions()==2);
    size_t index = schnakenbergModl.reactions.index("reactionA");
    ASSERT_TRUE(index==0);
    ASSERT_STREQ(schnakenbergModl.reactions[index].name.c_str(),"reactionA");
    ASSERT_DOUBLE_EQ(4e-5, schnakenbergModl.reactions[index].rate); 

    index = schnakenbergModl.reactions.index("reactionB");
    ASSERT_TRUE(index==1);
    ASSERT_STREQ(schnakenbergModl.reactions[index].name.c_str(),"reactionB");
    ASSERT_DOUBLE_EQ(10, schnakenbergModl.reactions[index].rate);
    ASSERT_TRUE( schnakenbergModl.reactions[index].left.size()==3);
    ASSERT_TRUE( schnakenbergModl.reactions[index].right.size()==4);  
    
    /// Check for reaction values
/*    std::multiset<int>::value_compare comparison = right.value_comp();
    std::multiset<int>::iterator it = right.begin();
    do {
        std::cout << ' ' << *it;
    } while (mycomp(*it++, */
}

/*
TEST(rdmini,reactions) {


}

TEST(rdmini,cell) {

}


TEST(rdmini,model) {


}
*/



