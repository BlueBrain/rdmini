/*
 * test_modelspec.cc: Unit testing of RDmini model specification
 * description: Test model specification defined in the equivalent 
 *              of a yaml file.  
 */

#include <string>

#include <gtest/gtest.h>

#include "rdmini/rdmodel.h"

std::string test_specification=
    "---\n"
    "model: modelTest\n"
    "cells:\n"
    "    wmvol:\n"
    "        volume: 1\n"
    "species:\n"
    "    name: A\n"
    "    concentration: 10\n"
    "    diffusivity: 1.0e-9\n"
    "species:\n"
    "    name: B\n"
    "    concentration: 10\n"
    "    diffusivity: 50.0\n"
    "reaction:\n"
    "    left: [ A, A, B ]\n"
    "    right: [ A, A, A ]\n"
    "    rate: 4e-5\n"
    "...\n";


TEST(yamlSpec,initialSpecification) {
    try {
        rdmini::rd_model M=rdmini::rd_model_read(test_specification,"modelTest");

        ASSERT_FALSE(M.species.empty());
        ASSERT_EQ(2,M.n_species());
        ASSERT_EQ(1,M.n_reactions());
        ASSERT_EQ(1,M.n_cells());

        /// Test model species specification in more detail
        ASSERT_EQ(0,M.species.index("A"));
        ASSERT_EQ(1,M.species.index("B"));
        ASSERT_TRUE(M.species.index("C")<0);
    }
    catch (rdmini::invalid_model &ex) {
        FAIL() << "unexpected invalid_model: " << ex.what();
    }
    catch (rdmini::model_io_error &ex) {
        FAIL() << "unexpected model_io_error: " << ex.what();
    }
}

TEST(yamlSpec, missingModelSpecification) {
    ASSERT_THROW(rdmini::rd_model_read(test_specification,"missingModel"),rdmini::model_io_error);
}

TEST(yamlSpec,negativeConcentrations) {
    std::string negative_concentration_spec1=
        "---\n"
        "model: modelTest1\n"
        "cells:\n"
        "    wmvol:\n"
        "        volume: 1\n"
        "species:\n"
        "    name: A\n"
        "    concentration: -10\n"
        "species:\n"
        "    name: B\n"
        "    concentration: 10\n"
        "reaction:\n"
        "    left: [ A, A, B ]\n"
        "    right: [ A, A, A ]\n"
        "    rate: 4e-5\n"
        "...\n";
 
    // Adding new species with negative values: Should throw 
    ASSERT_THROW(rdmini::rd_model_read(negative_concentration_spec1,"modelTest1"),rdmini::invalid_model);

    std::string negative_concentration_spec2=
        "---\n"
        "model: modelTest2\n"
        "cells:\n"
        "    wmvol:\n"
        "        volume: 1\n"
        "species:\n"
        "    name: A\n"
        "    concentration: 10\n"
        "species:\n"
        "    name: B\n"
        "    concentration: -10\n"
        "reaction:\n"
        "    left: [ A, A, B ]\n"
        "    right: [ A, A, A ]\n"
        "    rate: 4e-5\n"
        "...\n";

    ASSERT_THROW(rdmini::rd_model_read(negative_concentration_spec2,"modelTest2"),rdmini::invalid_model);
}

TEST(yamlSpec,negativeRate) {
    std::string negative_reaction_rate_spec=
        "---\n"
        "model: modelTest3\n"
        "cells:\n"
        "    wmvol:\n"
        "        volume: 1\n"
        "species:\n"
        "    name: A\n"
        "    concentration: -10\n"
        "species:\n"
        "    name: B\n"
        "    concentration: 10\n"
        "reaction:\n"
        "    left: [ A, A, B ]\n"
        "    right: [ A, A, A ]\n"
        "    rate: -4e-5\n"
        "...\n";

    ASSERT_THROW(rdmini::rd_model_read(negative_reaction_rate_spec,"modelTest3"),rdmini::invalid_model);
}

TEST(yamlSpec,zeroVolume) {
    std::string zero_volume_spec=
        "---\n"
        "model: modelTest4\n"
        "cells:\n"
        "    wmvol:\n"
        "        volume: 0\n"
        "species:\n"
        "    name: A\n"
        "    concentration: 10\n"
        "reaction:\n"
        "    right: [ A ]\n"
        "    rate: 4e-5\n"
        "...\n";

    ASSERT_THROW(rdmini::rd_model_read(zero_volume_spec,"modelTest4"),rdmini::invalid_model);
}

TEST(yamlSpec,negativeVolume) {
    std::string negative_volume_spec=
        "---\n"
        "model: modelTest5\n"
        "cells:\n"
        "    wmvol:\n"
        "        volume: -1\n"
        "species:\n"
        "    name: A\n"
        "    concentration: 10\n"
        "reaction:\n"
        "    right: [ A ]\n"
        "    rate: 4e-5\n"
        "...\n";

    ASSERT_THROW(rdmini::rd_model_read(negative_volume_spec,"modelTest5"),rdmini::invalid_model);
}

