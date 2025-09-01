#include <cmath>
#include <cstddef>
#include <doctest/doctest.h>
#include <vector>

#include <Chemistry/ChemistryModule.hpp>

TEST_CASE("Stats calculation")
{
    std::vector<double> real =
        {
            2, 2, 2, 2,              // species 1
            2.0, 0.01, 0.7, 0.5,     // species 2
            0.0, 0.0, 0.0, 0.0,      // species 3
            0.0, 0.0, 0.0, 0.0,      // species 4
            -2.5, -0.02, -0.7, -0.5, // species 5
            7.7, 6.01, 4.7, 0.5      // species 6
        };

    std::vector<double> pred =
        {
            2, 2, 2, 2,
            2.0, 0.02, 0.6, 0.5,
            0.1, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0,
            2.5, 0.01, 0.6, 0.5,
            2.8, 0.02, 0.7, 0.5
        };

    poet::ChemistryModule::error_stats stats(6, 5);
    poet::ChemistryModule::computeStats(real, pred, /*size_per_prop*/ 4, /*species_count*/ 6, stats);

    SUBCASE("Non-zero values")
    {

        // species 1 is ID, should stay 0
        CHECK_EQ(stats.mape[0], 0);
        CHECK_EQ(stats.rrsme[0], 0);

        /*
        mape species 2
        cell0: |(2.0 - 2.0)/2.0| = 0
        cell1: |(0.01 - 0.02)/0.01| = 1
        cell2: |(0.7 - 0.6)/0.7| = 0.142857143
        cell3: |(0.5 - 0.5)/0.5| = 0
        mape = 1.142857143/ 4 = 0.285714286 *100
        rrsme species 1
        squared err sum = 1.02040816
        rrsme = sqrt(1.02040816/4) = 0.50507627
        */

        CHECK_EQ(stats.mape[1], doctest::Approx(28.5714286).epsilon(1e-6));
        CHECK_EQ(stats.rrsme[1], doctest::Approx(0.50507627).epsilon(1e-6));
    }

    SUBCASE("Zero-denominator case")
    {
        /*
        mape species 3
        cell0: |(0.0 - 0.1)/0.0|
        cell1: |(0.0 - 0.0)/0.0|
        cell2: |(0.0 - 0.0)/0.0|
        cell3: |(0.0 - 0.0)/0.0|
        mape = 1 *100
        rrsme = 1
        */

        CHECK_EQ(stats.mape[2], 100.0);
        CHECK_EQ(stats.rrsme[2], 1.0);
    }

    SUBCASE("True and predicted values are zero")
    {
        /*
        mape species 4
        cell0: |(0.0 - 0.0)/0.0|
        cell1: |(0.0 - 0.0)/0.0|
        cell2: |(0.0 - 0.0)/0.0|
        cell3: |(0.0 - 0.0)/0.0|
        mape = 0.0
        rrsme = 0.0
        */

        CHECK_EQ(stats.mape[3], 0.0);
        CHECK_EQ(stats.rrsme[3], 0.0);
    }

    SUBCASE("Negative values")
    {
        /*
        mape species 5
        cell0: |(-2.5 - 2.5)/-2.5| = 2
        cell1: |(-0.02 - 0.01)/-0.02| = 1.5
        cell2: |(-0.7 - 0.6)/-0.7| = 1.85714286
        cell3: |(-0.5 - 0.5)/-0.5| = 2
        mape = (100.0 / 4) * 7.35714286 = 183.92857143
        rrsme = sqrt(13.6989796 / 4) = 1.85060663
        */

        CHECK_EQ(stats.mape[4], doctest::Approx(183.92857143).epsilon(1e-6));
        CHECK_EQ(stats.rrsme[4], doctest::Approx(1.85060663).epsilon(1e-6));
    }

    SUBCASE("Large differences")
    {
        /*
        mape species 6
        cell0: |(7.7 - 2.8)/7.7| = 0.63636364
        cell1: |(6.01 - 0.02)/6.01| = 0.99667221
        cell2: |(4.7 - 0.7)/4.7| = 0.85106383
        cell3: |(0.5 - 0.5)/0.5| = 0
        mape = (100.0 / 4) * 2.48409968 = 62.102492
        rrsme = sqrt(2,12262382 / 4) = 0.72846136
        */

        CHECK_EQ(stats.mape[5], doctest::Approx(62.102492).epsilon(1e-6));
        CHECK_EQ(stats.rrsme[5], doctest::Approx(0.72846136).epsilon(1e-6));
    }
}
