/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.apache.commons.math4.stat.inference;

import org.apache.commons.math4.exception.DimensionMismatchException;
import org.apache.commons.math4.exception.NotPositiveException;
import org.apache.commons.math4.exception.NotStrictlyPositiveException;
import org.apache.commons.math4.exception.OutOfRangeException;
import org.apache.commons.math4.exception.ZeroException;
import org.apache.commons.math4.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

/**
 * Test cases for the GTest class.
 *
 * Data for the tests are from p64-69 in: McDonald, J.H. 2009. Handbook of
 * Biological Statistics (2nd ed.). Sparky House Publishing, Baltimore,
 * Maryland.
 *
 */
public class GTestTest {

    protected GTest testStatistic = new GTest();

    @Test
    public void testGTestGoodnesOfFit1() {
        final double[] exp = new double[]{
            3d, 1d
        };

        final long[] obs = new long[]{
            423, 133
        };

        Assert.assertEquals("G test statistic",
                0.348721, testStatistic.g(exp, obs), 1E-6);
        final double p_gtgf = testStatistic.gTest(exp, obs);
        Assert.assertEquals("g-Test p-value", 0.55483, p_gtgf, 1E-5);

        Assert.assertFalse(testStatistic.gTest(exp, obs, 0.05));
    }

    @Test
    public void testGTestGoodnesOfFit2() {
        final double[] exp = new double[]{
            0.54d, 0.40d, 0.05d, 0.01d
        };

        final long[] obs = new long[]{
            70, 79, 3, 4
        };
        Assert.assertEquals("G test statistic",
                13.144799, testStatistic.g(exp, obs), 1E-6);
        final double p_gtgf = testStatistic.gTest(exp, obs);
        Assert.assertEquals("g-Test p-value", 0.004333, p_gtgf, 1E-5);

        Assert.assertTrue(testStatistic.gTest(exp, obs, 0.05));
    }

    @Test
    public void testGTestGoodnesOfFit3() {
        final double[] exp = new double[]{
            0.167d, 0.483d, 0.350d
        };

        final long[] obs = new long[]{
            14, 21, 25
        };

        Assert.assertEquals("G test statistic",
                4.5554, testStatistic.g(exp, obs), 1E-4);
        // Intrinisic (Hardy-Weinberg proportions) P-Value should be 0.033
        final double p_gtgf = testStatistic.gTestIntrinsic(exp, obs);
        Assert.assertEquals("g-Test p-value", 0.0328, p_gtgf, 1E-4);

        Assert.assertFalse(testStatistic.gTest(exp, obs, 0.05));
    }

    @Test
    public void testGTestIndependence4() {
        // Target values proofed using R version 3.6.3 Library DescTools and function GTest
        // Testdata of this test taken from similar chi square test

        long[][] counts = { {40, 22, 43}, {91, 21, 28}, {60, 10, 22}};
        Assert.assertEquals( "G statistic", 22.860523581, testStatistic.g(counts), 1E-9);
        Assert.assertEquals("g-test p-value", 0.000135014574778, testStatistic.gTest(counts), 1E-9);
        Assert.assertTrue("g-test test reject", testStatistic.gTest(counts, 0.0002));
        Assert.assertFalse("g-test test accept", testStatistic.gTest(counts, 0.0001));

        long[][] counts2 = {{10, 15}, {30, 40}, {60, 90} };
        Assert.assertEquals( "G statistic", 0.168482556026, testStatistic.g(counts2), 1E-9);
        Assert.assertEquals("g-test p-value",0.919209444309, testStatistic.gTest(counts2), 1E-9);
        Assert.assertFalse("g-test test accept", testStatistic.gTest(counts2, 0.1));

        // ragged input array
        long[][] counts3 = { {40, 22, 43}, {91, 21, 28}, {60, 10}};
        try {
            testStatistic.g(counts3);
            Assert.fail("Expecting DimensionMismatchException");
        } catch (DimensionMismatchException ex) {
            // expected
        }

        // insufficient data
        long[][] counts4 = {{40, 22, 43}};
        try {
            testStatistic.g(counts4);
            Assert.fail("Expecting DimensionMismatchException");
        } catch (DimensionMismatchException ex) {
            // expected
        }
        long[][] counts5 = {{40}, {40}, {30}, {10}};
        try {
            testStatistic.g(counts5);
            Assert.fail("Expecting DimensionMismatchException");
        } catch (DimensionMismatchException ex) {
            // expected
        }

        // negative counts
        long[][] counts6 = {{10, -2}, {30, 40}, {60, 90} };
        try {
            testStatistic.g(counts6);
            Assert.fail("Expecting NotPositiveException");
        } catch (NotPositiveException ex) {
            // expected
        }

        // bad alpha
        try {
            testStatistic.gTest(counts, 0);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
    }

    /** Contingency table containing zeros - PR # 32531 */
    @Test
    public void testGZeroCount() {
        // Testdata of this test taken from similar chi square test
        // Target values proofed using R version 3.6.3 Library DescTools and function GTest

        long[][] counts = { {40, 0, 4}, {91, 1, 2}, {60, 2, 0}};
        Assert.assertEquals( "G statistic", 9.99660546076,
                testStatistic.g(counts), 1E-9);
        Assert.assertEquals("g-test p-value", 0.0404849013949,
                testStatistic.gTest(counts), 1E-9);
    }

    @Test
    public void testGTestIndependence1() {
        long[][] counts = new long[][]{{268, 199, 42}, {807, 759, 184}};

        // Test for the 2 row DataSetsComparison Methods
        final double g = testStatistic.gDataSetsComparison(counts[0], counts[1]);
        Assert.assertEquals("G test statistic", 7.3008170, g, 1E-6);
        final double p_gti = testStatistic.gTestDataSetsComparison(counts[0], counts[1]);
        Assert.assertEquals("g-Test p-value", 0.0259805, p_gti, 1E-6);
        Assert.assertTrue(testStatistic.gTestDataSetsComparison(counts[0], counts[1], 0.05));

        // Same results for more flexible (rows>=2) array-of-array version:
        final double g2 = testStatistic.g(counts);
        Assert.assertEquals("G test statistic", 7.3008170, g2, 1E-6);
        final double p_gti2 = testStatistic.gTest(counts);
        Assert.assertEquals("g-Test p-value", 0.0259805, p_gti2, 1E-6);
        Assert.assertTrue(testStatistic.gTest(counts, 0.05));
    }

    @Test
    public void testGTestIndependence2() {
        long[][] counts = new long[][]{{127, 99, 264}, {116, 67, 161}};

        // Test for the 2 row DataSetsComparison Methods
        final double g = testStatistic.gDataSetsComparison(counts[0], counts[1]);
        Assert.assertEquals("G test statistic", 6.227288, g, 1E-6);
        final double p_gti = testStatistic.gTestDataSetsComparison(counts[0], counts[1]);
        Assert.assertEquals("g-Test p-value", 0.04443, p_gti, 1E-5);
        Assert.assertTrue(testStatistic.gTestDataSetsComparison(counts[0], counts[1], 0.05));

        // Same results for more flexible (rows>=2) array-of-array version:
        final double g2 = testStatistic.g(counts);
        Assert.assertEquals("G test statistic", 6.227288, g2, 1E-6);
        final double p_gti2 = testStatistic.gTest(counts);
        Assert.assertEquals("g-Test p-value", 0.04443, p_gti2, 1E-5);
        Assert.assertTrue(testStatistic.gTest(counts, 0.05));
    }

    @Test
    public void testGTestIndependence3() {
        long[][] counts = new long[][]{{190, 149}, {42, 49}};

        // Test for the 2 row DataSetsComparison Methods
        final double g = testStatistic.gDataSetsComparison(counts[0], counts[1]);
        Assert.assertEquals("G test statistic",2.8187, g, 1E-4);
        final double p_gti = testStatistic.gTestDataSetsComparison(counts[0], counts[1]);
        Assert.assertEquals("g-Test p-value", 0.09317325, p_gti, 1E-6);
        Assert.assertFalse(testStatistic.gTestDataSetsComparison(counts[0], counts[1], 0.05));

        // Same results for more flexible (rows>=2) array-of-array version:
        final double g2 = testStatistic.g(counts);
        Assert.assertEquals("G test statistic",2.8187, g2, 1E-4);
        final double p_gti2 = testStatistic.gTest(counts);
        Assert.assertEquals("g-Test p-value", 0.09317325, p_gti2, 1E-6);
        Assert.assertFalse(testStatistic.gTest(counts, 0.05));
    }

    @Test
    public void testGTestSetsComparisonBadCounts() {
        long[] observed1 = {10, -1, 12, 10, 15};
        long[] observed2 = {15, 10, 10, 15, 5};
        try {
            testStatistic.gTest(
                    new long[][]{observed1, observed2});
            Assert.fail("Expecting NotPositiveException - negative count");
        } catch (NotPositiveException ex) {
            // expected
        }
        long[] observed3 = {10, 0, 12, 10, 15};
        long[] observed4 = {15, 0, 10, 15, 5};
        try {
            testStatistic.gTest(
                    new long[][]{observed3, observed4});
            Assert.fail("Expecting ZeroException - double 0's");
        } catch (ZeroException ex) {
            // expected
        }
        long[] observed5 = {10, 10, 12, 10, 15};
        long[] observed6 = {0, 0, 0, 0, 0};
        try {
            testStatistic.gTestDataSetsComparison(
                    observed5, observed6);
            Assert.fail("Expecting ZeroException - vanishing counts");
        } catch (ZeroException ex) {
            // expected
        }
    }

    @Test
    public void testUnmatchedArrays() {
        final long[] observed = { 0, 1, 2, 3 };
        final double[] expected = { 1, 1, 2 };
        final long[] observed2 = {3, 4};
        try {
            testStatistic.gTest(expected, observed);
            Assert.fail("arrays have different lengths, DimensionMismatchException expected");
        } catch (DimensionMismatchException ex) {
            // expected
        }
        try {
            testStatistic.gTestDataSetsComparison(observed, observed2);
            Assert.fail("arrays have different lengths, DimensionMismatchException expected");
        } catch (DimensionMismatchException ex) {
            // expected
        }
    }

    @Test
    public void testNegativeObservedCounts() {
        final long[] observed = { 0, 1, 2, -3 };
        final double[] expected = { 1, 1, 2, 3};
        final long[] observed2 = {3, 4, 5, 0};
        try {
            testStatistic.gTest(expected, observed);
            Assert.fail("negative observed count, NotPositiveException expected");
        } catch (NotPositiveException ex) {
            // expected
        }
        try {
            testStatistic.gTestDataSetsComparison(observed, observed2);
            Assert.fail("negative observed count, NotPositiveException expected");
        } catch (NotPositiveException ex) {
            // expected
        }
    }

    @Test
    public void testZeroExpectedCounts() {
        final long[] observed = { 0, 1, 2, -3 };
        final double[] expected = { 1, 0, 2, 3};
        try {
            testStatistic.gTest(expected, observed);
            Assert.fail("zero expected count, NotStrictlyPositiveException expected");
        } catch (NotStrictlyPositiveException ex) {
            // expected
        }
    }

    @Test
    public void testBadAlpha() {
        final long[] observed = { 0, 1, 2, 3 };
        final double[] expected = { 1, 2, 2, 3};
        final long[] observed2 = { 0, 2, 2, 3 };
        try {
            testStatistic.gTest(expected, observed, 0.8);
            Assert.fail("zero expected count, NotStrictlyPositiveException expected");
        } catch (OutOfRangeException ex) {
            // expected
        }
        try {
            testStatistic.gTestDataSetsComparison(observed, observed2, -0.5);
            Assert.fail("zero expected count, NotStrictlyPositiveException expected");
        } catch (OutOfRangeException ex) {
            // expected
        }
    }

    @Test
    public void testScaling() {
      final long[] observed = {9, 11, 10, 8, 12};
      final double[] expected1 = {10, 10, 10, 10, 10};
      final double[] expected2 = {1000, 1000, 1000, 1000, 1000};
      final double[] expected3 = {1, 1, 1, 1, 1};
      final double tol = 1E-15;
      Assert.assertEquals(
              testStatistic.gTest(expected1, observed),
              testStatistic.gTest(expected2, observed),
              tol);
      Assert.assertEquals(
              testStatistic.gTest(expected1, observed),
              testStatistic.gTest(expected3, observed),
              tol);
    }

    @Test
    public void testRootLogLikelihood() {
        // positive where k11 is bigger than expected.
        Assert.assertTrue(testStatistic.rootLogLikelihoodRatio(904, 21060, 1144, 283012) > 0.0);

        // negative because k11 is lower than expected
        Assert.assertTrue(testStatistic.rootLogLikelihoodRatio(36, 21928, 60280, 623876) < 0.0);

        Assert.assertEquals(FastMath.sqrt(2.772589), testStatistic.rootLogLikelihoodRatio(1, 0, 0, 1), 0.000001);
        Assert.assertEquals(-FastMath.sqrt(2.772589), testStatistic.rootLogLikelihoodRatio(0, 1, 1, 0), 0.000001);
        Assert.assertEquals(FastMath.sqrt(27.72589), testStatistic.rootLogLikelihoodRatio(10, 0, 0, 10), 0.00001);

        Assert.assertEquals(FastMath.sqrt(39.33052), testStatistic.rootLogLikelihoodRatio(5, 1995, 0, 100000), 0.00001);
        Assert.assertEquals(-FastMath.sqrt(39.33052), testStatistic.rootLogLikelihoodRatio(0, 100000, 5, 1995), 0.00001);

        Assert.assertEquals(FastMath.sqrt(4730.737), testStatistic.rootLogLikelihoodRatio(1000, 1995, 1000, 100000), 0.001);
        Assert.assertEquals(-FastMath.sqrt(4730.737), testStatistic.rootLogLikelihoodRatio(1000, 100000, 1000, 1995), 0.001);

        Assert.assertEquals(FastMath.sqrt(5734.343), testStatistic.rootLogLikelihoodRatio(1000, 1000, 1000, 100000), 0.001);
        Assert.assertEquals(FastMath.sqrt(5714.932), testStatistic.rootLogLikelihoodRatio(1000, 1000, 1000, 99000), 0.001);
    }
}
