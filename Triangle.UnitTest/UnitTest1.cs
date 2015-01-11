using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace Triangle.UnitTest
{
    [TestClass]
    public class UnitTest1
    {
        private static double x11, x12, x21, x22, x31, x32, k1, k2;

        private static double k(double x1, double x2)
        {
            double[] abcab = {x11, x12, x21, x22, x31, x32, x11, x12, x21, x22};
            int i;
            for (i = 0; i < 3; i++)
            {
                double[] v0 = {abcab[2*i + 2] - abcab[2*i], abcab[2*i + 3] - abcab[2*i + 1]};
                double[] v1 = {abcab[2*i + 4] - abcab[2*i], abcab[2*i + 5] - abcab[2*i + 1]};
                double[] v2 = {x1 - abcab[2*i], x2 - abcab[2*i + 1]};
                double sq1 = v0[0]*v1[1] - v0[1]*v1[0];
                double sq2 = v0[0]*v2[1] - v0[1]*v2[0];
                if (sq1*sq2 < 0.0) return k2;
            }
            return k1;
        }

        [TestMethod]
        public void TestTriangleMethod()
        {
            x11 = 0.2;
            x12 = 0.5;
            x21 = 0.7;
            x22 = 0.2;
            x31 = 0.5;
            x32 = 0.8;
            k1 = 25.0;
            k2 = 1.0;

            Assert.IsTrue(Math.Abs(k(0.0, 0.0) - 1.0) < 0.000001);
            Assert.IsTrue(Math.Abs(k(1.0, 1.0) - 1.0) < 0.000001);
            Assert.IsTrue(Math.Abs(k(1.0, 0.0) - 1.0) < 0.000001);
            Assert.IsTrue(Math.Abs(k(0.0, 1.0) - 1.0) < 0.000001);
            Assert.IsTrue(Math.Abs(k(0.5, 0.5) - 25.0) < 0.000001);
        }
    }
}