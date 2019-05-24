using System;
using System.Collections.Generic;
using System.IO;
using Microsoft.VisualStudio.TestTools.UnitTesting;

internal struct Vector2 {
    public double X;
    public double Y;
    public Vector2(double x, double y) {
        this.X = x;
        this.Y = y;
    }
}

[TestClass()]
public class Tests {

    [TestMethod()]
    public void TriangulatesArray() {
        var points = LoadDefaultPoints();
        var concatenated = new List<double>(points.Count * 2);
        points.ForEach(p => {
            concatenated.Add(p.X);
            concatenated.Add(p.Y);
        });
        var d1 = new Delaunator.Triangulation(concatenated);
        var d2 = Delaunator.Triangulation.From(points, p => p.X, p => p.Y);
        CollectionAssert.AreEqual(d1.triangles, d2.triangles);
    }

    [TestMethod()]
    public void ProducesCorrectTriangulation() {
        Validate(LoadDefaultPoints());
    }

    [TestMethod()]
    public void TestIssue11() {
        Validate(new List<Vector2>() {
            new Vector2(516, 661),
            new Vector2(369, 793),
            new Vector2(426, 539),
            new Vector2(273, 525),
            new Vector2(204, 694),
            new Vector2(747, 750),
            new Vector2(454, 390)
        });
    }

    [TestMethod()]
    public void TestIssue24() {
        Validate(new List<Vector2>() {
            new Vector2(382, 302),
            new Vector2(382, 328),
            new Vector2(382, 205),
            new Vector2(623, 175),
            new Vector2(382, 188),
            new Vector2(382, 284),
            new Vector2(623, 87),
            new Vector2(623, 341),
            new Vector2(141, 227)
        });
    }

    [TestMethod()]
    public void TestIssue13() {
        Validate(new List<Vector2>() {
            new Vector2(4, 1),
            new Vector2(3.7974166882130675, 2.0837249985614585),
            new Vector2(3.2170267516619773, 3.0210869309396715),
            new Vector2(2.337215067329615, 3.685489874065187),
            new Vector2(1.276805078389906, 3.9872025288851036),
            new Vector2(0.17901102978375127, 3.885476929518457),
            new Vector2(-0.8079039091377689, 3.3940516818407187),
            new Vector2(-1.550651407188842, 2.5792964886320684),
            new Vector2(-1.9489192990517052, 1.5512485534497125),
            new Vector2(-1.9489192990517057, 0.44875144655029087),
            new Vector2(-1.5506514071888438, -0.5792964886320653),
            new Vector2(-0.8079039091377715, -1.394051681840717),
            new Vector2(0.17901102978374794, -1.8854769295184561),
            new Vector2(1.276805078389902, -1.987202528885104),
            new Vector2(2.337215067329611, -1.6854898740651891),
            new Vector2(3.217026751661974, -1.021086930939675),
            new Vector2(3.7974166882130653, -0.08372499856146409)
        });
    }

    [TestMethod()]
    public void TestRobustness() {
        var points = LoadPoints("robustness1.txt", 79);
        Validate(points);
        Validate(points.ConvertAll(p => new Vector2(p.X / 1e9, p.Y / 1e9)));
        Validate(points.ConvertAll(p => new Vector2(p.X / 100, p.Y / 100)));
        Validate(points.ConvertAll(p => new Vector2(p.X * 100, p.Y * 100)));
        Validate(points.ConvertAll(p => new Vector2(p.X * 1e9, p.Y * 1e9)));
        Validate(LoadPoints("robustness2.txt", 100));
        Validate(LoadPoints("robustness3.txt", 1000));
    }

    [TestMethod()]
    public void ThrowsOnSmallNumberOfPoints() {
        var points = LoadDefaultPoints();
        try {
            Delaunator.Triangulation.From(points.GetRange(0, 1), v => v.X, v => v.Y);
            Assert.Fail("Should have caught");
        }
        catch (Exception) {
        }
        try {
            Delaunator.Triangulation.From(points.GetRange(0, 2), v => v.X, v => v.Y);
            Assert.Fail("Should have caught");
        }
        catch (Exception) {
        }
    }

    [TestMethod()]
    public void ThrowsOnAllColinearPoints() {
        try {
            Delaunator.Triangulation.From(new List<Vector2>() {
                new Vector2(0, 0),
                new Vector2(1, 0),
                new Vector2(2, 0),
                new Vector2(3, 0)
            }, v => v.X, v => v.Y);
            Assert.Fail("Should have caught");
        }
        catch (Exception) {
        }
    }

    [TestMethod()]
    public void SupportsCustomPointFormat() {
        Func<string, double> x = (string s) => {
            double.TryParse(s.Split(',')[0], out double result);
            return result;
        };
        Func<string, double> y = (string s) => {
            double.TryParse(s.Split(',')[1], out double result);
            return result;
        };
        Delaunator.Triangulation.From(new List<string>() {
            "5,5", "7,5", "7,6"
        }, v => x(v), v => y(v));
    }

    private static void Validate(List<Vector2> points) {
        var d = Delaunator.Triangulation.From(points, v => v.X, v => v.Y);

        // validate halfedges
        for (int i = 0; i < d.halfedges.Count; i++) {
            int i2 = d.halfedges[i];
            if (i2 != -1 && d.halfedges[i2] != i) {
                Assert.Fail("invalid halfedge connection");
            }
        }

        // validate triangulation
        List<double> hullAreas = new List<double>();
        for (int i = 0, len = d.hull.Count, j = len - 1; i < len; j = i++) {
            Vector2 p0 = points[d.hull[j]];
            Vector2 p = points[d.hull[i]];
            hullAreas.Add((p.X - p0.X) * (p.Y + p0.Y));
        }
        double hullArea = Sum(hullAreas);

        List<double> triangleAreas = new List<double>();
        for (int i = 0; i < d.triangles.Count; i += 3) {
            Vector2 a = points[d.triangles[i]];
            Vector2 b = points[d.triangles[i + 1]];
            Vector2 c = points[d.triangles[i + 2]];
            triangleAreas.Add(Math.Abs((b.Y - a.Y) * (c.X - b.X) - (b.X - a.X) * (c.Y - b.Y)));
        }
        double trianglesArea = Sum(triangleAreas);

        double err = Math.Abs((hullArea - trianglesArea) / hullArea);
        if (err <= Math.Pow(2, -51)) {
        }
        else {
            Assert.Fail("triangulation is broken: " + err + " error");
        }
    }

    // Kahan and Babuska summation, Neumaier variant; accumulates less FP error
    public static double Sum(List<double> x) {
        double sum = x[0];
        double err = 0;
        for (int i = 1; i < x.Count; i++) {
            double k = x[i];
            double m = sum + k;
            err += Math.Abs(sum) >= Math.Abs(k) ? sum - m + k : k - m + sum;
            sum = m;
        }
        return sum + err;
    }

    private static List<Vector2> LoadDefaultPoints() {
        return LoadPoints("points.txt", 874);
    }

    private static List<Vector2> LoadPoints(string file, int count) {
        List<Vector2> result = new List<Vector2>();
        string[] lines = File.ReadAllLines("../../" + file);
        foreach (string line in lines) {
            var parts = line.Split(',');
            if (parts.Length != 2) {
                continue;
            }
            string xString = parts[0].Trim().Replace("[", "");
            string yString = parts[1].Trim().Replace("]", "");
            if (!double.TryParse(xString, out double x)) {
                continue;
            }
            if (!double.TryParse(yString, out double y)) {
                continue;
            }
            result.Add(new Vector2(x, y));
        }
        if (result.Count != count) {
            Assert.AreEqual(count, result.Count);
        }
        return result;
    }
}
