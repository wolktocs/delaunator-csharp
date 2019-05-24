using System;
using System.Collections.Generic;

public struct Vector {
    public double x;
    public double y;
    public Vector(double x, double y) {
        this.x = x;
        this.y = y;
    }
}

public class Benchmark {

    private Random random = new Random();
    private System.Diagnostics.Stopwatch stopwatch = new System.Diagnostics.Stopwatch();

    static void Main(string[] args) {
        Benchmark program = new Benchmark();
        var Distributions = new List<Tuple<string, Func<int, List<Vector>>>>() {
            new Tuple<string, Func<int, List<Vector>>>("Uniform", program.Uniform),
            new Tuple<string, Func<int, List<Vector>>>("Gaussian", program.Gaussian),
            new Tuple<string, Func<int, List<Vector>>>("Grid", program.Grid),
            new Tuple<string, Func<int, List<Vector>>>("Degenerate", program.Degenerate)
        };
        List<int> Counts = new List<int>() { 20000, 100000, 200000, 500000, 1000000 };

        // The project by default does not add references to other Delaunay triangulation
        // libraries for benchmarking, but you can add those references yourself and 
        // uncomment the code below accordingly.

        // Benchmark Delaunator
        Action<List<Vector>, int, bool> triangulate = program.TriangulateDelaunator;

        // Benchmark Triangle.NET
        //Action<List<Vector>, int, bool> triangulate = program.TriangulateWithTriangleNet;

        // Benchmark MLConvexHull
        //Action<List<Vector>, int, bool> triangulate = program.TriangulateWithMIConvexHull;

        foreach (var generatorPair in Distributions) {
            // warmup
            triangulate(generatorPair.Item2(Counts[0]), Counts[0], false);
            triangulate(generatorPair.Item2(Counts[1]), Counts[1], false);

            Console.Out.WriteLine(string.Format("{0}:", generatorPair.Item1));
            for (var i = 0; i < Counts.Count; i++) {
                var c = Counts[i];
                var points = generatorPair.Item2(c);
                triangulate(points, c, true);
            }
        }

        Console.Out.WriteLine("Press return/enter to finish");
        Console.In.ReadLine();
    }

    private void WriteResult(int count, System.Diagnostics.Stopwatch stopwatch) {
        Console.Out.WriteLine(string.Format("{0,10:N0}: {1,8:N0}ms", count, stopwatch.ElapsedMilliseconds));
    }

    public void TriangulateDelaunator(List<Vector> points, int count, bool showResult = true) {
        stopwatch.Restart();
        var d = Delaunator.Triangulation.From(points, (Vector v) => { return v.x; }, (Vector v) => { return v.y; });
        stopwatch.Stop();
        if (showResult) {
            WriteResult(count, stopwatch);
        }
    }

    // To benchmark Triangle.NET, uncomment and add the necessary reference
    //public void TriangulateWithTriangleNet(List<Vector> points, int count, bool showResult = true) {
    //    var input = new TriangleNet.Geometry.Polygon(points.Count);
    //    points.ForEach(p => {
    //        input.Add(new TriangleNet.Geometry.Vertex(p.x, p.y));
    //    });
    //    var options = new TriangleNet.Meshing.ConstraintOptions();
    //    var quality = new TriangleNet.Meshing.QualityOptions();
    //    stopwatch.Restart();
    //    var delaunay = TriangleNet.Geometry.ExtensionMethods.Triangulate(input, options, quality);
    //    stopwatch.Stop();
    //    if (showResult) {
    //        WriteResult(count, stopwatch);
    //    }
    //}

    // To benchmark MIConvexHull, uncomment and add the necessary reference
    //public void TriangulateWithMIConvexHull(List<Vector> points, int count, bool showResult = true) {
    //    var input = points.ConvertAll(p => new MIConvexHullVertex2(p.x, p.y));
    //    stopwatch.Restart();
    //    var result = MIConvexHull.VoronoiMesh.Create<MIConvexHullVertex2, MIConvexHullCell2>(input);
    //    stopwatch.Stop();
    //    if (showResult) {
    //        WriteResult(count, stopwatch);
    //    }
    //}
    //private class MIConvexHullFace2 : MIConvexHull.ConvexFace<MIConvexHullVertex2, MIConvexHullFace2> { }
    //private class MIConvexHullCell2 : MIConvexHull.TriangulationCell<MIConvexHullVertex2, MIConvexHullCell2> { }
    //private struct MIConvexHullVertex2 : MIConvexHull.IVertex {
    //    public MIConvexHullVertex2(double x, double y) {
    //        Position = new double[] { x, y };
    //    }
    //    public double[] Position {
    //        get;
    //        set;
    //    }
    //}

    private List<Vector> Uniform(int count) {
        List<Vector> points = new List<Vector>(count);
        for (int i=0; i< count; i++) {
            points.Add(new Vector(random.NextDouble() * 1e3, random.NextDouble() * 1e3));
        }
        return points;
    }

    private List<Vector> Grid(int count) {
        List<Vector> points = new List<Vector>(count);
        int size = (int)Math.Sqrt(count);
        for (var i = 0; i<size; i++) {
            for (var j = 0; j < size; j++) {
                points.Add(new Vector(i, j));
            }
        }
        return points;
    }

    private List<Vector> Gaussian(int count) {
        List<Vector> points = new List<Vector>(count);
        for (int i = 0; i < count; i++) {
            points.Add(new Vector(PseudoNormal() * 1e3, PseudoNormal() * 1e3));
        }
        return points;
    }

    private List<Vector> Degenerate(int count) {
        List<Vector> points = new List<Vector>() { new Vector(0, 0) };
        for (int i = 0; i < count; i++) {
            var angle = 2.0 * Math.PI * (double)i / (double)count;
            points.Add(new Vector(1e10 * Math.Sin(angle), 1e10 * Math.Cos(angle)));
        }
        return points;
    }

    private double PseudoNormal() {
        var v = random.NextDouble() + random.NextDouble() + random.NextDouble() + random.NextDouble() + random.NextDouble() + random.NextDouble();
        return Math.Min(0.5 * (v - 3.0) / 3.0, 1.0);
    }

}
