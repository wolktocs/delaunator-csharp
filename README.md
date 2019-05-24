# delaunator-csharp
A C# library for Delaunay triangulation of 2D points.  A port of [Delaunator](https://github.com/mapbox/delaunator).

While looking for a library to help me with experiments in generating Delaunay triangulations and their accompanying Voronoi diagrams, I came across some excellent Javascript projects that made use of Delaunator.  I was a little unhappy with the available C# options, so I decided to do a port.

## Building

The solution and projects were created with Visual Studio 2019.  The projects include:
 + *Delaunator*: the library itself, builds a DLL
 + *DelaunatorBenchmarks*: a console application that runs a simple benchmark
 + *DelaunatorTests*: a unit test class that makes use of the MSTest framework
 
## Usage

The library closely follows the original implementation.  The constructor for the `Delaunator` class accepts a `List` of `double`'s where every two elements of the list are an `x` and `y` coordinate of the input point set.

```
using Delaunator;
using System;
using System.Collections.Generic;

public class Program {
  static void Main(string[] args) {
    var inputPoints = new List<double>() { 
      66.103648384371410, 68.588612471664760,
      146.680713462100413, 121.680713462100428,
      128.868896560467447, 117.261797559041411,
      66.103648384371439, 68.588612471664774,
      169.552139667571992, 146.133776538276890,
      126.629392246050883, 181.111404660392082,
      74.434448280233709, 78.630898779520691,
      121.111404660392054, 153.370607753949116,
      98.888595339607888, 186.629392246050855,
      52.660668968140221, 63.178539267712423
    }
    Delaunator d = new Delaunator(inputPoints);
  }
}
```

Alternately, you can call the static `Delaunator.From()` method with a `List` of a custom objects and provide functions that extract the x- and y- coordinates.
```
using Delaunator;
using System;
using System.Collections.Generic;

public class Program {
  public struct Vector {
    public double X;
    public double Y;
    public Vector(double x, double y) {
      X = x;
      Y = y;
    }
  }
  static void Main(string[] args) {
    var inputPoints = new List<Vector>() { 
      new Vector(66.103648384371410, 68.588612471664760),
      new Vector(146.680713462100413, 121.680713462100428),
      new Vector(128.868896560467447, 117.261797559041411),
      new Vector(66.103648384371439, 68.588612471664774),
      new Vector(169.552139667571992, 146.133776538276890),
      new Vector(126.629392246050883, 181.111404660392082),
      new Vector(74.434448280233709, 78.630898779520691),
      new Vector(121.111404660392054, 153.370607753949116),
      new Vector(98.888595339607888, 186.629392246050855),
      new Vector(52.660668968140221, 63.178539267712423)
    };
    Delaunator d = Delaunator.From(inputPoints,
      p => p.X, // extracts the x-coordinate from a Vector
      p => p.Y // extracts the y-coordinate from a Vector
    );
  }
}
```

## Performance

The C# port does not perform as well the original Javascript version running with nodejs, and I'm not entirely sure why.  However, it does seem to run quite a bit faster than the C# alternatives that I tried out before doing the port.

&nbsp; | uniform 100k | gauss 100k | grid 100k | degen 100k | uniform 1&nbsp;million | gauss 1&nbsp;million | grid 1&nbsp;million | degen 1&nbsp;million
:-- | --: | --: | --: | --: | --: | --: | --: | --:
[delaunator](https://github.com/mapbox/delaunator) | 75ms | 74ms | 117ms | 30ms | 1.36s | 1.26s | 1.07s | 323ms
**delaunator-csharp** | 110ms | 107ms | 90ms | 42ms | 1.65s | 1.63s | 1.30s | 469ms
[Triangle.NET](https://github.com/eppz/Triangle.NET) | 340ms | 350ms | 545ms | 500ms | 4.10s | 4.10s | 5.00s | 5.03s
[MIConvexHull](https://github.com/DesignEngrLab/MIConvexHull) | 1.23s | 1.24s | 1.03s | 532ms | 15.80s | 16.23s | 11.01s | 5.26s

**Note:** all of these libraries typically require you to convert your point/vector data structures into corresponding data structures used by the libraries.  I believe that real-world usage will typically require those conversions, so the benchmark timings include the conversion time.

Benchmarks were run on Windows 10 on an x64-based processor, Intel Core i7-4790S CPU @ 3.20GHz, 16 GB of memory.
