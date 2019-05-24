using System;
using System.Collections.Generic;
using System.Runtime.CompilerServices;

namespace Delaunator {

    public class Triangulation {
        private static readonly double Epsilon = Math.Pow(2, -52);
        private List<int> EdgeStack;
        public List<double> coords;
        public List<int> triangles;
        public List<int> halfedges;
        private List<int> hullPrev;
        private List<int> hullNext;
        private List<int> hullTri;
        private int _hashSize;
        private double _cx;
        private double _cy;
        private int hullStart;
        private int trianglesLen;
        public List<int> hull;

        public static Triangulation From<T>(List<T> points, Func<T, double> xExtractor, Func<T, double> yExtractor) {
            int n = points.Count;
            var coords = new List<double>(n * 2);
            for (int i = 0; i < n; i++) {
                var p = points[i];
                coords.Add(xExtractor(p));
                coords.Add(yExtractor(p));
            }
            return new Triangulation(coords);
        }

        public Triangulation(List<double> coords) {
            EdgeStack = new List<int>(512).Fill();
            int n = coords.Count >> 1;
            this.coords = coords;

            // arrays that will store the triangulation graph
            var maxTriangles = 2 * n - 5;
            var triangles = this.triangles = new List<int>(maxTriangles * 3);
            var halfedges = this.halfedges = new List<int>(maxTriangles * 3);

            // temporary arrays for tracking the edges of the advancing convex hull
            this._hashSize = (int)Math.Ceiling(Math.Sqrt(n));
            var hullPrev = this.hullPrev = new List<int>(n); // edge to prev edge
            var hullNext = this.hullNext = new List<int>(n); // edge to next edge
            var hullTri = this.hullTri = new List<int>(n); // edge to adjacent triangle
            var hullHash = new List<int>(this._hashSize).Fill(-1); // angular edge hash

            // populate an array of point indices; calculate input data bbox
            var ids = new List<int>(n);
            double minX = double.PositiveInfinity;
            double minY = double.PositiveInfinity;
            double maxX = double.NegativeInfinity;
            double maxY = double.NegativeInfinity;

            for (int i = 0; i < n; i++) {
                double x = coords[2 * i];
                double y = coords[2 * i + 1];
                if (x < minX) {
                    minX = x;
                }
                if (y < minY) {
                    minY = y;
                }
                if (x > maxX) {
                    maxX = x;
                }
                if (y > maxY) {
                    maxY = y;
                }
                ids.Add(i); //ids.SetSafely(i, i); //ids[i] = i;
            }
            double cx = (minX + maxX) / 2.0;
            double cy = (minY + maxY) / 2.0;

            double minDist = double.PositiveInfinity;
            int i0 = 0, i1 = 0, i2 = 0;

            // pick a seed point close to the center
            for (int i = 0; i < n; i++) {
                double d = Distance(cx, cy, coords[2 * i], coords[2 * i + 1]);
                if (d < minDist) {
                    i0 = i;
                    minDist = d;
                }
            }
            double i0x = coords[2 * i0];
            double i0y = coords[2 * i0 + 1];

            minDist = double.PositiveInfinity;

            // find the point closest to the seed
            for (int i = 0; i < n; i++) {
                if (i == i0) {
                    continue;
                }
                double d = Distance(i0x, i0y, coords[2 * i], coords[2 * i + 1]);
                if (d < minDist && d > 0) {
                    i1 = i;
                    minDist = d;
                }
            }
            double i1x = coords[2 * i1];
            double i1y = coords[2 * i1 + 1];

            double minRadius = double.PositiveInfinity;

            // find the third point which forms the smallest circumcircle with the first two
            for (int i = 0; i < n; i++) {
                if (i == i0 || i == i1) {
                    continue;
                }
                double r = Circumradius(i0x, i0y, i1x, i1y, coords[2 * i], coords[2 * i + 1]);
                if (r < minRadius) {
                    i2 = i;
                    minRadius = r;
                }
            }
            double i2x = coords[2 * i2];
            double i2y = coords[2 * i2 + 1];

            if (minRadius == double.PositiveInfinity) {
                throw new Exception("No Delaunay triangulation exists for this input");
            }

            // swap the order of the seed points for counter-clockwise orientation
            if (Orient(i0x, i0y, i1x, i1y, i2x, i2y)) {
                int i = i1;
                double x = i1x;
                double y = i1y;
                i1 = i2;
                i1x = i2x;
                i1y = i2y;
                i2 = i;
                i2x = x;
                i2y = y;
            }

            ValueTuple<double, double> center = Circumcenter(i0x, i0y, i1x, i1y, i2x, i2y);
            this._cx = center.Item1;
            this._cy = center.Item2;

            List<double> dists = new List<double>(n);
            for (int i = 0; i < n; i++) {
                dists.Add(Distance(coords[2 * i], coords[2 * i + 1], center.Item1, center.Item2));
            }

            // sort the points by distance from the seed triangle circumcenter
            Quicksort(ids, dists, 0, n - 1);

            // set up the seed triangle as the starting hull
            this.hullStart = i0;
            int hullSize = 3;

            //hullNext[i0] = hullPrev[i2] = i1;
            //hullNext[i1] = hullPrev[i0] = i2;
            //hullNext[i2] = hullPrev[i1] = i0;
            hullNext.SetSafely(i0, i1);
            hullNext.SetSafely(i1, i2);
            hullNext.SetSafely(i2, i0);
            hullPrev.SetSafely(i2, i1);
            hullPrev.SetSafely(i0, i2);
            hullPrev.SetSafely(i1, i0);

            //hullTri[i0] = 0;
            //hullTri[i1] = 1;
            //hullTri[i2] = 2;
            hullTri.SetSafely(i0, 0);
            hullTri.SetSafely(i1, 1);
            hullTri.SetSafely(i2, 2);

            hullHash[this.HashKey(i0x, i0y)] = i0;
            hullHash[this.HashKey(i1x, i1y)] = i1;
            hullHash[this.HashKey(i2x, i2y)] = i2;

            this.trianglesLen = 0;
            this.AddTriangle(i0, i1, i2, -1, -1, -1);

            double xp = 0;
            double yp = 0;
            for (int k = 0; k < ids.Count; k++) {
                int i = ids[k];
                double x = coords[2 * i];
                double y = coords[2 * i + 1];

                // skip near-duplicate points
                if (k > 0 && Math.Abs(x - xp) <= Epsilon && Math.Abs(y - yp) <= Epsilon) {
                    continue;
                }
                xp = x;
                yp = y;

                // skip seed triangle points
                if (i == i0 || i == i1 || i == i2) {
                    continue;
                }

                // find a visible edge on the convex hull using edge hash
                int start = 0;
                for (int j = 0, key = this.HashKey(x, y); j < this._hashSize; j++) {
                    start = hullHash[(key + j) % this._hashSize];
                    if (start != -1 && start != hullNext[start])
                        break;
                }

                start = hullPrev[start];
                int e = start;
                int q = hullNext[e];
                while (!Orient(x, y, coords[2 * e], coords[2 * e + 1], coords[2 * q], coords[2 * q + 1])) {
                    e = q;
                    if (e == start) {
                        e = -1;
                        break;
                    }
                    q = hullNext[e];
                }
                if (e == -1) {
                    continue; // likely a near-duplicate point; skip it
                }

                // add the first triangle from the point
                int t = this.AddTriangle(e, i, hullNext[e], -1, -1, hullTri[e]);

                // recursively flip triangles from the point until they satisfy the Delaunay condition
                //hullTri[i] = this.Legalize(t + 2);
                //hullTri[e] = t; // keep track of boundary triangles on the hull
                hullTri.SetSafely(i, this.Legalize(t + 2));
                hullTri.SetSafely(e, t); // keep track of boundary triangles on the hull
                hullSize++;

                // walk forward through the hull, adding more triangles and flipping recursively
                n = hullNext[e];
                q = hullNext[n];
                while (Orient(x, y, coords[2 * n], coords[2 * n + 1], coords[2 * q], coords[2 * q + 1])) {
                    t = this.AddTriangle(n, i, q, hullTri[i], -1, hullTri[n]);
                    //hullTri[i] = this.Legalize(t + 2);
                    //hullNext[n] = n; // mark as removed
                    hullTri.SetSafely(i, this.Legalize(t + 2));
                    hullNext.SetSafely(n, n); // mark as removed
                    hullSize--;
                    n = q;
                    q = hullNext[n];
                }

                // walk backward from the other side, adding more triangles and flipping
                if (e == start) {
                    q = hullPrev[e];
                    while (Orient(x, y, coords[2 * q], coords[2 * q + 1], coords[2 * e], coords[2 * e + 1])) {
                        t = this.AddTriangle(q, i, e, -1, hullTri[e], hullTri[q]);
                        this.Legalize(t + 2);
                        //hullTri[q] = t;
                        //hullNext[e] = e; // mark as removed
                        hullTri.SetSafely(q, t);
                        hullNext.SetSafely(e, e); // mark as removed
                        hullSize--;
                        e = q;
                        q = hullPrev[e];
                    }
                }

                // update the hull indices
                //this.hullStart = hullPrev[i] = e;
                //hullNext[e] = hullPrev[n] = i;
                //hullNext[i] = n;
                this.hullStart = e;
                hullPrev.SetSafely(i, e);
                hullNext.SetSafely(e, i);
                hullPrev.SetSafely(n, i);
                hullNext.SetSafely(i, n);

                // save the two new edges in the hash table
                //hullHash[this.HashKey(x, y)] = i;
                //hullHash[this.HashKey(coords[2 * e], coords[2 * e + 1])] = e;
                hullHash.SetSafely(this.HashKey(x, y), i);
                hullHash.SetSafely(this.HashKey(coords[2 * e], coords[2 * e + 1]), e);
            }

            this.hull = new List<int>(hullSize);
            for (int i = 0, e = this.hullStart; i < hullSize; i++) {
                this.hull.SetSafely(i, e);
                e = hullNext[e];
            }
            this.hullPrev = this.hullNext = this.hullTri = null; // get rid of temporary arrays

            // trim typed triangle mesh arrays
            this.triangles = triangles.GetRange(0, this.trianglesLen);
            this.halfedges = halfedges.GetRange(0, this.trianglesLen);
        }

        private int AddTriangle(int i0, int i1, int i2, int a, int b, int c) {
            int t = this.trianglesLen;

            //this.triangles[t] = i0;
            //this.triangles[t + 1] = i1;
            //this.triangles[t + 2] = i2;
            this.triangles.SetSafely(t, i0);
            this.triangles.SetSafely(t + 1, i1);
            this.triangles.SetSafely(t + 2, i2);

            this.Link(t, a);
            this.Link(t + 1, b);
            this.Link(t + 2, c);

            this.trianglesLen += 3;

            return t;
        }

        private int Legalize(int a) {
            var triangles = this.triangles;
            var coords = this.coords;
            var halfedges = this.halfedges;

            int i = 0;
            int ar = 0;

            // recursion eliminated with a fixed-size stack
            while (true) {
                int b = halfedges[a];

                int a0 = a - a % 3;
                ar = a0 + (a + 2) % 3;

                if (b == -1) { // convex hull edge
                    if (i == 0) {
                        break;
                    }
                    a = EdgeStack[--i];
                    continue;
                }

                int b0 = b - b % 3;
                int al = a0 + (a + 1) % 3;
                int bl = b0 + (b + 2) % 3;

                int p0 = triangles[ar];
                int pr = triangles[a];
                int pl = triangles[al];
                int p1 = triangles[bl];

                bool illegal = InCircle(
                    coords[2 * p0], coords[2 * p0 + 1],
                    coords[2 * pr], coords[2 * pr + 1],
                    coords[2 * pl], coords[2 * pl + 1],
                    coords[2 * p1], coords[2 * p1 + 1]);

                if (illegal) {
                    triangles[a] = p1;
                    triangles[b] = p0;

                    int hbl = halfedges[bl];

                    // edge swapped on the other side of the hull (rare); fix the halfedge reference
                    if (hbl == -1) {
                        int e = this.hullStart;
                        do {
                            if (this.hullTri[e] == bl) {
                                this.hullTri[e] = a;
                                break;
                            }
                            e = this.hullNext[e];
                        } while (e != this.hullStart);
                    }
                    this.Link(a, hbl);
                    this.Link(b, halfedges[ar]);
                    this.Link(ar, bl);

                    int br = b0 + (b + 1) % 3;

                    // don't worry about hitting the cap: it can only happen on extremely degenerate input
                    if (i < EdgeStack.Count) {
                        EdgeStack[i++] = br;
                    }
                }
                else {
                    if (i == 0) {
                        break;
                    }
                    a = EdgeStack[--i];
                }
            }

            return ar;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private void Link(int a, int b) {
            //this.halfedges[a] = b;
            this.halfedges.SetSafely(a, b);
            if (b != -1) {
                //this.halfedges[b] = a;
                this.halfedges.SetSafely(b, a);
            }
        }

        private static double Circumradius(double ax, double ay, double bx, double by, double cx, double cy) {
            double dx = bx - ax;
            double dy = by - ay;
            double ex = cx - ax;
            double ey = cy - ay;

            double bl = dx * dx + dy * dy;
            double cl = ex * ex + ey * ey;
            double d = 0.5 / (dx * ey - dy * ex);

            double x = (ey * bl - dy * cl) * d;
            double y = (dx * cl - ex * bl) * d;

            return x * x + y * y;
        }

        private static ValueTuple<double, double> Circumcenter(double ax, double ay, double bx, double by, double cx, double cy) {
            double dx = bx - ax;
            double dy = by - ay;
            double ex = cx - ax;
            double ey = cy - ay;

            double bl = dx * dx + dy * dy;
            double cl = ex * ex + ey * ey;
            double d = 0.5 / (dx * ey - dy * ex);

            double x = ax + (ey * bl - dy * cl) * d;
            double y = ay + (dx * cl - ex * bl) * d;

            return ValueTuple.Create(x, y);
        }

        private static void Quicksort(List<int> ids, List<double> dists, int left, int right) {
            if (right - left <= 20) {
                for (int i = left + 1; i <= right; i++) {
                    int temp = ids[i];
                    double tempDist = dists[temp];
                    int j = i - 1;
                    while (j >= left && dists[ids[j]] > tempDist) {
                        ids[j + 1] = ids[j--];
                    }
                    ids[j + 1] = temp;
                }
            }
            else {
                int median = (left + right) >> 1;
                int i = left + 1;
                int j = right;
                Swap(ids, median, i);
                if (dists[ids[left]] > dists[ids[right]]) {
                    Swap(ids, left, right);
                }
                if (dists[ids[i]] > dists[ids[right]]) {
                    Swap(ids, i, right);
                }
                if (dists[ids[left]] > dists[ids[i]]) {
                    Swap(ids, left, i);
                }

                int temp = ids[i];
                double tempDist = dists[temp];
                while (true) {
                    do {
                        i++;
                    }
                    while (dists[ids[i]] < tempDist);
                    do {
                        j--;
                    }
                    while (dists[ids[j]] > tempDist);
                    if (j < i) {
                        break;
                    }
                    Swap(ids, i, j);
                }
                ids[left + 1] = ids[j];
                ids[j] = temp;

                if (right - i + 1 >= j - left) {
                    Quicksort(ids, dists, i, right);
                    Quicksort(ids, dists, left, j - 1);
                }
                else {
                    Quicksort(ids, dists, left, j - 1);
                    Quicksort(ids, dists, i, right);
                }
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double Distance(double ax, double ay, double bx, double by) {
            double dx = ax - bx;
            double dy = ay - by;
            return dx * dx + dy * dy;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static bool Orient(double px, double py, double qx, double qy, double rx, double ry) {
            return (qy - py) * (rx - qx) - (qx - px) * (ry - qy) < 0;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static void Swap<T>(List<T> arr, int i, int j) {
            var tmp = arr[i];
            arr[i] = arr[j];
            arr[j] = tmp;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private int HashKey(double x, double y) {
            return (int)Math.Floor(PseudoAngle(x - this._cx, y - this._cy) * this._hashSize) % this._hashSize;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double PseudoAngle(double dx, double dy) {
            double p = dx / (Math.Abs(dx) + Math.Abs(dy));
            return (dy > 0 ? 3 - p : 1 + p) / 4; // [0..1]
        }

        private static bool InCircle(double ax, double ay, double bx, double by, double cx, double cy, double px, double py) {
            double dx = ax - px;
            double dy = ay - py;
            double ex = bx - px;
            double ey = by - py;
            double fx = cx - px;
            double fy = cy - py;

            double ap = dx * dx + dy * dy;
            double bp = ex * ex + ey * ey;
            double cp = fx * fx + fy * fy;

            return dx * (ey * cp - bp * fy) -
                    dy * (ex * cp - bp * fx) +
                    ap * (ex * fy - ey * fx) < 0;
        }
    }
}