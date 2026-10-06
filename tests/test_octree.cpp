// Octree correctness tests. Every test compares against a brute-force
// reference and reports accuracy to the terminal.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <set>
#include <string>
#include <vector>

#include "octree.h"

using octree::AABB;
using octree::Octree;
using octree::Point3D;

namespace {

std::mt19937_64 rng(20261003);

double randUniform(double lo, double hi) {
    std::uniform_real_distribution<double> d(lo, hi);
    return d(rng);
}

Point3D randPoint(double lo, double hi) {
    return {randUniform(lo, hi), randUniform(lo, hi), randUniform(lo, hi)};
}

// Canonical key so multisets of points can be compared exactly.
std::multiset<std::string> keysOf(const std::vector<Point3D>& pts) {
    std::multiset<std::string> s;
    char buf[128];
    for (const auto& p : pts) {
        std::snprintf(buf, sizeof(buf), "%.17g|%.17g|%.17g", p.x, p.y, p.z);
        s.insert(buf);
    }
    return s;
}

// ---- brute-force references -------------------------------------------------

std::vector<Point3D> bfRange(const std::vector<Point3D>& pts, const AABB& b) {
    std::vector<Point3D> out;
    for (const auto& p : pts)
        if (b.contains(p)) out.push_back(p);
    return out;
}

std::vector<Point3D> bfRadius(const std::vector<Point3D>& pts,
                              const Point3D& c, double r) {
    std::vector<Point3D> out;
    const double rSq = r * r;
    for (const auto& p : pts)
        if (octree::distanceSq(p, c) <= rSq) out.push_back(p);
    return out;
}

bool bfNearest(const std::vector<Point3D>& pts, const Point3D& q,
               Point3D& out) {
    if (pts.empty()) return false;
    double best = INFINITY;
    for (const auto& p : pts) {
        const double d = octree::distanceSq(p, q);
        if (d < best) {
            best = d;
            out = p;
        }
    }
    return true;
}

// ---- test harness -----------------------------------------------------------

int g_pass = 0, g_fail = 0;

void report(const std::string& name, size_t correct, size_t total) {
    const bool ok = correct == total;
    if (ok) ++g_pass; else ++g_fail;
    std::printf("[%s] %-34s accuracy: %zu/%zu (%.2f%%)\n",
                ok ? "PASS" : "FAIL", name.c_str(), correct, total,
                total ? 100.0 * correct / total : 100.0);
}

void expectTrue(const std::string& name, bool cond) {
    report(name, cond ? 1 : 0, 1);
}

// ---- tests ------------------------------------------------------------------

void testRandomCloud() {
    const size_t kPoints = 5000, kQueries = 200;
    std::vector<Point3D> pts;
    pts.reserve(kPoints);
    Octree tree(8, 10);
    for (size_t i = 0; i < kPoints; ++i) {
        Point3D p = randPoint(-100.0, 100.0);
        pts.push_back(p);
        tree.insert(p);
    }
    expectTrue("random: size", tree.size() == kPoints);

    size_t ok = 0;
    for (size_t i = 0; i < kQueries; ++i) {
        Point3D a = randPoint(-120.0, 60.0);
        Point3D b = randPoint(a.x, a.x + 80.0);
        AABB box{{a.x, a.y, a.z}, {b.x, b.y, b.z}};
        // normalize (a may exceed b componentwise after randPoint above)
        box.min.x = std::min(a.x, b.x); box.max.x = std::max(a.x, b.x);
        box.min.y = std::min(a.y, b.y); box.max.y = std::max(a.y, b.y);
        box.min.z = std::min(a.z, b.z); box.max.z = std::max(a.z, b.z);
        if (keysOf(tree.rangeQuery(box)) == keysOf(bfRange(pts, box))) ++ok;
    }
    report("random: range query", ok, kQueries);

    ok = 0;
    for (size_t i = 0; i < kQueries; ++i) {
        Point3D c = randPoint(-100.0, 100.0);
        const double r = randUniform(0.0, 60.0);
        if (keysOf(tree.radiusSearch(c, r)) == keysOf(bfRadius(pts, c, r)))
            ++ok;
    }
    report("random: radius search", ok, kQueries);

    ok = 0;
    for (size_t i = 0; i < kQueries; ++i) {
        Point3D q = randPoint(-150.0, 150.0);
        Point3D got, want;
        const bool found = tree.nearestNeighbor(q, got);
        bfNearest(pts, q, want);
        if (found && octree::distanceSq(got, q) ==
                         octree::distanceSq(want, q))
            ++ok;
    }
    report("random: nearest neighbor", ok, kQueries);
}

void testBoundaryPoints() {
    // Points exactly on split planes, box edges and corners.
    Octree tree(4, 6, nullptr);
    std::vector<Point3D> pts;
    // First point defines the root cube [-0.5,0.5]^3 around origin.
    std::vector<Point3D> seeds = {
        {0, 0, 0},       {0.5, 0.5, 0.5},   {-0.5, -0.5, -0.5},
        {0.25, 0, 0},    {0, 0.25, 0},      {0, 0, 0.25},
        {0.25, 0.25, 0}, {0.25, 0.25, 0.25}, {-0.25, 0.25, -0.25},
        {0.5, 0, 0},     {0, 0.5, 0.5},     {-0.5, 0.5, -0.5},
    };
    for (const auto& p : seeds) {
        pts.push_back(p);
        tree.insert(p);
    }
    // Query box exactly matching root bounds must return everything.
    AABB rootBox = tree.bounds();
    expectTrue("boundary: full-box query",
               keysOf(tree.rangeQuery(rootBox)) == keysOf(pts));
    // Degenerate box on a split plane.
    AABB plane{{0.25, -1, -1}, {0.25, 1, 1}};
    expectTrue("boundary: split-plane box",
               keysOf(tree.rangeQuery(plane)) == keysOf(bfRange(pts, plane)));
    // Zero-size box exactly on a point.
    AABB dot{{0.25, 0.25, 0.25}, {0.25, 0.25, 0.25}};
    expectTrue("boundary: zero-size box",
               keysOf(tree.rangeQuery(dot)) == keysOf(bfRange(pts, dot)));
    // Radius search with r=0 centered on a point.
    expectTrue("boundary: zero radius",
               keysOf(tree.radiusSearch({0.25, 0.25, 0.25}, 0.0)) ==
                   keysOf(bfRadius(pts, {0.25, 0.25, 0.25}, 0.0)));
    // Sphere tangent to split planes.
    expectTrue("boundary: tangent sphere",
               keysOf(tree.radiusSearch({0, 0, 0}, 0.25)) ==
                   keysOf(bfRadius(pts, {0, 0, 0}, 0.25)));
}

void testDuplicatePoints() {
    Octree tree(2, 8);
    std::vector<Point3D> pts;
    const Point3D dup{1.0, 2.0, 3.0};
    for (int i = 0; i < 50; ++i) {
        pts.push_back(dup);
        tree.insert(dup);
    }
    for (int i = 0; i < 30; ++i) {
        Point3D p = randPoint(-5, 5);
        pts.push_back(p);
        tree.insert(p);
    }
    expectTrue("duplicates: size", tree.size() == pts.size());
    AABB box{{0.5, 1.5, 2.5}, {1.5, 2.5, 3.5}};
    expectTrue("duplicates: range keeps all copies",
               keysOf(tree.rangeQuery(box)) == keysOf(bfRange(pts, box)));
    expectTrue("duplicates: radius keeps all copies",
               keysOf(tree.radiusSearch(dup, 0.0)) ==
                   keysOf(bfRadius(pts, dup, 0.0)));
    Point3D nn;
    expectTrue("duplicates: nearest is the duplicate",
               tree.nearestNeighbor({1.001, 2.0, 3.0}, nn) && nn == dup);
}

void testEmptyTree() {
    Octree tree;
    AABB box{{-1, -1, -1}, {1, 1, 1}};
    expectTrue("empty: range query empty", tree.rangeQuery(box).empty());
    expectTrue("empty: radius search empty",
               tree.radiusSearch({0, 0, 0}, 10.0).empty());
    Point3D out;
    expectTrue("empty: nearest returns false",
               !tree.nearestNeighbor({0, 0, 0}, out));
    expectTrue("empty: size is zero", tree.size() == 0);
}

void testDepthLimit() {
    const size_t kMaxDepth = 4;
    const size_t kPoints = 3000;
    Octree tree(1, kMaxDepth);  // capacity 1 forces maximal splitting
    std::vector<Point3D> pts;
    for (size_t i = 0; i < kPoints; ++i) {
        Point3D p = randPoint(-1.0, 1.0);
        pts.push_back(p);
        tree.insert(p);
    }
    // Even with a shallow depth limit every point must still be found.
    size_t ok = 0;
    const size_t kQueries = 100;
    for (size_t i = 0; i < kQueries; ++i) {
        Point3D c = randPoint(-1.0, 1.0);
        const double r = randUniform(0.1, 1.5);
        if (keysOf(tree.radiusSearch(c, r)) == keysOf(bfRadius(pts, c, r)))
            ++ok;
    }
    report("depth-limit: radius search", ok, kQueries);
    ok = 0;
    for (size_t i = 0; i < kQueries; ++i) {
        Point3D q = randPoint(-2.0, 2.0);
        Point3D got, want;
        tree.nearestNeighbor(q, got);
        bfNearest(pts, q, want);
        if (octree::distanceSq(got, q) == octree::distanceSq(want, q)) ++ok;
    }
    report("depth-limit: nearest neighbor", ok, kQueries);
}

void testPersistence(const std::string& tmpPath) {
    std::vector<Point3D> pts;
    Octree tree(8, 10);
    for (int i = 0; i < 2000; ++i) {
        Point3D p = randPoint(-50.0, 50.0);
        pts.push_back(p);
        tree.insert(p);
    }
    // add duplicates and boundary points
    for (int i = 0; i < 20; ++i) {
        pts.push_back({7.0, 7.0, 7.0});
        tree.insert({7.0, 7.0, 7.0});
    }
    expectTrue("persist: save succeeds", tree.save(tmpPath));

    Octree loaded;
    expectTrue("persist: load succeeds", loaded.load(tmpPath));
    expectTrue("persist: size matches", loaded.size() == tree.size());

    size_t ok = 0;
    const size_t kQueries = 100;
    for (size_t i = 0; i < kQueries; ++i) {
        Point3D a = randPoint(-60, 0), b = randPoint(0, 60);
        AABB box{{a.x, a.y, a.z}, {b.x, b.y, b.z}};
        if (keysOf(loaded.rangeQuery(box)) == keysOf(tree.rangeQuery(box)))
            ++ok;
    }
    report("persist: range identical after load", ok, kQueries);

    ok = 0;
    for (size_t i = 0; i < kQueries; ++i) {
        Point3D c = randPoint(-50, 50);
        const double r = randUniform(1, 40);
        if (keysOf(loaded.radiusSearch(c, r)) ==
            keysOf(tree.radiusSearch(c, r)))
            ++ok;
    }
    report("persist: radius identical after load", ok, kQueries);

    ok = 0;
    for (size_t i = 0; i < kQueries; ++i) {
        Point3D q = randPoint(-80, 80), g1, g2;
        loaded.nearestNeighbor(q, g1);
        tree.nearestNeighbor(q, g2);
        if (octree::distanceSq(g1, q) == octree::distanceSq(g2, q)) ++ok;
    }
    report("persist: nearest identical after load", ok, kQueries);

    // Loading garbage must fail cleanly.
    expectTrue("persist: reject invalid file", !loaded.load("/dev/null"));
}

void testSingleAndFarPoints() {
    Octree tree;
    tree.insert({1e6, -1e6, 1e6});  // far from origin, root must expand
    tree.insert({-1e6, 1e6, -1e6});
    tree.insert({0, 0, 0});
    Point3D nn;
    expectTrue("far: nearest among 3",
               tree.nearestNeighbor({1e6, -1e6, 1e6 + 1}, nn) &&
                   nn == (Point3D{1e6, -1e6, 1e6}));
    AABB all{{-2e6, -2e6, -2e6}, {2e6, 2e6, 2e6}};
    expectTrue("far: range finds all 3", tree.rangeQuery(all).size() == 3);
}

}  // namespace

int main() {
    testRandomCloud();
    testBoundaryPoints();
    testDuplicatePoints();
    testEmptyTree();
    testDepthLimit();
    testPersistence("/tmp/octree_test_dump.bin");
    testSingleAndFarPoints();

    const int total = g_pass + g_fail;
    std::printf("\n==== summary: %d/%d test groups passed (%.2f%%) ====\n",
                g_pass, total, 100.0 * g_pass / total);
    return g_fail == 0 ? 0 : 1;
}
