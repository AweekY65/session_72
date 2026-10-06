#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

#include "octree.h"

namespace {

using octree::Point3d;

int g_checks = 0;
int g_failures = 0;

void report(const std::string& name, std::size_t matched, std::size_t total) {
    ++g_checks;
    const double accuracy = total == 0 ? 100.0 : 100.0 * matched / total;
    if (matched != total) ++g_failures;
    std::printf("  [%s] %-44s %zu/%zu  accuracy = %.2f%%\n",
                matched == total ? "PASS" : "FAIL", name.c_str(), matched,
                total, accuracy);
}

void expect(bool condition, const std::string& name) {
    report(name, condition ? 1 : 0, 1);
}

// ---------- brute-force reference implementations ----------

std::vector<Point3d> bruteRange(const std::vector<Point3d>& pts,
                                const Point3d& lo, const Point3d& hi) {
    std::vector<Point3d> out;
    for (const Point3d& p : pts) {
        if (p.x >= lo.x && p.x <= hi.x && p.y >= lo.y && p.y <= hi.y &&
            p.z >= lo.z && p.z <= hi.z)
            out.push_back(p);
    }
    return out;
}

std::vector<Point3d> bruteRadius(const std::vector<Point3d>& pts,
                                 const Point3d& c, double r) {
    std::vector<Point3d> out;
    const double r2 = r * r;
    for (const Point3d& p : pts) {
        if (octree::squaredDistance(c, p) <= r2) out.push_back(p);
    }
    return out;
}

bool bruteNearest(const std::vector<Point3d>& pts, const Point3d& q,
                  double& bestDist) {
    if (pts.empty()) return false;
    double best = HUGE_VAL;
    for (const Point3d& p : pts) best = std::min(best, octree::squaredDistance(q, p));
    bestDist = std::sqrt(best);
    return true;
}

// Multiset comparison of two point clouds (exact double equality).
bool sameMultiset(std::vector<Point3d> a, std::vector<Point3d> b) {
    if (a.size() != b.size()) return false;
    auto cmp = [](const Point3d& u, const Point3d& v) {
        if (u.x != v.x) return u.x < v.x;
        if (u.y != v.y) return u.y < v.y;
        return u.z < v.z;
    };
    std::sort(a.begin(), a.end(), cmp);
    std::sort(b.begin(), b.end(), cmp);
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a[i].x != b[i].x || a[i].y != b[i].y || a[i].z != b[i].z)
            return false;
    }
    return true;
}

struct Rng {
    std::mt19937_64 gen;
    explicit Rng(std::uint64_t seed) : gen(seed) {}
    double uniform(double lo, double hi) {
        return std::uniform_real_distribution<double>(lo, hi)(gen);
    }
};

// ---------- test cases ----------

void testRandomCloud() {
    std::printf("[test] random point cloud vs brute force\n");
    Rng rng(42);
    std::vector<Point3d> pts;
    const std::size_t n = 5000;
    pts.reserve(n);
    octree::Octree tree({-100, -100, -100}, {100, 100, 100}, 8, 16);
    for (std::size_t i = 0; i < n; ++i) {
        Point3d p{rng.uniform(-100, 100), rng.uniform(-100, 100),
                  rng.uniform(-100, 100)};
        pts.push_back(p);
        tree.insert(p);
    }
    expect(tree.size() == n, "insert count matches");

    std::size_t ok = 0, total = 0;
    for (int i = 0; i < 200; ++i) {
        Point3d a{rng.uniform(-100, 100), rng.uniform(-100, 100),
                  rng.uniform(-100, 100)};
        Point3d b{rng.uniform(-100, 100), rng.uniform(-100, 100),
                  rng.uniform(-100, 100)};
        Point3d lo{std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z)};
        Point3d hi{std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z)};
        ++total;
        ok += sameMultiset(tree.rangeQuery(lo, hi), bruteRange(pts, lo, hi));
    }
    report("random range queries", ok, total);

    ok = 0; total = 0;
    for (int i = 0; i < 200; ++i) {
        Point3d c{rng.uniform(-100, 100), rng.uniform(-100, 100),
                  rng.uniform(-100, 100)};
        const double r = rng.uniform(0, 60);
        ++total;
        ok += sameMultiset(tree.radiusSearch(c, r), bruteRadius(pts, c, r));
    }
    report("random radius searches", ok, total);

    ok = 0; total = 0;
    for (int i = 0; i < 500; ++i) {
        Point3d q{rng.uniform(-120, 120), rng.uniform(-120, 120),
                  rng.uniform(-120, 120)};
        Point3d nearest;
        double dist = 0.0, refDist = 0.0;
        ++total;
        if (tree.nearestNeighbor(q, nearest, &dist) &&
            bruteNearest(pts, q, refDist) && std::fabs(dist - refDist) < 1e-9)
            ++ok;
    }
    report("random nearest neighbors", ok, total);
}

void testBoundaryPoints() {
    std::printf("[test] split-plane and boundary points\n");
    octree::Octree tree({0, 0, 0}, {8, 8, 8}, 2, 8);
    std::vector<Point3d> pts = {
        {4, 4, 4},  // exactly on the root split plane
        {2, 2, 2}, {6, 6, 6},          // deeper split planes
        {0, 0, 0}, {8, 8, 8},          // root box corners
        {0, 8, 0}, {8, 0, 8},          // edges of the root box
        {4, 0, 8}, {8, 4, 0},          // mixed boundary coords
    };
    for (const Point3d& p : pts) tree.insert(p);

    expect(tree.size() == pts.size(), "all boundary points inserted");

    // Whole-space query must return everything.
    expect(sameMultiset(tree.rangeQuery({0, 0, 0}, {8, 8, 8}), pts),
           "full-range query returns all");

    // Box whose faces lie exactly on split planes.
    expect(sameMultiset(tree.rangeQuery({4, 4, 4}, {8, 8, 8}),
                        bruteRange(pts, {4, 4, 4}, {8, 8, 8})),
           "range on split plane (upper octant)");
    expect(sameMultiset(tree.rangeQuery({0, 0, 0}, {4, 4, 4}),
                        bruteRange(pts, {0, 0, 0}, {4, 4, 4})),
           "range on split plane (lower octant)");

    // Degenerate box: a single coordinate.
    expect(sameMultiset(tree.rangeQuery({4, 4, 4}, {4, 4, 4}),
                        bruteRange(pts, {4, 4, 4}, {4, 4, 4})),
           "degenerate box query");

    // Radius reaching exactly a split-plane point (distance 4 from origin).
    expect(sameMultiset(tree.radiusSearch({0, 0, 0}, 4.0),
                        bruteRadius(pts, {0, 0, 0}, 4.0)),
           "radius touches split-plane point");

    // Nearest neighbor of a boundary point is itself.
    Point3d nearest;
    double dist = -1.0;
    expect(tree.nearestNeighbor({4, 4, 4}, nearest, &dist) && dist == 0.0 &&
               nearest.x == 4 && nearest.y == 4 && nearest.z == 4,
           "nearest of split-plane point is itself");
}

void testDuplicatePoints() {
    std::printf("[test] duplicate / coincident points\n");
    octree::Octree tree({-1, -1, -1}, {1, 1, 1}, 4, 10);
    std::vector<Point3d> pts;
    const int copies = 100;
    for (int i = 0; i < copies; ++i) {
        pts.push_back({0.5, 0.5, 0.5});   // all identical, on a split plane
        tree.insert(pts.back());
    }
    pts.push_back({-0.5, -0.5, -0.5});
    tree.insert(pts.back());

    expect(tree.size() == pts.size(), "duplicate insert count");

    auto inBox = tree.rangeQuery({0.5, 0.5, 0.5}, {0.5, 0.5, 0.5});
    expect(inBox.size() == static_cast<std::size_t>(copies),
           "degenerate range returns every duplicate");
    expect(sameMultiset(tree.rangeQuery({-1, -1, -1}, {1, 1, 1}), pts),
           "full range keeps all duplicates");

    auto ball = tree.radiusSearch({0.5, 0.5, 0.5}, 0.0);
    expect(ball.size() == static_cast<std::size_t>(copies),
           "zero-radius search returns every duplicate");

    Point3d nearest;
    double dist = -1.0;
    expect(tree.nearestNeighbor({0.5, 0.5, 0.5}, nearest, &dist) &&
               dist == 0.0,
           "nearest among duplicates has distance 0");
}

void testEmptyTree() {
    std::printf("[test] empty tree\n");
    octree::Octree tree({-10, -10, -10}, {10, 10, 10});
    expect(tree.empty(), "tree starts empty");
    expect(tree.rangeQuery({-10, -10, -10}, {10, 10, 10}).empty(),
           "range query on empty tree");
    expect(tree.radiusSearch({0, 0, 0}, 5).empty(),
           "radius search on empty tree");
    Point3d nearest;
    expect(!tree.nearestNeighbor({0, 0, 0}, nearest),
           "nearest neighbor on empty tree fails");
    expect(!tree.insert({100, 0, 0}), "insert outside root box rejected");
    expect(tree.empty(), "rejected insert leaves tree empty");
}

void testDepthLimit() {
    std::printf("[test] max depth limit\n");
    const std::size_t maxDepth = 6;
    octree::Octree tree({0, 0, 0}, {1, 1, 1}, 1, maxDepth);
    Rng rng(7);
    std::vector<Point3d> pts;
    // Cluster points tightly so subdivision wants to go very deep.
    for (int i = 0; i < 2000; ++i) {
        const double t = rng.uniform(0, 1e-3);
        pts.push_back({t, t * 0.5, t * 0.25});
        tree.insert(pts.back());
    }
    expect(tree.depth() <= maxDepth, "tree depth respects maxDepth");
    expect(tree.size() == pts.size(), "all clustered points stored");

    std::size_t ok = 0, total = 0;
    for (int i = 0; i < 50; ++i) {
        Point3d c{rng.uniform(0, 1e-3), rng.uniform(0, 5e-4),
                  rng.uniform(0, 2.5e-4)};
        const double r = rng.uniform(0, 1e-3);
        ++total;
        ok += sameMultiset(tree.radiusSearch(c, r), bruteRadius(pts, c, r));
    }
    report("clustered radius searches", ok, total);

    Point3d nearest;
    double dist = 0.0, refDist = 0.0;
    expect(tree.nearestNeighbor({1e-4, 1e-4, 1e-4}, nearest, &dist) &&
               bruteNearest(pts, {1e-4, 1e-4, 1e-4}, refDist) &&
               std::fabs(dist - refDist) < 1e-12,
           "nearest neighbor in deep cluster");
}

void testSaveLoad() {
    std::printf("[test] save / load consistency\n");
    const std::string path = "octree_test_dump.bin";
    Rng rng(123);
    std::vector<Point3d> pts;
    octree::Octree tree({-50, -50, -50}, {50, 50, 50}, 6, 12);
    for (int i = 0; i < 3000; ++i) {
        Point3d p{rng.uniform(-50, 50), rng.uniform(-50, 50),
                  rng.uniform(-50, 50)};
        pts.push_back(p);
        tree.insert(p);
    }
    for (int i = 0; i < 50; ++i) {  // duplicates on split planes
        pts.push_back({0, 0, 0});
        tree.insert(pts.back());
    }

    expect(tree.save(path), "save succeeds");

    octree::Octree loaded({0, 0, 0}, {1, 1, 1});
    expect(loaded.load(path), "load succeeds");
    expect(loaded.size() == tree.size(), "size preserved");
    expect(loaded.depth() == tree.depth(), "depth preserved");

    std::size_t ok = 0, total = 0;
    for (int i = 0; i < 100; ++i) {
        Point3d a{rng.uniform(-50, 50), rng.uniform(-50, 50),
                  rng.uniform(-50, 50)};
        Point3d b{rng.uniform(-50, 50), rng.uniform(-50, 50),
                  rng.uniform(-50, 50)};
        Point3d lo{std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z)};
        Point3d hi{std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z)};
        ++total;
        ok += sameMultiset(tree.rangeQuery(lo, hi), loaded.rangeQuery(lo, hi));
    }
    report("range queries identical after reload", ok, total);

    ok = 0; total = 0;
    for (int i = 0; i < 100; ++i) {
        Point3d c{rng.uniform(-50, 50), rng.uniform(-50, 50),
                  rng.uniform(-50, 50)};
        const double r = rng.uniform(0, 40);
        ++total;
        ok += sameMultiset(tree.radiusSearch(c, r), loaded.radiusSearch(c, r));
    }
    report("radius searches identical after reload", ok, total);

    ok = 0; total = 0;
    for (int i = 0; i < 200; ++i) {
        Point3d q{rng.uniform(-60, 60), rng.uniform(-60, 60),
                  rng.uniform(-60, 60)};
        Point3d n1, n2;
        double d1 = 0.0, d2 = 0.0;
        ++total;
        if (tree.nearestNeighbor(q, n1, &d1) &&
            loaded.nearestNeighbor(q, n2, &d2) && d1 == d2)
            ++ok;
    }
    report("nearest neighbors identical after reload", ok, total);

    std::remove(path.c_str());
}

}  // namespace

int main() {
    std::printf("==== Octree test suite (brute-force verified) ====\n");
    testRandomCloud();
    testBoundaryPoints();
    testDuplicatePoints();
    testEmptyTree();
    testDepthLimit();
    testSaveLoad();
    std::printf("==== summary: %d checks, %d failed, overall accuracy %.2f%% ====\n",
                g_checks, g_failures,
                g_checks == 0
                    ? 0.0
                    : 100.0 * (g_checks - g_failures) / g_checks);
    return g_failures == 0 ? 0 : 1;
}
