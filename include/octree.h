#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <memory>
#include <string>
#include <vector>

namespace octree {

struct Point3d {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

struct Aabb {
    Point3d min;
    Point3d max;
};

double squaredDistance(const Point3d& a, const Point3d& b);
bool contains(const Aabb& box, const Point3d& p);
bool intersects(const Aabb& a, const Aabb& b);
double squaredDistanceToBox(const Aabb& box, const Point3d& p);

class Octree {
public:
    // minCorner/maxCorner define the (inclusive) root bounding box.
    // capacity: max points stored in a leaf before it subdivides.
    // maxDepth: maximum tree depth; leaves at maxDepth keep all points
    // (this is what makes duplicate / coincident points safe).
    Octree(const Point3d& minCorner, const Point3d& maxCorner,
           std::size_t capacity = 8, std::size_t maxDepth = 16);

    // Returns false when the point lies outside the root bounding box.
    bool insert(const Point3d& p);

    // Inclusive axis-aligned box range query. Returns every stored point
    // with min <= p <= max component-wise (duplicates included).
    std::vector<Point3d> rangeQuery(const Point3d& minCorner,
                                    const Point3d& maxCorner) const;

    // Euclidean ball query: all points with distance(center, p) <= radius.
    std::vector<Point3d> radiusSearch(const Point3d& center, double radius) const;

    // Best-first nearest neighbor using node bounding boxes for pruning.
    // Returns false on an empty tree. On ties any nearest point may be returned.
    bool nearestNeighbor(const Point3d& query, Point3d& nearest,
                         double* distance = nullptr) const;

    // Binary serialization to / from a local file. The full node structure
    // is stored, so a loaded tree answers queries identically.
    bool save(const std::string& path) const;
    bool load(const std::string& path);

    std::size_t size() const { return size_; }
    bool empty() const { return size_ == 0; }
    std::size_t depth() const;
    const Aabb& bounds() const { return root_->box; }

private:
    struct Node {
        Aabb box;
        bool isLeaf = true;
        std::vector<Point3d> points;
        std::array<std::unique_ptr<Node>, 8> children;
    };

    std::unique_ptr<Node> root_;
    std::size_t capacity_;
    std::size_t maxDepth_;
    std::size_t size_ = 0;

    static std::unique_ptr<Node> makeNode(const Aabb& box);
    static int childIndex(const Node& node, const Point3d& p);
    static Aabb childBox(const Node& node, int index);
    void subdivide(Node& node);
    void insertInto(Node& node, const Point3d& p, std::size_t depth);
    static void rangeQueryNode(const Node& node, const Aabb& query,
                               std::vector<Point3d>& out);
    static void radiusSearchNode(const Node& node, const Point3d& center,
                                 double squaredRadius, std::vector<Point3d>& out);
    static std::size_t depthOf(const Node& node);
    static bool saveNode(const Node& node, std::ostream& os);
    static bool loadNode(Node& node, std::istream& is);
};

}  // namespace octree
