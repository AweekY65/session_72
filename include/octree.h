#pragma once

#include <array>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace octree {

struct Point3D {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    bool operator==(const Point3D& o) const {
        return x == o.x && y == o.y && z == o.z;
    }
};

struct AABB {
    Point3D min;
    Point3D max;

    bool contains(const Point3D& p) const;
    bool intersects(const AABB& o) const;
    bool intersectsSphere(const Point3D& c, double r) const;
    double minDistanceSq(const Point3D& p) const;
};

double distanceSq(const Point3D& a, const Point3D& b);

class Octree {
public:
    // If `bounds` is null the root bounds grow dynamically to contain
    // every inserted point (always kept cubic).
    explicit Octree(size_t capacity = 8, size_t maxDepth = 8,
                    const AABB* bounds = nullptr);

    void insert(const Point3D& p);

    // All points inside (inclusive) the axis-aligned box.
    std::vector<Point3D> rangeQuery(const AABB& box) const;

    // All points within Euclidean distance `radius` of `center`.
    std::vector<Point3D> radiusSearch(const Point3D& center, double radius) const;

    // Nearest neighbor by Euclidean distance; returns false on empty tree.
    bool nearestNeighbor(const Point3D& query, Point3D& result) const;

    bool save(const std::string& path) const;
    bool load(const std::string& path);

    size_t size() const { return size_; }
    size_t capacity() const { return capacity_; }
    size_t maxDepth() const { return maxDepth_; }
    AABB bounds() const;

private:
    struct Node {
        AABB box;
        bool isLeaf = true;
        std::vector<Point3D> points;
        std::array<std::unique_ptr<Node>, 8> children;

        explicit Node(const AABB& b) : box(b) {}
    };

    void insert(Node* node, const Point3D& p, size_t depth);
    int childIndex(const Node* node, const Point3D& p) const;
    AABB childBox(const Node* node, int index) const;
    void ensureContains(const Point3D& p);

    void rangeQuery(const Node* node, const AABB& box,
                    std::vector<Point3D>& out) const;
    void radiusSearch(const Node* node, const Point3D& c, double rSq,
                      std::vector<Point3D>& out) const;

    static void writeNode(std::ostream& os, const Node* node);
    static std::unique_ptr<Node> readNode(std::istream& is, bool& ok);

    size_t capacity_;
    size_t maxDepth_;
    std::unique_ptr<Node> root_;
    size_t size_ = 0;
};

}  // namespace octree
