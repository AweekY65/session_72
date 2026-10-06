#include "octree.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <functional>
#include <limits>
#include <queue>
#include <utility>

namespace octree {

namespace {
constexpr std::uint32_t kMagic = 0x4f435431;  // "OCT1"
constexpr std::uint8_t kLeafTag = 0;
constexpr std::uint8_t kInternalTag = 1;
}  // namespace

double squaredDistance(const Point3d& a, const Point3d& b) {
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    const double dz = a.z - b.z;
    return dx * dx + dy * dy + dz * dz;
}

bool contains(const Aabb& box, const Point3d& p) {
    return p.x >= box.min.x && p.x <= box.max.x && p.y >= box.min.y &&
           p.y <= box.max.y && p.z >= box.min.z && p.z <= box.max.z;
}

bool intersects(const Aabb& a, const Aabb& b) {
    return a.min.x <= b.max.x && a.max.x >= b.min.x && a.min.y <= b.max.y &&
           a.max.y >= b.min.y && a.min.z <= b.max.z && a.max.z >= b.min.z;
}

double squaredDistanceToBox(const Aabb& box, const Point3d& p) {
    double sum = 0.0;
    const double px[3] = {p.x, p.y, p.z};
    const double lo[3] = {box.min.x, box.min.y, box.min.z};
    const double hi[3] = {box.max.x, box.max.y, box.max.z};
    for (int axis = 0; axis < 3; ++axis) {
        double d = 0.0;
        if (px[axis] < lo[axis]) d = lo[axis] - px[axis];
        else if (px[axis] > hi[axis]) d = px[axis] - hi[axis];
        sum += d * d;
    }
    return sum;
}

Octree::Octree(const Point3d& minCorner, const Point3d& maxCorner,
               std::size_t capacity, std::size_t maxDepth)
    : capacity_(capacity == 0 ? 1 : capacity), maxDepth_(maxDepth) {
    root_ = makeNode(Aabb{minCorner, maxCorner});
}

std::unique_ptr<Octree::Node> Octree::makeNode(const Aabb& box) {
    auto node = std::make_unique<Node>();
    node->box = box;
    return node;
}

int Octree::childIndex(const Node& node, const Point3d& p) {
    // Points exactly on a split plane go to the upper (>= mid) child.
    const double midX = 0.5 * (node.box.min.x + node.box.max.x);
    const double midY = 0.5 * (node.box.min.y + node.box.max.y);
    const double midZ = 0.5 * (node.box.min.z + node.box.max.z);
    int index = 0;
    if (p.x >= midX) index |= 1;
    if (p.y >= midY) index |= 2;
    if (p.z >= midZ) index |= 4;
    return index;
}

Aabb Octree::childBox(const Node& node, int index) {
    Aabb box = node.box;
    const double midX = 0.5 * (box.min.x + box.max.x);
    const double midY = 0.5 * (box.min.y + box.max.y);
    const double midZ = 0.5 * (box.min.z + box.max.z);
    if (index & 1) box.min.x = midX; else box.max.x = midX;
    if (index & 2) box.min.y = midY; else box.max.y = midY;
    if (index & 4) box.min.z = midZ; else box.max.z = midZ;
    return box;
}

void Octree::subdivide(Node& node) {
    for (int i = 0; i < 8; ++i) node.children[i] = makeNode(childBox(node, i));
    std::vector<Point3d> pending;
    pending.swap(node.points);
    node.isLeaf = false;
    for (const Point3d& p : pending) {
        node.children[childIndex(node, p)]->points.push_back(p);
    }
}

void Octree::insertInto(Node& node, const Point3d& p, std::size_t depth) {
    if (!node.isLeaf) {
        insertInto(*node.children[childIndex(node, p)], p, depth + 1);
        return;
    }
    if (node.points.size() < capacity_ || depth >= maxDepth_) {
        // At maxDepth we always append, so coincident duplicates can never
        // trigger unbounded subdivision.
        node.points.push_back(p);
        return;
    }
    subdivide(node);
    insertInto(node, p, depth);
}

bool Octree::insert(const Point3d& p) {
    if (!contains(root_->box, p)) return false;
    insertInto(*root_, p, 0);
    ++size_;
    return true;
}

std::vector<Point3d> Octree::rangeQuery(const Point3d& minCorner,
                                        const Point3d& maxCorner) const {
    std::vector<Point3d> out;
    rangeQueryNode(*root_, Aabb{minCorner, maxCorner}, out);
    return out;
}

void Octree::rangeQueryNode(const Node& node, const Aabb& query,
                            std::vector<Point3d>& out) {
    if (!intersects(node.box, query)) return;
    if (node.isLeaf) {
        for (const Point3d& p : node.points) {
            if (contains(query, p)) out.push_back(p);
        }
        return;
    }
    for (const auto& child : node.children) rangeQueryNode(*child, query, out);
}

std::vector<Point3d> Octree::radiusSearch(const Point3d& center,
                                          double radius) const {
    std::vector<Point3d> out;
    if (radius < 0.0) return out;
    radiusSearchNode(*root_, center, radius * radius, out);
    return out;
}

void Octree::radiusSearchNode(const Node& node, const Point3d& center,
                              double squaredRadius, std::vector<Point3d>& out) {
    if (squaredDistanceToBox(node.box, center) > squaredRadius) return;
    if (node.isLeaf) {
        for (const Point3d& p : node.points) {
            if (squaredDistance(center, p) <= squaredRadius) out.push_back(p);
        }
        return;
    }
    for (const auto& child : node.children) {
        radiusSearchNode(*child, center, squaredRadius, out);
    }
}

bool Octree::nearestNeighbor(const Point3d& query, Point3d& nearest,
                             double* distance) const {
    if (size_ == 0) return false;

    // Best-first search: nodes are expanded in order of their bounding-box
    // lower-bound distance, so any node whose bound exceeds the current best
    // (and everything still queued behind it) is pruned.
    using Entry = std::pair<double, const Node*>;
    std::priority_queue<Entry, std::vector<Entry>, std::greater<Entry>> queue;
    queue.emplace(squaredDistanceToBox(root_->box, query), root_.get());

    double best = std::numeric_limits<double>::infinity();
    bool found = false;
    while (!queue.empty()) {
        const Entry top = queue.top();
        queue.pop();
        if (top.first > best) break;  // remaining nodes cannot improve
        const Node& node = *top.second;
        if (node.isLeaf) {
            for (const Point3d& p : node.points) {
                const double d = squaredDistance(query, p);
                if (d < best) {
                    best = d;
                    nearest = p;
                    found = true;
                }
            }
        } else {
            for (const auto& child : node.children) {
                const double bound = squaredDistanceToBox(child->box, query);
                if (bound <= best) queue.emplace(bound, child.get());
            }
        }
    }
    if (found && distance) *distance = std::sqrt(best);
    return found;
}

std::size_t Octree::depth() const { return depthOf(*root_); }

std::size_t Octree::depthOf(const Node& node) {
    if (node.isLeaf) return 0;
    std::size_t deepest = 0;
    for (const auto& child : node.children) {
        deepest = std::max(deepest, depthOf(*child));
    }
    return deepest + 1;
}

bool Octree::save(const std::string& path) const {
    std::ofstream os(path, std::ios::binary);
    if (!os) return false;
    os.write(reinterpret_cast<const char*>(&kMagic), sizeof(kMagic));
    const std::uint64_t capacity = capacity_;
    const std::uint64_t maxDepth = maxDepth_;
    const std::uint64_t count = size_;
    os.write(reinterpret_cast<const char*>(&capacity), sizeof(capacity));
    os.write(reinterpret_cast<const char*>(&maxDepth), sizeof(maxDepth));
    os.write(reinterpret_cast<const char*>(&count), sizeof(count));
    os.write(reinterpret_cast<const char*>(&root_->box), sizeof(Aabb));
    return saveNode(*root_, os) && os.good();
}

bool Octree::saveNode(const Node& node, std::ostream& os) {
    if (node.isLeaf) {
        os.write(reinterpret_cast<const char*>(&kLeafTag), sizeof(kLeafTag));
        const std::uint64_t count = node.points.size();
        os.write(reinterpret_cast<const char*>(&count), sizeof(count));
        if (count > 0) {
            os.write(reinterpret_cast<const char*>(node.points.data()),
                     static_cast<std::streamsize>(count * sizeof(Point3d)));
        }
    } else {
        os.write(reinterpret_cast<const char*>(&kInternalTag),
                 sizeof(kInternalTag));
        for (const auto& child : node.children) {
            if (!saveNode(*child, os)) return false;
        }
    }
    return os.good();
}

bool Octree::load(const std::string& path) {
    std::ifstream is(path, std::ios::binary);
    if (!is) return false;
    std::uint32_t magic = 0;
    is.read(reinterpret_cast<char*>(&magic), sizeof(magic));
    if (!is || magic != kMagic) return false;
    std::uint64_t capacity = 0, maxDepth = 0, count = 0;
    is.read(reinterpret_cast<char*>(&capacity), sizeof(capacity));
    is.read(reinterpret_cast<char*>(&maxDepth), sizeof(maxDepth));
    is.read(reinterpret_cast<char*>(&count), sizeof(count));
    Aabb box{};
    is.read(reinterpret_cast<char*>(&box), sizeof(box));
    if (!is || capacity == 0) return false;

    auto root = makeNode(box);
    if (!loadNode(*root, is)) return false;

    root_ = std::move(root);
    capacity_ = static_cast<std::size_t>(capacity);
    maxDepth_ = static_cast<std::size_t>(maxDepth);
    size_ = static_cast<std::size_t>(count);
    return true;
}

bool Octree::loadNode(Node& node, std::istream& is) {
    std::uint8_t tag = 0;
    is.read(reinterpret_cast<char*>(&tag), sizeof(tag));
    if (!is) return false;
    if (tag == kLeafTag) {
        node.isLeaf = true;
        std::uint64_t count = 0;
        is.read(reinterpret_cast<char*>(&count), sizeof(count));
        if (!is) return false;
        node.points.resize(static_cast<std::size_t>(count));
        if (count > 0) {
            is.read(reinterpret_cast<char*>(node.points.data()),
                    static_cast<std::streamsize>(count * sizeof(Point3d)));
        }
        return is.good();
    }
    if (tag == kInternalTag) {
        node.isLeaf = false;
        for (auto& child : node.children) {
            child = makeNode(Aabb{});
        }
        // Child boxes are derived from the parent box instead of trusting the
        // file for geometry (only the root box is persisted). They must be
        // assigned before recursing so deeper levels derive correct boxes.
        for (int i = 0; i < 8; ++i) {
            node.children[i]->box = childBox(node, i);
            if (!loadNode(*node.children[i], is)) return false;
        }
        return true;
    }
    return false;
}

}  // namespace octree
