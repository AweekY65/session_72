#include "octree.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <limits>
#include <queue>

namespace octree {

namespace {
constexpr char kMagic[8] = {'O', 'C', 'T', 'R', 'E', 'E', '0', '1'};
constexpr double kMinExtent = 1e-9;
}

bool AABB::contains(const Point3D& p) const {
    return p.x >= min.x && p.x <= max.x && p.y >= min.y && p.y <= max.y &&
           p.z >= min.z && p.z <= max.z;
}

bool AABB::intersects(const AABB& o) const {
    return min.x <= o.max.x && max.x >= o.min.x && min.y <= o.max.y &&
           max.y >= o.min.y && min.z <= o.max.z && max.z >= o.min.z;
}

bool AABB::intersectsSphere(const Point3D& c, double r) const {
    return minDistanceSq(c) <= r * r;
}

double AABB::minDistanceSq(const Point3D& p) const {
    double d = 0.0;
    const double pv[3] = {p.x, p.y, p.z};
    const double lo[3] = {min.x, min.y, min.z};
    const double hi[3] = {max.x, max.y, max.z};
    for (int i = 0; i < 3; ++i) {
        if (pv[i] < lo[i]) {
            d += (lo[i] - pv[i]) * (lo[i] - pv[i]);
        } else if (pv[i] > hi[i]) {
            d += (pv[i] - hi[i]) * (pv[i] - hi[i]);
        }
    }
    return d;
}

double distanceSq(const Point3D& a, const Point3D& b) {
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    const double dz = a.z - b.z;
    return dx * dx + dy * dy + dz * dz;
}

Octree::Octree(size_t capacity, size_t maxDepth, const AABB* bounds)
    : capacity_(capacity == 0 ? 1 : capacity), maxDepth_(maxDepth) {
    if (bounds) {
        root_ = std::make_unique<Node>(*bounds);
    }
}

AABB Octree::bounds() const {
    if (root_) return root_->box;
    return AABB{};
}

int Octree::childIndex(const Node* node, const Point3D& p) const {
    const Point3D& lo = node->box.min;
    const Point3D& hi = node->box.max;
    const double mx = (lo.x + hi.x) * 0.5;
    const double my = (lo.y + hi.y) * 0.5;
    const double mz = (lo.z + hi.z) * 0.5;
    // Points exactly on a split plane go to the upper (>=) child.
    int idx = 0;
    if (p.x >= mx) idx |= 1;
    if (p.y >= my) idx |= 2;
    if (p.z >= mz) idx |= 4;
    return idx;
}

AABB Octree::childBox(const Node* node, int index) const {
    const Point3D& lo = node->box.min;
    const Point3D& hi = node->box.max;
    const double mx = (lo.x + hi.x) * 0.5;
    const double my = (lo.y + hi.y) * 0.5;
    const double mz = (lo.z + hi.z) * 0.5;
    AABB b;
    b.min.x = (index & 1) ? mx : lo.x;
    b.max.x = (index & 1) ? hi.x : mx;
    b.min.y = (index & 2) ? my : lo.y;
    b.max.y = (index & 2) ? hi.y : my;
    b.min.z = (index & 4) ? mz : lo.z;
    b.max.z = (index & 4) ? hi.z : mz;
    return b;
}

void Octree::ensureContains(const Point3D& p) {
    if (!root_) {
        AABB b{{p.x - 0.5, p.y - 0.5, p.z - 0.5},
               {p.x + 0.5, p.y + 0.5, p.z + 0.5}};
        root_ = std::make_unique<Node>(b);
        return;
    }
    if (root_->box.contains(p)) return;

    // Grow the root box (kept cubic) until it contains p.
    AABB b = root_->box;
    double cx = (b.min.x + b.max.x) * 0.5;
    double cy = (b.min.y + b.max.y) * 0.5;
    double cz = (b.min.z + b.max.z) * 0.5;
    double half = std::max({b.max.x - b.min.x, b.max.y - b.min.y,
                            b.max.z - b.min.z, kMinExtent}) * 0.5;
    const double need = std::max({std::fabs(p.x - cx), std::fabs(p.y - cy),
                                  std::fabs(p.z - cz)});
    while (half < need) half *= 2.0;
    b.min = {cx - half, cy - half, cz - half};
    b.max = {cx + half, cy + half, cz + half};
    root_->box = b;
}

void Octree::insert(const Point3D& p) {
    ensureContains(p);
    insert(root_.get(), p, 0);
    ++size_;
}

void Octree::insert(Node* node, const Point3D& p, size_t depth) {
    if (!node->isLeaf) {
        insert(node->children[childIndex(node, p)].get(), p, depth + 1);
        return;
    }
    if (node->points.size() < capacity_ || depth >= maxDepth_) {
        node->points.push_back(p);
        return;
    }
    // Leaf is full and we may still split: subdivide, then relocate.
    const double ex = node->box.max.x - node->box.min.x;
    const double ey = node->box.max.y - node->box.min.y;
    const double ez = node->box.max.z - node->box.min.z;
    if (std::min({ex, ey, ez}) <= kMinExtent) {
        // Box too small to split meaningfully; just store.
        node->points.push_back(p);
        return;
    }
    std::vector<Point3D> relocated;
    relocated.swap(node->points);
    for (int i = 0; i < 8; ++i) {
        node->children[i] = std::make_unique<Node>(childBox(node, i));
    }
    node->isLeaf = false;
    for (const Point3D& q : relocated) {
        insert(node->children[childIndex(node, q)].get(), q, depth + 1);
    }
    insert(node->children[childIndex(node, p)].get(), p, depth + 1);
}

std::vector<Point3D> Octree::rangeQuery(const AABB& box) const {
    std::vector<Point3D> out;
    if (root_) rangeQuery(root_.get(), box, out);
    return out;
}

void Octree::rangeQuery(const Node* node, const AABB& box,
                        std::vector<Point3D>& out) const {
    if (!node->box.intersects(box)) return;
    if (node->isLeaf) {
        for (const Point3D& p : node->points) {
            if (box.contains(p)) out.push_back(p);
        }
        return;
    }
    for (const auto& child : node->children) {
        if (child) rangeQuery(child.get(), box, out);
    }
}

std::vector<Point3D> Octree::radiusSearch(const Point3D& center,
                                          double radius) const {
    std::vector<Point3D> out;
    if (root_ && radius >= 0.0) {
        radiusSearch(root_.get(), center, radius * radius, out);
    }
    return out;
}

void Octree::radiusSearch(const Node* node, const Point3D& c, double rSq,
                          std::vector<Point3D>& out) const {
    if (node->box.minDistanceSq(c) > rSq) return;
    if (node->isLeaf) {
        for (const Point3D& p : node->points) {
            if (distanceSq(p, c) <= rSq) out.push_back(p);
        }
        return;
    }
    for (const auto& child : node->children) {
        if (child) radiusSearch(child.get(), c, rSq, out);
    }
}

bool Octree::nearestNeighbor(const Point3D& query, Point3D& result) const {
    if (!root_ || size_ == 0) return false;

    // Best-first search: expand the node whose bounding box is closest to
    // the query. Any node whose box is farther than the current best point
    // is pruned, so most of the tree is never visited.
    struct Entry {
        double distSq;
        const Node* node;
        bool operator>(const Entry& o) const { return distSq > o.distSq; }
    };
    std::priority_queue<Entry, std::vector<Entry>, std::greater<Entry>> pq;
    pq.push({root_->box.minDistanceSq(query), root_.get()});

    double bestSq = std::numeric_limits<double>::infinity();
    const Point3D* best = nullptr;

    while (!pq.empty()) {
        Entry e = pq.top();
        pq.pop();
        if (e.distSq > bestSq) break;  // everything left is worse
        const Node* node = e.node;
        if (node->isLeaf) {
            for (const Point3D& p : node->points) {
                const double d = distanceSq(p, query);
                if (d < bestSq) {
                    bestSq = d;
                    best = &p;
                }
            }
        } else {
            for (const auto& child : node->children) {
                if (!child) continue;
                const double d = child->box.minDistanceSq(query);
                if (d <= bestSq) pq.push({d, child.get()});
            }
        }
    }
    if (!best) return false;
    result = *best;
    return true;
}

bool Octree::save(const std::string& path) const {
    std::ofstream os(path, std::ios::binary);
    if (!os) return false;
    os.write(kMagic, sizeof(kMagic));
    const uint64_t cap = capacity_, depth = maxDepth_, n = size_;
    os.write(reinterpret_cast<const char*>(&cap), sizeof(cap));
    os.write(reinterpret_cast<const char*>(&depth), sizeof(depth));
    os.write(reinterpret_cast<const char*>(&n), sizeof(n));
    const uint8_t hasRoot = root_ ? 1 : 0;
    os.write(reinterpret_cast<const char*>(&hasRoot), sizeof(hasRoot));
    if (root_) writeNode(os, root_.get());
    return os.good();
}

void Octree::writeNode(std::ostream& os, const Node* node) {
    const double box[6] = {node->box.min.x, node->box.min.y, node->box.min.z,
                           node->box.max.x, node->box.max.y, node->box.max.z};
    os.write(reinterpret_cast<const char*>(box), sizeof(box));
    const uint8_t leaf = node->isLeaf ? 1 : 0;
    os.write(reinterpret_cast<const char*>(&leaf), sizeof(leaf));
    if (node->isLeaf) {
        const uint64_t count = node->points.size();
        os.write(reinterpret_cast<const char*>(&count), sizeof(count));
        for (const Point3D& p : node->points) {
            os.write(reinterpret_cast<const char*>(&p.x), sizeof(p.x));
            os.write(reinterpret_cast<const char*>(&p.y), sizeof(p.y));
            os.write(reinterpret_cast<const char*>(&p.z), sizeof(p.z));
        }
    } else {
        for (const auto& child : node->children) {
            const uint8_t present = child ? 1 : 0;
            os.write(reinterpret_cast<const char*>(&present), sizeof(present));
            if (child) writeNode(os, child.get());
        }
    }
}

bool Octree::load(const std::string& path) {
    std::ifstream is(path, std::ios::binary);
    if (!is) return false;
    char magic[8];
    is.read(magic, sizeof(magic));
    if (!is || std::memcmp(magic, kMagic, sizeof(kMagic)) != 0) return false;
    uint64_t cap = 0, depth = 0, n = 0;
    is.read(reinterpret_cast<char*>(&cap), sizeof(cap));
    is.read(reinterpret_cast<char*>(&depth), sizeof(depth));
    is.read(reinterpret_cast<char*>(&n), sizeof(n));
    uint8_t hasRoot = 0;
    is.read(reinterpret_cast<char*>(&hasRoot), sizeof(hasRoot));
    if (!is) return false;

    std::unique_ptr<Node> newRoot;
    if (hasRoot) {
        bool ok = true;
        newRoot = readNode(is, ok);
        if (!ok || !newRoot) return false;
    }
    capacity_ = cap == 0 ? 1 : static_cast<size_t>(cap);
    maxDepth_ = static_cast<size_t>(depth);
    size_ = static_cast<size_t>(n);
    root_ = std::move(newRoot);
    return true;
}

std::unique_ptr<Octree::Node> Octree::readNode(std::istream& is, bool& ok) {
    double box[6];
    is.read(reinterpret_cast<char*>(box), sizeof(box));
    uint8_t leaf = 0;
    is.read(reinterpret_cast<char*>(&leaf), sizeof(leaf));
    if (!is) {
        ok = false;
        return nullptr;
    }
    AABB b{{box[0], box[1], box[2]}, {box[3], box[4], box[5]}};
    auto node = std::make_unique<Node>(b);
    node->isLeaf = leaf != 0;
    if (node->isLeaf) {
        uint64_t count = 0;
        is.read(reinterpret_cast<char*>(&count), sizeof(count));
        if (!is || count > (1ull << 32)) {
            ok = false;
            return nullptr;
        }
        node->points.reserve(static_cast<size_t>(count));
        for (uint64_t i = 0; i < count; ++i) {
            Point3D p;
            is.read(reinterpret_cast<char*>(&p.x), sizeof(p.x));
            is.read(reinterpret_cast<char*>(&p.y), sizeof(p.y));
            is.read(reinterpret_cast<char*>(&p.z), sizeof(p.z));
            if (!is) {
                ok = false;
                return nullptr;
            }
            node->points.push_back(p);
        }
    } else {
        for (int i = 0; i < 8; ++i) {
            uint8_t present = 0;
            is.read(reinterpret_cast<char*>(&present), sizeof(present));
            if (!is) {
                ok = false;
                return nullptr;
            }
            if (present) {
                node->children[i] = readNode(is, ok);
                if (!ok) return nullptr;
            }
        }
    }
    return node;
}

}  // namespace octree
