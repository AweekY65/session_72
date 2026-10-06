# Local 3D Octree Spatial Index

一个纯本地的三维 Octree 空间索引实现（C++17）。所有点云、树节点、查询状态与序列化数据只存在于**内存或本地文件**中，不依赖 PCL、数据库、ROS 或任何外部服务。

## 功能

- **点插入**：`insert(p)`，叶子节点超过 `capacity` 且未达 `maxDepth` 时递归划分为 8 个子空间。
- **AABB 范围查询**：`rangeQuery(lo, hi)`，闭区间，返回范围内全部点（含重复点），不重不漏。
- **半径搜索**：`radiusSearch(center, r)`，返回欧氏距离 `<= r` 的全部点。
- **最近邻查询**：`nearestNeighbor(q, out, &dist)`，基于节点 bounding box 的 best-first 剪枝搜索（非全点遍历）。
- **序列化**：`save(path)` / `load(path)`，二进制本地文件保存完整树结构，加载后查询结果与原树完全一致。

## 空间划分规则

- 每个节点持有轴对齐包围盒（AABB）。内部节点在三个轴的中点处切分，产生 8 个全等子空间。
- 恰好落在划分平面上的点按 `>= mid` 规则归入**上侧**子节点，规则确定且无歧义，保证插入与查询走同一条路径。
- 根节点包围盒在构造时指定（闭区间），落在根盒外的点插入返回 `false`。
- 叶子达到 `capacity` 后划分；到达 `maxDepth` 的叶子不再划分、直接追加。因此完全相同位置的重复点不会导致无限递归，任意数量重复点都可安全存储。

## 剪枝规则

- **范围查询**：节点盒与查询盒不相交（`intersects` 为假）时整棵子树跳过；叶子内逐点做闭区间包含测试。
- **半径搜索**：查询点到节点盒的最短平方距离 `> r²` 时整棵子树跳过；叶子内逐点比较平方距离，避免开方。
- **最近邻**：优先队列按"查询点到节点盒的最短距离"做 best-first 展开；当队首下界超过当前最优距离时终止，其余节点全部剪枝。只访问可能含更优解的节点，绝不简单遍历全部点。

## 复杂度

设 `n` 为点数、`h` 为树高（`h <= maxDepth`）、`k` 为结果数：

| 操作 | 平均 | 最坏 |
|---|---|---|
| 插入 | O(h) | O(h) |
| 范围查询 | O(h + k) | O(n) |
| 半径搜索 | O(h + k) | O(n) |
| 最近邻 | O(h) 期望 | O(n)（点分布退化时） |
| 空间 | O(n) | O(n) |

## 构建与测试

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/test_octree      # 直接运行，逐项打印 PASS/FAIL 与正确率
# 或
cd build && ctest --output-on-failure
```

测试覆盖：随机点云（5000 点，200 次范围查询 / 200 次半径搜索 / 500 次最近邻）、划分平面与空间边界点、完全相同位置的重复点、空树、最大深度限制、保存/加载一致性，全部与 brute-force 结果做多重集合（multiset）对照，终端输出每项正确率与总体正确率（当前为 100%）。

## 序列化格式

二进制本地文件：`magic("OCT1") | capacity | maxDepth | pointCount | rootBox | 递归节点流`。节点流中叶子记录点数与点数组，内部节点递归记录 8 个子节点；子节点包围盒加载时由父盒推导，不依赖文件中的几何数据。

## 接口速览

```cpp
#include "octree.h"

octree::Octree tree({-100, -100, -100}, {100, 100, 100}, /*capacity=*/8, /*maxDepth=*/16);
tree.insert({1, 2, 3});
auto inBox  = tree.rangeQuery({0, 0, 0}, {5, 5, 5});
auto inBall = tree.radiusSearch({0, 0, 0}, 3.0);
octree::Point3d nn; double dist;
bool ok = tree.nearestNeighbor({1, 1, 1}, nn, &dist);
tree.save("tree.bin");
octree::Octree copy({0,0,0},{1,1,1});
copy.load("tree.bin");
```

## 目录结构

```
include/octree.h      接口与几何工具
src/octree.cpp        插入/划分/查询/序列化实现
tests/test_octree.cpp 自动化测试（brute-force 对照 + 正确率输出）
CMakeLists.txt        构建脚本
```
