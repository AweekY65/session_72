# Local 3D Octree Spatial Index

一个纯本地的三维 Octree 空间索引（C++17，零第三方依赖）。所有点云、
树节点、查询状态与序列化数据只存在于内存或本地文件中，不依赖 PCL
服务、数据库、ROS 或任何外部服务。

## 功能

- **点插入**：`insert(Point3D)`，按节点容量 `capacity` 与最大深度
  `maxDepth` 递归划分为 8 个子空间（octant）。
- **范围查询**：`rangeQuery(AABB)` 返回轴对齐盒内全部点，不重不漏。
- **半径搜索**：`radiusSearch(center, r)` 返回欧氏距离 ≤ r 的全部点。
- **最近邻**：`nearestNeighbor(query)` 基于节点 bounding box 的
  best-first 剪枝搜索，绝不退化为全点遍历。
- **序列化**：`save(path)` / `load(path)` 以二进制格式保存/加载整棵树，
  加载后查询结果与原树完全一致。

## 空间划分规则

- 每个节点持有一个 axis-aligned bounding box（AABB）。叶节点点数超过
  `capacity` 且深度小于 `maxDepth` 时，沿 x/y/z 三个中位面划分为 8 个
  全等子盒，点被重新分配到对应子节点。
- **划分平面归属**：位于中位面上的点统一归入上侧（`>=` 中点）子节点，
  规则确定，因此查询时按闭区间相交判断即可保证不重不漏。
- **重复点**：完全相同坐标的点全部保留（多重集合语义），它们会沿同一
  路径下沉；到达 `maxDepth` 或子盒尺寸低于 `1e-9` 时停止划分，直接存于
  当前叶节点。
- **根边界**：构造时可指定固定边界；否则根盒从第一个点出发动态立方
  扩张（每次边长翻倍），保证任意位置的点都能插入。

## 剪枝规则

- **Range query**：仅递归进入与查询盒相交（闭区间 `intersects`）的
  节点；叶节点内逐点用 `contains` 精确过滤。
- **Radius search**：计算查询点到节点 AABB 的最小距离平方
  `minDistanceSq`，若大于 r² 则整棵子树被剪掉。
- **Nearest neighbor**：使用最小堆按"节点 AABB 到查询点的最小距离"
  best-first 展开；当堆顶距离超过当前最优点的距离时终止，堆中其余
  节点全部被剪枝。点与 AABB 距离均为下界，因此结果精确。

## 复杂度

设 n 为点数，h 为树高（h ≤ maxDepth），k 为结果数：

| 操作 | 平均 | 最坏 |
|---|---|---|
| 插入 | O(h) | O(maxDepth) |
| 范围/半径查询 | O(h + k) | O(n)（点极度聚集时） |
| 最近邻 | O(h) 期望 | O(n)（退化分布时） |
| 序列化 | O(n) | O(n) |

空间占用 O(n)。相比 brute-force 的 O(n) 最近邻/查询，剪枝后只访问与
查询区域相交的少量节点。

## 构建与测试

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
./build/test_octree        # 直接运行，终端输出每组测试正确率
ctest --test-dir build --output-on-failure   # 或通过 CTest 运行
```

测试覆盖（均与 brute-force 结果逐项对照，输出 accuracy）：

- 5000 个随机点的 range / radius / nearest 各 200 次随机查询
- 位于划分平面、盒边界、角点上的点，零尺寸查询盒、零半径搜索
- 50 个完全相同坐标的重复点
- 空树的全部查询行为
- `capacity=1, maxDepth=4` 的强制浅树（深度限制）下结果正确性
- 保存→重新加载后 range / radius / nearest 结果与原树一致
- 远离原点（±1e6）的点触发根盒动态扩张

## 序列化格式

二进制小端布局：8 字节 magic `OCTREE01`，随后是 capacity、maxDepth、
点总数（uint64），是否含根节点（uint8），然后递归写入节点：节点 AABB
（6 个 double）+ 是否叶节点（uint8）+ 叶节点点数与坐标，或 8 个子节点
存在标志及子节点数据。加载时校验 magic 与读取完整性，非法文件返回
`false` 且不破坏现有树。

## 文件结构

```
include/octree.h    接口：Point3D / AABB / Octree
src/octree.cpp      实现：插入、划分、三类查询、序列化
tests/test_octree.cpp  自动化测试（brute-force 对照 + 正确率输出）
CMakeLists.txt
```
