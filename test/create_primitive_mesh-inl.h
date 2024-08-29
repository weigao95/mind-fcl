#pragma once

#include <cstdlib>

template <typename BV>
::fcl::BVHModel<BV>* fcl::test::makeCubeBVH(typename BV::S size_x,
                                            typename BV::S size_y,
                                            typename BV::S size_z) {
  using Scalar = typename BV::S;
  int faces[6][4] = {{0, 1, 2, 3}, {3, 2, 6, 7}, {7, 6, 5, 4},
                     {4, 5, 1, 0}, {5, 6, 2, 1}, {7, 4, 0, 3}};
  Scalar v[8][3];

  v[0][0] = v[1][0] = v[2][0] = v[3][0] = -size_x / 2;
  v[4][0] = v[5][0] = v[6][0] = v[7][0] = size_x / 2;
  v[0][1] = v[1][1] = v[4][1] = v[5][1] = -size_y / 2;
  v[2][1] = v[3][1] = v[6][1] = v[7][1] = size_y / 2;
  v[0][2] = v[3][2] = v[4][2] = v[7][2] = -size_z / 2;
  v[1][2] = v[2][2] = v[5][2] = v[6][2] = size_z / 2;

  ::fcl::BVHModel<BV>* model = new ::fcl::BVHModel<BV>;
  fcl::Vector3<Scalar> p1, p2, p3;
  model->beginModel();

  for (int i = 0; i < 6; i++) {
    p1 = fcl::Vector3<Scalar>(v[faces[i][0]][0], v[faces[i][0]][1],
                              v[faces[i][0]][2]);
    p2 = fcl::Vector3<Scalar>(v[faces[i][1]][0], v[faces[i][1]][1],
                              v[faces[i][1]][2]);
    p3 = fcl::Vector3<Scalar>(v[faces[i][2]][0], v[faces[i][2]][1],
                              v[faces[i][2]][2]);
    model->addTriangle(p1, p2, p3);

    p1 = fcl::Vector3<Scalar>(v[faces[i][0]][0], v[faces[i][0]][1],
                              v[faces[i][0]][2]);
    p2 = fcl::Vector3<Scalar>(v[faces[i][2]][0], v[faces[i][2]][1],
                              v[faces[i][2]][2]);
    p3 = fcl::Vector3<Scalar>(v[faces[i][3]][0], v[faces[i][3]][1],
                              v[faces[i][3]][2]);
    model->addTriangle(p1, p2, p3);
  }
  model->endModel();
  return model;
}

template <typename BV>
::fcl::BVHModel<BV>* fcl::test::makeEllipsoidBVH(typename BV::S size_x,
                                                 typename BV::S size_y,
                                                 typename BV::S size_z) {
  using Scalar = typename BV::S;

  // Make the triangles
  std::vector<std::array<Scalar, 3>> points;
  std::vector<std::array<int, 3>> triangles;
  meshEllipsoidTriangles(size_x, size_y, size_z, &points, &triangles);

  // Then make bvh
  fcl::BVHModel<BV>* model = new ::fcl::BVHModel<BV>;
  fcl::Vector3<Scalar> p1, p2, p3;
  model->beginModel();

  for (std::size_t i = 0; i < triangles.size(); i++) {
    int p1_idx = triangles[i][0];
    int p2_idx = triangles[i][1];
    int p3_idx = triangles[i][2];

    // Make points
    // clang-format off
    p1 = fcl::Vector3<Scalar>(points[p1_idx][0], points[p1_idx][1], points[p1_idx][2]);
    p2 = fcl::Vector3<Scalar>(points[p2_idx][0], points[p2_idx][1], points[p2_idx][2]);
    p3 = fcl::Vector3<Scalar>(points[p3_idx][0], points[p3_idx][1], points[p3_idx][2]);
    // clang-format on

    model->addTriangle(p1, p2, p3);
  }

  model->endModel();
  return model;
}

template <typename T>
::fcl::Convex<T> fcl::test::makeEllipsoidConvex(T size_x, T size_y, T size_z) {
  // Make the triangles
  std::vector<std::array<T, 3>> points;
  std::vector<std::array<int, 3>> triangles;
  meshEllipsoidTriangles(size_x, size_y, size_z, &points, &triangles);

  // The format used by fcl::Convex
  std::shared_ptr<std::vector<fcl::Vector3<T>>> vertices =
      std::make_shared<std::vector<fcl::Vector3<T>>>();
  int num_faces = triangles.size();
  std::shared_ptr<std::vector<int>> faces =
      std::make_shared<std::vector<int>>();

  // The vertex is simple
  for (std::size_t i = 0; i < points.size(); i++) {
    vertices->push_back(
        fcl::Vector3<T>(points[i][0], points[i][1], points[i][2]));
  }

  // The face is more complex
  for (std::size_t i = 0; i < triangles.size(); i++) {
    faces->push_back(3);
    faces->push_back(triangles[i][0]);
    faces->push_back(triangles[i][1]);
    faces->push_back(triangles[i][2]);
  }

  fcl::Convex<T> convex(vertices, num_faces, faces, true);
  return convex;
}

template <typename T>
::fcl::Convex<T> fcl::test::makeEllipsoidConvexUV(T size_x, T size_y, T size_z,
                                                  int n_slices, int n_stacks) {
  // Make the triangles
  std::vector<std::array<T, 3>> points;
  std::vector<std::array<int, 3>> triangles;
  meshEllipsoidUV(size_x, size_y, size_z, n_slices, n_stacks, &points,
                  &triangles);

  // The format used by fcl::Convex
  std::shared_ptr<std::vector<fcl::Vector3<T>>> vertices =
      std::make_shared<std::vector<fcl::Vector3<T>>>();
  int num_faces = triangles.size();
  std::shared_ptr<std::vector<int>> faces =
      std::make_shared<std::vector<int>>();

  // The vertex is simple
  for (std::size_t i = 0; i < points.size(); i++) {
    vertices->push_back(
        fcl::Vector3<T>(points[i][0], points[i][1], points[i][2]));
  }

  // The face is more complex
  for (std::size_t i = 0; i < triangles.size(); i++) {
    faces->push_back(3);
    faces->push_back(triangles[i][0]);
    faces->push_back(triangles[i][1]);
    faces->push_back(triangles[i][2]);
  }

  fcl::Convex<T> convex(vertices, num_faces, faces, true);
  return convex;
}

template <typename T>
void fcl::test::meshEllipsoidTriangles(
    T size_x, T size_y, T size_z, std::vector<std::array<T, 3>>* points,
    std::vector<std::array<int, 3>>* triangles) {
  T v[58][3] = {
      {0.135299, -0.461940, -0.135299},  {0.000000, -0.461940, -0.191342},
      {-0.135299, -0.461940, -0.135299}, {-0.191342, -0.461940, 0.000000},
      {-0.135299, -0.461940, 0.135299},  {0.000000, -0.461940, 0.191342},
      {0.135299, -0.461940, 0.135299},   {0.191342, -0.461940, 0.000000},
      {0.250000, -0.353553, -0.250000},  {0.000000, -0.353553, -0.353553},
      {-0.250000, -0.353553, -0.250000}, {-0.353553, -0.353553, 0.000000},
      {-0.250000, -0.353553, 0.250000},  {0.000000, -0.353553, 0.353553},
      {0.250000, -0.353553, 0.250000},   {0.353553, -0.353553, 0.000000},
      {0.326641, -0.191342, -0.326641},  {0.000000, -0.191342, -0.461940},
      {-0.326641, -0.191342, -0.326641}, {-0.461940, -0.191342, 0.000000},
      {-0.326641, -0.191342, 0.326641},  {0.000000, -0.191342, 0.461940},
      {0.326641, -0.191342, 0.326641},   {0.461940, -0.191342, 0.000000},
      {0.353553, 0.000000, -0.353553},   {0.000000, 0.000000, -0.500000},
      {-0.353553, 0.000000, -0.353553},  {-0.500000, 0.000000, 0.000000},
      {-0.353553, 0.000000, 0.353553},   {0.000000, 0.000000, 0.500000},
      {0.353553, 0.000000, 0.353553},    {0.500000, 0.000000, 0.000000},
      {0.326641, 0.191342, -0.326641},   {0.000000, 0.191342, -0.461940},
      {-0.326641, 0.191342, -0.326641},  {-0.461940, 0.191342, 0.000000},
      {-0.326641, 0.191342, 0.326641},   {0.000000, 0.191342, 0.461940},
      {0.326641, 0.191342, 0.326641},    {0.461940, 0.191342, 0.000000},
      {0.250000, 0.353553, -0.250000},   {0.000000, 0.353553, -0.353553},
      {-0.250000, 0.353553, -0.250000},  {-0.353553, 0.353553, 0.000000},
      {-0.250000, 0.353553, 0.250000},   {0.000000, 0.353553, 0.353553},
      {0.250000, 0.353553, 0.250000},    {0.353553, 0.353553, 0.000000},
      {0.135299, 0.461940, -0.135299},   {0.000000, 0.461940, -0.191342},
      {-0.135299, 0.461940, -0.135299},  {-0.191342, 0.461940, 0.000000},
      {-0.135299, 0.461940, 0.135299},   {0.000000, 0.461940, 0.191342},
      {0.135299, 0.461940, 0.135299},    {0.191342, 0.461940, 0.000000},
      {0.000000, -0.500000, 0.000000},   {0.000000, 0.500000, 0.000000}};

  int f[112][3] = {
      {1, 2, 9},    {9, 2, 10},   {2, 3, 10},   {10, 3, 11},  {3, 4, 11},
      {11, 4, 12},  {4, 5, 12},   {12, 5, 13},  {5, 6, 13},   {13, 6, 14},
      {6, 7, 14},   {14, 7, 15},  {7, 8, 15},   {15, 8, 16},  {8, 1, 16},
      {16, 1, 9},   {9, 10, 17},  {17, 10, 18}, {10, 11, 18}, {18, 11, 19},
      {11, 12, 19}, {19, 12, 20}, {12, 13, 20}, {20, 13, 21}, {13, 14, 21},
      {21, 14, 22}, {14, 15, 22}, {22, 15, 23}, {15, 16, 23}, {23, 16, 24},
      {16, 9, 24},  {24, 9, 17},  {17, 18, 25}, {25, 18, 26}, {18, 19, 26},
      {26, 19, 27}, {19, 20, 27}, {27, 20, 28}, {20, 21, 28}, {28, 21, 29},
      {21, 22, 29}, {29, 22, 30}, {22, 23, 30}, {30, 23, 31}, {23, 24, 31},
      {31, 24, 32}, {24, 17, 32}, {32, 17, 25}, {25, 26, 33}, {33, 26, 34},
      {26, 27, 34}, {34, 27, 35}, {27, 28, 35}, {35, 28, 36}, {28, 29, 36},
      {36, 29, 37}, {29, 30, 37}, {37, 30, 38}, {30, 31, 38}, {38, 31, 39},
      {31, 32, 39}, {39, 32, 40}, {32, 25, 40}, {40, 25, 33}, {33, 34, 41},
      {41, 34, 42}, {34, 35, 42}, {42, 35, 43}, {35, 36, 43}, {43, 36, 44},
      {36, 37, 44}, {44, 37, 45}, {37, 38, 45}, {45, 38, 46}, {38, 39, 46},
      {46, 39, 47}, {39, 40, 47}, {47, 40, 48}, {40, 33, 48}, {48, 33, 41},
      {41, 42, 49}, {49, 42, 50}, {42, 43, 50}, {50, 43, 51}, {43, 44, 51},
      {51, 44, 52}, {44, 45, 52}, {52, 45, 53}, {45, 46, 53}, {53, 46, 54},
      {46, 47, 54}, {54, 47, 55}, {47, 48, 55}, {55, 48, 56}, {48, 41, 56},
      {56, 41, 49}, {2, 1, 57},   {3, 2, 57},   {4, 3, 57},   {5, 4, 57},
      {6, 5, 57},   {7, 6, 57},   {8, 7, 57},   {1, 8, 57},   {49, 50, 58},
      {50, 51, 58}, {51, 52, 58}, {52, 53, 58}, {53, 54, 58}, {54, 55, 58},
      {55, 56, 58}, {56, 49, 58}};

  // Append to result
  points->clear();
  for (auto i = 0; i < 58; i++) {
    std::array<T, 3> point_i{v[i][0] * size_x, v[i][1] * size_y,
                             v[i][2] * size_z};
    points->emplace_back(point_i);
  }

  triangles->clear();
  for (auto i = 0; i < 112; i++) {
    std::array<int, 3> triangle_i{f[i][0] - 1, f[i][1] - 1, f[i][2] - 1};
    triangles->emplace_back(triangle_i);
  }
}

template <typename T>
void fcl::test::meshEllipsoidUV(T size_x, T size_y, T size_z, int n_slices,
                                int n_stacks,
                                std::vector<std::array<T, 3>>* points,
                                std::vector<std::array<int, 3>>* triangles) {
  constexpr T local_pi = 3.1415926;
  using PointT = std::array<T, 3>;
  using TriangleT = std::array<int, 3>;
  points->clear();
  points->emplace_back(PointT{0, 0, T(1.0)});

  // generate vertices per stack / slice
  for (int i = 0; i < n_stacks - 1; i++) {
    T phi = local_pi * T(i + 1) / T(n_stacks);
    for (int j = 0; j < n_slices; j++) {
      const T theta = 2.0 * local_pi * T(j) / T(n_slices);
      T x = std::sin(phi) * std::cos(theta);
      T y = std::sin(phi) * std::sin(theta);
      T z = std::cos(phi);
      points->emplace_back(PointT{x, y, z});
    }
  }

  // The bottom
  points->emplace_back(PointT{0, 0, T(-1.0)});
  int bottom_idx = static_cast<int>(points->size() - 1);

  // Scale the points
  for (std::size_t i = 0; i < points->size(); i++) {
    points->at(i)[0] *= (0.5 * size_x);
    points->at(i)[1] *= (0.5 * size_y);
    points->at(i)[2] *= (0.5 * size_z);
  }

  // add top / bottom triangles
  triangles->clear();
  for (int i = 0; i < n_slices; ++i) {
    auto i0 = i + 1;
    auto i1 = (i + 1) % n_slices + 1;
    // mesh.add_triangle(v0, Vertex(i1), Vertex(i0));
    triangles->emplace_back(TriangleT{0, i1, i0});
    i0 = i + n_slices * (n_stacks - 2) + 1;
    i1 = (i + 1) % n_slices + n_slices * (n_stacks - 2) + 1;
    // mesh.add_triangle(v1, Vertex(i0), Vertex(i1));
    triangles->emplace_back(TriangleT{bottom_idx, i0, i1});
  }

  // add quads per stack / slice
  for (int j = 0; j < n_stacks - 2; j++) {
    auto j0 = j * n_slices + 1;
    auto j1 = (j + 1) * n_slices + 1;
    for (int i = 0; i < n_slices; i++) {
      auto i0 = j0 + i;
      auto i1 = j0 + (i + 1) % n_slices;
      auto i2 = j1 + (i + 1) % n_slices;
      auto i3 = j1 + i;
      // mesh.add_quad(Vertex(i0), Vertex(i1), Vertex(i2), Vertex(i3));
      triangles->emplace_back(TriangleT{i0, i1, i2});
      triangles->emplace_back(TriangleT{i2, i3, i0});
    }
  }
}