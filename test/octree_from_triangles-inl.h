#pragma once

namespace fcl {

template <typename S>
void samplePointsOnTriangles(
    const TriangleP<S>& triangle, S resolution,
    const std::function<void(const Vector3<S>&)>& sampled_point_processor) {
  const auto& raw_a = triangle.a;
  const auto& raw_b = triangle.b;
  const auto& raw_c = triangle.c;

  // Reorder such that ab > bc and ac > bc
  std::array<Vector3<S>, 3> triangle_points{raw_a, raw_b, raw_c};
  std::array<Vector3<S>, 3> triangle_edges{raw_b - raw_a, raw_c - raw_b,
                                           raw_a - raw_c};

  // Compute the edge with smallest length
  int min_length_edge_idx = -1;
  S min_edge_length = -1;
  for (auto i = 0; i < 3; i++) {
    const auto edge_length_i = triangle_edges[i].norm();
    if (i == 0 || edge_length_i < min_edge_length) {
      min_edge_length = edge_length_i;
      min_length_edge_idx = i;
    }
  }

  // Make b be the one corresponds to min_length_edge
  // such that ab >= bc, ac >= bc
  const auto a_index = (min_length_edge_idx + 2) % 3;
  const auto b_index = (min_length_edge_idx + 0) % 3;
  const auto c_index = (min_length_edge_idx + 1) % 3;
  const auto& a = triangle_points[a_index];
  const auto& b = triangle_points[b_index];
  const auto& c = triangle_points[c_index];
  const auto ab_length = (a - b).norm();
  const auto ac_length = (a - c).norm();
  const auto bc_length = min_edge_length;
  assert(std::abs((c - b).norm() - min_edge_length) <= 1e-5);
  assert(ab_length >= min_edge_length);
  assert(ac_length >= min_edge_length);

  const int n_ab_division =
      std::max<int>(2, static_cast<int>(std::ceil(ab_length / resolution)));
  const S ab_delta = 1.0 / (n_ab_division - 1);
  const int n_ac_division =
      std::max<int>(2, static_cast<int>(std::ceil(ac_length / resolution)));
  const S ac_delta = 1.0 / (n_ac_division - 1);
  const int n_bc_division =
      std::max<int>(2, static_cast<int>(std::ceil(bc_length / resolution)));
  const S bc_delta = 1.0 / (n_bc_division - 1);

  // If min_edge_length is too small, then do not sample in interior
  // Instead, iterate through the two longer edge
  if (min_edge_length <= resolution) {
    for (int i = 0; i < n_ab_division; i++) {
      const S alpha_i = i * ab_delta;
      const auto p = alpha_i * b + (1.0 - alpha_i) * a;
      sampled_point_processor(p);
    }

    for(int j = 1; j < n_ac_division; j++) {
      const S beta_j = j * ac_delta;
      const auto p = beta_j * c + (1.0 - beta_j) * a;
      sampled_point_processor(p);
    }

    // Done
    return;
  }

  // Iterate into the interior, use coordinate that contains the smallest edge and largest edge
  if (ab_length > ac_length) {
    const auto n_ba_division = n_ab_division;
    const auto ba_delta = ab_delta;
    for (int i = 0; i < n_ba_division; i++) {
      const S alpha_i = i * ba_delta;
      for (int j = 0; j < n_bc_division; j++) {
        const S beta_j = j * bc_delta;
        const S gamma_ij = 1.0 - alpha_i - beta_j;
        if (gamma_ij < 0) {
          continue;
        }

        // Now, compute the point
        assert(gamma_ij >= 0);
        const auto p_ij = alpha_i * a + beta_j * c + gamma_ij * b;
        sampled_point_processor(p_ij);
      }
    }
  } else {
    const auto n_ca_division = n_ac_division;
    const auto ca_delta = ac_delta;
    const auto n_cb_division = n_bc_division;
    const auto cb_delta = bc_delta;
    for (int i = 0; i < n_ca_division; i++) {
      const S alpha_i = i * ca_delta;
      for (int j = 0; j < n_cb_division; j++) {
        const S beta_j = j * cb_delta;
        const S gamma_ij = 1.0 - alpha_i - beta_j;
        if (gamma_ij < 0) {
          continue;
        }

        // Now, compute the point
        assert(gamma_ij >= 0);
        const auto p_ij = alpha_i * a + beta_j * b + gamma_ij * c;
        sampled_point_processor(p_ij);
      }
    }
  }
}

template <typename S>
void makeTriangleSoupSurfaceVoxel(
    const std::vector<TriangleP<S>>& triangles, S resolution,
    std::unordered_set<std::int64_t>& encoded_voxel_set) {
  encoded_voxel_set.clear();
  auto functor = [&encoded_voxel_set,
                  &resolution](const Vector3<S>& point) -> void {
    const auto px = point[0] / resolution;
    const auto py = point[1] / resolution;
    const auto pz = point[2] / resolution;
    const std::int16_t px_int16 =
        static_cast<std::int16_t>(px >= 0 ? px : (px - 1));
    const std::int16_t py_int16 =
        static_cast<std::int16_t>(py >= 0 ? py : (py - 1));
    const std::int16_t pz_int16 =
        static_cast<std::int16_t>(pz >= 0 ? pz : (pz - 1));
    const std::int64_t encoded = (static_cast<std::int64_t>(px_int16) << 0) +
                                 (static_cast<std::int64_t>(py_int16) << 16) +
                                 (static_cast<std::int64_t>(pz_int16) << 32);
    encoded_voxel_set.insert(encoded);
  };

  for (auto iter = triangles.cbegin(); iter != triangles.cend(); iter++) {
    samplePointsOnTriangles<S>(*iter, resolution, functor);
  }
}

}  // namespace fcl