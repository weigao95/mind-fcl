//
// Created by wei on 22-6-5.
//

#pragma once

#include "epa_polytope_utils.h"
#include "minkowski_diff.h"

namespace fcl {
namespace cvx_collide {

enum class PolytopeElementType {
  Invalid,
  ListHeader,
  FreeElementHeader,
  Vertex,
  Edge,
  Face
};

// For visibility computation
enum class EdgeVisibilityStatus {
  Unknown,  // not computed yet
  Border,   // edges for which only *one* adjacent face is in the patch
  Internal  // the edges for which *both* adjacent faces are in the patch
};
enum class FaceVisibilityStatus {
  Unknown,  // not computed yet
  Visible,  // visible from a point (thus, should be removed)
  Hidden    // not visible
};

// Base class for all polytope vertex/edge/face
struct PolytopeElementBase {
  // Basic meta
  int id;
  PolytopeElementType type;

  // The (single) linked list element for members with same class
  PolytopeElementBase* next;

  // The (single) linked list element for malloc
  PolytopeElementBase* malloc_next;

  // Construct as invalid
  explicit PolytopeElementBase()
      : id(-1),
        type(PolytopeElementType::Invalid),
        next(nullptr),
        malloc_next(nullptr) {}

  // Clear a potentially valid record
  void Clear() {
    id = -1;
    type = PolytopeElementType::Invalid;
    next = nullptr;
    // Do NOT touch malloc_next
  }
};

template <typename T>
struct PolytopeVertex : public PolytopeElementBase {
  Vector3<T> vertex_location;
  Vector3<T> minkowski_diff_direction;
  MinDistanceToSimplex<T> min_distance;

  // For visible computation cache
  mutable bool cached_vertex_should_remove{false};
  mutable PolytopeElementBase* cached_vertex2new_v_edge{nullptr};
};

template <typename T>
struct PolytopeFace;

template <typename T>
struct PolytopeEdge : public PolytopeElementBase {
  std::array<PolytopeVertex<T>*, 2> vertices_of_edge;
  std::array<PolytopeFace<T>*, 2> faces_of_edge;
  MinDistanceToSimplex<T> min_distance;

  // For visible computation cache
  mutable EdgeVisibilityStatus cached_visibility{EdgeVisibilityStatus::Unknown};
};

template <typename T>
struct PolytopeFace : public PolytopeElementBase {
  std::array<PolytopeEdge<T>*, 3> edges_of_face;
  MinDistanceToSimplex<T> min_distance;

  // For visible computation cache
  mutable FaceVisibilityStatus cached_visibility{FaceVisibilityStatus::Unknown};

  // Helpers to get the vertices of the face
  bool GetVertices(Vector3<T>& a, Vector3<T>& b, Vector3<T>& c) const;
  bool GetVertices(MinkowskiDiffVertex<T>& a, MinkowskiDiffVertex<T>& b,
                   MinkowskiDiffVertex<T>& c) const;

 private:
  bool getRawVertices(PolytopeVertex<T>** a, PolytopeVertex<T>** b,
                      PolytopeVertex<T>** c) const;
};

template <typename T>
class Polytope {
 private:
  // The valid element as a list
  using PolytopeElementList = PolytopeElementBase;
  PolytopeElementList vertex_list_;
  PolytopeElementList edge_list_;
  PolytopeElementList face_list_;

  // The buffer in the polytope
  std::vector<PolytopeVertex<T>> vertex_buffer_;
  std::vector<PolytopeEdge<T>> edge_buffer_;
  std::vector<PolytopeFace<T>> face_buffer_;

  // The id for vertex/edge/face, they are not necessarily index in buffer
  std::size_t next_vertex_id_{0};
  std::size_t next_edge_id_{0};
  std::size_t next_face_id_{0};

  // The freed buffer list
  using PolytopeElementMallocList = PolytopeElementBase;
  PolytopeElementMallocList vertex_malloc_list_;
  PolytopeElementMallocList edge_malloc_list_;
  PolytopeElementMallocList face_malloc_list_;

  // The construction/re-initialization of the polytope
 public:
  explicit Polytope(std::size_t max_n_faces) { Reset(max_n_faces); }
  ~Polytope() = default;
  // No copy/assign
  Polytope(const Polytope<T>&) = delete;
  Polytope& operator=(const Polytope<T>&) = delete;

  void Reset(std::size_t max_n_faces);
  std::size_t face_capacity() const { return face_buffer_.size(); }
  bool empty() const { return isEmptyList(vertex_list_); }

 private:
  // The management of list
  static bool isEmptyList(const PolytopeElementList& list);
  static void resetElementList(PolytopeElementList& list);
  template <typename Derived>
  static void resetBufferAndMallocList(std::vector<Derived>& buffer_vector,
                                       PolytopeElementMallocList& malloc_list);
  static void addElementToList(PolytopeElementList& list,
                               PolytopeElementBase* element);
  static void visitList(
      PolytopeElementList& list,
      const std::function<void(PolytopeElementBase*)>& visitor);
  static void visitListConst(
      const PolytopeElementList& list,
      const std::function<void(const PolytopeElementBase*)>& visitor);
  bool VisitElementList(
      PolytopeElementType element_type,
      const std::function<void(const PolytopeElementBase*)>& visitor) const;
  static std::size_t removeIf(
      PolytopeElementList& list,
      const std::function<bool(PolytopeElementBase*)>& remove_if_true,
      PolytopeElementMallocList& reclaim_removed_list);

  // Add/remove a new vertex/edge/face
 public:
  PolytopeVertex<T>* AddNewVertex(const MinkowskiDiffVertex<T>& vertex);
  PolytopeEdge<T>* AddNewEdge(PolytopeVertex<T>* v1, PolytopeVertex<T>* v2);
  PolytopeFace<T>* AddNewFace(PolytopeEdge<T>* e1, PolytopeEdge<T>* e2,
                              PolytopeEdge<T>* e3);

 public:
  // Iterate over point/segment/triangle, get the one with min-distance to
  // the origin. Prefer point to segment, and segment to triangle.
  PolytopeElementBase* ComputeMinDistanceToOrigin(bool exclude_vertex = false);

  // Compute the normal of face direction that is pointed outside
  bool ComputeFaceNormalPointingOutward(const PolytopeFace<T>* face,
                                        Vector3<T>& normal,
                                        T* area = nullptr) const;

  // Expand the polytope
  enum class PolytopeExpansionStatus { OK, Failed, MallocFailed };
  PolytopeExpansionStatus ExpandPolytope(const MinkowskiDiffVertex<T>& next_v,
                                         PolytopeFace<T>* start_face);
  bool IsPointOutsidePolytopeFace(
      const PolytopeFace<T>* face, const Vector3<T>& point,
      const T point2plane_should_greater_than = T(0.0)) const;

 private:
  void initVisibilityCacheVariables();
  bool computeVisiblePatch(PolytopeFace<T>* start_face,
                           const Vector3<T>& next_v);

  // Given visibility result of edge/face, compute the vertices to remove
  void updateVertexRemoveFlag();
  void removeAccordingToVisibility();
};

}  // namespace cvx_collide
}  // namespace fcl

#include "epa_polytope.hpp"
#include "epa_polytope_expand.hpp"
