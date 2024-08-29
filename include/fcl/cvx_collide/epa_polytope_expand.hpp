//
// Created by wei on 22-6-5.
//

#pragma once

#include "epa_polytope_utils.h"

namespace fcl {
namespace cvx_collide {

template <typename T>
bool Polytope<T>::IsPointOutsidePolytopeFace(
    const PolytopeFace<T>* face, const Vector3<T>& point,
    const T point2plane_should_greater_than) const {
  Vector3<T> face_normal;
  T area{0};
  auto ok = ComputeFaceNormalPointingOutward(face, face_normal, &area);
  if (!ok) {
    // If area is zero, then also remove it
    if (area <= T(0.0))
      return true;
    else
      return false;
  }
  const Vector3<T>& point_on_face =
      face->edges_of_face[0]->vertices_of_edge[0]->vertex_location;
  const T dot_value = face_normal.dot(point - point_on_face);
  return dot_value >= point2plane_should_greater_than;
}

template <typename T>
typename Polytope<T>::PolytopeExpansionStatus Polytope<T>::ExpandPolytope(
    const MinkowskiDiffVertex<T>& next_v, PolytopeFace<T>* start_face) {
  // If we are here, we should have a valid start face
  assert(start_face != nullptr);

  // First mark edge/face visibility as unknown
  initVisibilityCacheVariables();
  auto ok = computeVisiblePatch(start_face, next_v.vertex);
  if (!ok) return PolytopeExpansionStatus::Failed;

  // Given visibility of faces/edges, remove the vertices
  updateVertexRemoveFlag();

  // Now, remove all visible faces and internal edges
  removeAccordingToVisibility();

  // Add new vertex
  PolytopeVertex<T>* new_v = AddNewVertex(next_v);
  if (new_v == nullptr) return PolytopeExpansionStatus::MallocFailed;

  // Add new faces
  bool new_face_ok = true;
  auto add_new_edge_face_functor =
      [this, &new_v, &new_face_ok](PolytopeElementBase* element) -> void {
    auto type = element->type;
    if (type != PolytopeElementType::Edge) return;
    auto* edge = static_cast<PolytopeEdge<T>*>(element);
    if (edge->cached_visibility != EdgeVisibilityStatus::Border) return;

    // Two other edges
    std::array<PolytopeEdge<T>*, 2> e{nullptr, nullptr};
    for (auto i = 0; i < 2; i++) {
      PolytopeVertex<T>* v_i = edge->vertices_of_edge[i];
      if (v_i->cached_vertex2new_v_edge == nullptr) {
        v_i->cached_vertex2new_v_edge = AddNewEdge(new_v, v_i);
      }

      // Now the edge should be valid
      e[i] = static_cast<PolytopeEdge<T>*>(v_i->cached_vertex2new_v_edge);
    }

    // Make the face
    auto* f = AddNewFace(edge, e[0], e[1]);
    if (f == nullptr) {
      new_face_ok = false;
      return;
    }
  };
  visitList(edge_list_, add_new_edge_face_functor);

  // Return flag
  return new_face_ok ? PolytopeExpansionStatus::OK
                     : PolytopeExpansionStatus::MallocFailed;
}

template <typename T>
void Polytope<T>::initVisibilityCacheVariables() {
  auto visit_element_functor = [](PolytopeElementBase* element) -> void {
    // Now it holds a distance
    switch (element->type) {
      case PolytopeElementType::Vertex: {
        auto* vertex = static_cast<PolytopeVertex<T>*>(element);
        // mark as true, invalidate later
        vertex->cached_vertex_should_remove = true;
        vertex->cached_vertex2new_v_edge = nullptr;
        return;
      }
      case PolytopeElementType::Edge: {
        auto* edge = static_cast<PolytopeEdge<T>*>(element);
        edge->cached_visibility = EdgeVisibilityStatus::Unknown;
        return;
      }
      default: {
        assert(element->type == PolytopeElementType::Face);
        auto* face = static_cast<PolytopeFace<T>*>(element);
        face->cached_visibility = FaceVisibilityStatus::Unknown;
        return;
      }
    }
  };

  // Go
  visitList(vertex_list_, visit_element_functor);
  visitList(edge_list_, visit_element_functor);
  visitList(face_list_, visit_element_functor);
}

template <typename T>
bool Polytope<T>::computeVisiblePatch(PolytopeFace<T>* start_face,
                                      const Vector3<T>& next_v) {
  struct Task {
    PolytopeFace<T>* face;
    int edge_idx_of_face;
  };
  std::stack<Task> task_stack;

  // Init
  start_face->cached_visibility = FaceVisibilityStatus::Visible;
  for (auto edge_idx = 0; edge_idx < 3; edge_idx++) {
    task_stack.push({start_face, edge_idx});
  }

  // Start the loop
  while (!task_stack.empty()) {
    // Get this task
    const auto& task = task_stack.top();
    PolytopeFace<T>* f = task.face;
    int edge_index = task.edge_idx_of_face;
    task_stack.pop();

    assert(f->cached_visibility == FaceVisibilityStatus::Visible);
    PolytopeEdge<T>* edge = f->edges_of_face[edge_index];
    if (edge == nullptr) return false;
    PolytopeFace<T>* g = (edge->faces_of_edge[0] == f) ? edge->faces_of_edge[1]
                                                       : edge->faces_of_edge[0];
    if (g == nullptr) return false;
    assert(g != nullptr);
    bool is_visible = (g->cached_visibility == FaceVisibilityStatus::Visible);
    bool is_hidden = (g->cached_visibility == FaceVisibilityStatus::Hidden);
    assert(!(is_visible && is_hidden));

    if (is_visible) {
      // Both f and g are visible
      edge->cached_visibility = EdgeVisibilityStatus::Internal;
      continue;
    }

    if (is_hidden) {
      edge->cached_visibility = EdgeVisibilityStatus::Border;
      continue;
    }

    // face g is not classified yet
    is_visible = IsPointOutsidePolytopeFace(g, next_v);
    if (is_visible) {
      g->cached_visibility = FaceVisibilityStatus::Visible;
      edge->cached_visibility = EdgeVisibilityStatus::Internal;
      for (auto i = 0; i < 3; i++) {
        task_stack.push({g, i});
      }
    } else {
      g->cached_visibility = FaceVisibilityStatus::Hidden;
      edge->cached_visibility = EdgeVisibilityStatus::Border;
    }
  }

  // Expand ok
  return true;
}

template <typename T>
void Polytope<T>::updateVertexRemoveFlag() {
  auto visit_element_functor = [](PolytopeElementBase* element) -> void {
    // Now it holds a distance
    switch (element->type) {
      case PolytopeElementType::Edge: {
        auto* edge = static_cast<PolytopeEdge<T>*>(element);
        if (edge->cached_visibility != EdgeVisibilityStatus::Internal) {
          // do not remove the edge, thus must keep all vertices
          edge->vertices_of_edge[0]->cached_vertex_should_remove = false;
          edge->vertices_of_edge[1]->cached_vertex_should_remove = false;
        }
        return;
      }
      default: {
        return;
      }
    }
  };

  // Go
  visitList(edge_list_, visit_element_functor);
}

template <typename T>
std::size_t Polytope<T>::removeIf(
    PolytopeElementList& list,
    const std::function<bool(PolytopeElementBase*)>& remove_if_true,
    PolytopeElementMallocList& reclaim_removed_list) {
  // Nothing removed
  if (isEmptyList(list)) return 0;

  // Now, at least one element
  std::size_t remove_count = 0;
  PolytopeElementBase* prev = &list;
  PolytopeElementBase* current = prev->next;
  assert(current != nullptr);
  while (current != nullptr) {
    bool should_remove = remove_if_true(current);
    if (should_remove) {
      // Update the connection
      prev->next = current->next;

      // Clear and reclaim this one
      current->Clear();
      current->malloc_next = reclaim_removed_list.malloc_next;
      reclaim_removed_list.malloc_next = current;

      // Move to next
      remove_count += 1;
      current = prev->next;
    } else {
      prev = current;
      current = prev->next;
    }
  }

  return remove_count;
}

template <typename T>
void Polytope<T>::removeAccordingToVisibility() {
  // First remove face
  auto should_face_removed = [](PolytopeElementBase* element) -> bool {
    assert(element->type == PolytopeElementType::Face);
    auto* face = static_cast<const PolytopeFace<T>*>(element);
    auto should_remove =
        (face->cached_visibility == FaceVisibilityStatus::Visible);

    // Need to remove the reference in edges
    if (should_remove) {
      for (auto i = 0; i < 3; i++) {
        PolytopeEdge<T>* edge_i = face->edges_of_face[i];
        if (edge_i->faces_of_edge[0] == face) {
          edge_i->faces_of_edge[0] = edge_i->faces_of_edge[1];
        }
        edge_i->faces_of_edge[1] = nullptr;
      }
    }
    return should_remove;
  };
  removeIf(face_list_, should_face_removed, face_malloc_list_);

  // Next remove edge
  auto should_edge_removed = [](PolytopeElementBase* element) -> bool {
    assert(element->type == PolytopeElementType::Edge);
    auto* edge = static_cast<const PolytopeEdge<T>*>(element);
    return edge->cached_visibility == EdgeVisibilityStatus::Internal;
  };
  removeIf(edge_list_, should_edge_removed, edge_malloc_list_);

  // Finally, remove vertex
  auto should_vertex_removed = [](PolytopeElementBase* element) -> bool {
    assert(element->type == PolytopeElementType::Vertex);
    auto* vertex = static_cast<const PolytopeVertex<T>*>(element);
    return vertex->cached_vertex_should_remove;
  };
  removeIf(vertex_list_, should_vertex_removed, vertex_malloc_list_);
}

}  // namespace cvx_collide
}  // namespace fcl