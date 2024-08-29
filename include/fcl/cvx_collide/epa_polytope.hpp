//
// Created by wei on 22-6-5.
//

#pragma once

namespace fcl {
namespace cvx_collide {

template <typename T>
bool PolytopeFace<T>::getRawVertices(PolytopeVertex<T>** a,
                                     PolytopeVertex<T>** b,
                                     PolytopeVertex<T>** c) const {
  if (edges_of_face[0] == nullptr || edges_of_face[1] == nullptr ||
      edges_of_face[2] == nullptr)
    return false;

  // Get two edges
  PolytopeEdge<T>* e0 = edges_of_face[0];
  PolytopeVertex<T>* v0 = e0->vertices_of_edge[0];
  PolytopeVertex<T>* v1 = e0->vertices_of_edge[1];
  if (v0 == nullptr || v1 == nullptr) return false;

  // The remaining vertices
  PolytopeEdge<T>* e1 = edges_of_face[1];
  PolytopeVertex<T>* v2 = nullptr;
  if (e1->vertices_of_edge[0] != v0 && e1->vertices_of_edge[0] != v1) {
    v2 = e1->vertices_of_edge[0];
  } else {
    v2 = e1->vertices_of_edge[1];
  }
  if (v2 == nullptr) return false;

  // Assign output
  *a = v0;
  *b = v1;
  *c = v2;
  return true;
}

template <typename T>
bool PolytopeFace<T>::GetVertices(Vector3<T>& a, Vector3<T>& b,
                                  Vector3<T>& c) const {
  PolytopeVertex<T>* v0;
  PolytopeVertex<T>* v1;
  PolytopeVertex<T>* v2;
  auto ok = getRawVertices(&v0, &v1, &v2);
  if (!ok) return false;

  // Get vertex location
  a = v0->vertex_location;
  b = v1->vertex_location;
  c = v2->vertex_location;
  return true;
}

template <typename T>
bool PolytopeFace<T>::GetVertices(MinkowskiDiffVertex<T>& a,
                                  MinkowskiDiffVertex<T>& b,
                                  MinkowskiDiffVertex<T>& c) const {
  PolytopeVertex<T>* v0;
  PolytopeVertex<T>* v1;
  PolytopeVertex<T>* v2;
  auto ok = getRawVertices(&v0, &v1, &v2);
  if (!ok) return false;

  // Get vertex location/direction
  a.vertex = v0->vertex_location;
  b.vertex = v1->vertex_location;
  c.vertex = v2->vertex_location;
  a.direction = v0->minkowski_diff_direction;
  b.direction = v1->minkowski_diff_direction;
  c.direction = v2->minkowski_diff_direction;
  return true;
}

template <typename T>
void Polytope<T>::Reset(std::size_t max_n_faces) {
  // Reset the list element
  resetElementList(vertex_list_);
  resetElementList(edge_list_);
  resetElementList(face_list_);

  // Clear the current elements
  next_vertex_id_ = 0;
  next_edge_id_ = 0;
  next_face_id_ = 0;

  // Reserve new elements
  const auto max_n_edges = static_cast<std::size_t>(1.51 * double(max_n_faces));
  const auto max_n_vertices =
      static_cast<std::size_t>(1.51 * double(max_n_faces));
  vertex_buffer_.resize(max_n_vertices);
  edge_buffer_.resize(max_n_edges);
  face_buffer_.resize(max_n_faces);

  // Add to malloc list
  resetBufferAndMallocList<PolytopeVertex<T>>(vertex_buffer_,
                                              vertex_malloc_list_);
  resetBufferAndMallocList<PolytopeEdge<T>>(edge_buffer_, edge_malloc_list_);
  resetBufferAndMallocList<PolytopeFace<T>>(face_buffer_, face_malloc_list_);
}

template <typename T>
bool Polytope<T>::isEmptyList(const PolytopeElementList& list) {
  return (list.next == nullptr);
}

template <typename T>
void Polytope<T>::resetElementList(PolytopeElementList& list) {
  list.id = -1;
  list.type = PolytopeElementType::ListHeader;
  list.next = nullptr;
}

template <typename T>
template <typename Derived>
void Polytope<T>::resetBufferAndMallocList(
    std::vector<Derived>& buffer_vector,
    PolytopeElementMallocList& malloc_list) {
  assert(!buffer_vector.empty());
  malloc_list.Clear();
  malloc_list.type = PolytopeElementType::FreeElementHeader;
  malloc_list.malloc_next = nullptr;
  for (std::size_t i = 0; i < buffer_vector.size(); i++) {
    std::size_t reversed_idx = buffer_vector.size() - i - 1;
    buffer_vector[reversed_idx].Clear();
    buffer_vector[reversed_idx].malloc_next = malloc_list.malloc_next;
    malloc_list.malloc_next = &buffer_vector[reversed_idx];
  }
}

template <typename T>
void Polytope<T>::addElementToList(PolytopeElementList& list,
                                   PolytopeElementBase* element) {
  element->next = list.next;
  list.next = element;
}

template <typename T>
void Polytope<T>::visitList(
    PolytopeElementList& list,
    const std::function<void(PolytopeElementBase*)>& visitor) {
  PolytopeElementBase* current = list.next;
  while (current != nullptr) {
    visitor(current);
    current = current->next;
  }
}

template <typename T>
void Polytope<T>::visitListConst(
    const PolytopeElementList& list,
    const std::function<void(const PolytopeElementBase*)>& visitor) {
  PolytopeElementBase* current = list.next;
  while (current != nullptr) {
    visitor(current);
    current = current->next;
  }
}

template <typename T>
bool Polytope<T>::VisitElementList(
    PolytopeElementType element_type,
    const std::function<void(const PolytopeElementBase*)>& visitor) const {
  if (element_type == PolytopeElementType::Vertex) {
    visitListConst(vertex_list_, visitor);
    return true;
  } else if (element_type == PolytopeElementType::Edge) {
    visitListConst(edge_list_, visitor);
    return true;
  } else if (element_type == PolytopeElementType::Face) {
    visitListConst(face_list_, visitor);
    return true;
  } else {
    return false;
  }
}

template <typename T>
PolytopeVertex<T>* Polytope<T>::AddNewVertex(
    const MinkowskiDiffVertex<T>& vertex) {
  // Do malloc
  if (vertex_malloc_list_.malloc_next == nullptr) return nullptr;
  auto* v = static_cast<PolytopeVertex<T>*>(vertex_malloc_list_.malloc_next);
  vertex_malloc_list_.malloc_next = v->malloc_next;

  // Init the vertex
  v->id = next_vertex_id_;
  next_vertex_id_ += 1;
  v->type = PolytopeElementType::Vertex;
  v->vertex_location = vertex.vertex;
  v->minkowski_diff_direction = vertex.direction;

  // The min-distance
  v->min_distance.witness_in_simplex = true;
  v->min_distance.witness = vertex.vertex;
  v->min_distance.min_distance_square = vertex.vertex.squaredNorm();

  // Add to vertex
  addElementToList(vertex_list_, v);

  // Claim this one
  return v;
}

template <typename T>
PolytopeEdge<T>* Polytope<T>::AddNewEdge(PolytopeVertex<T>* v1,
                                         PolytopeVertex<T>* v2) {
  // Do not continue if input invalid
  if (v1 == nullptr || v2 == nullptr) return nullptr;

  // Do malloc
  if (edge_malloc_list_.malloc_next == nullptr) return nullptr;
  auto* edge = static_cast<PolytopeEdge<T>*>(edge_malloc_list_.malloc_next);
  edge_malloc_list_.malloc_next = edge->malloc_next;

  // Init the edge
  edge->id = next_edge_id_;
  next_edge_id_ += 1;
  edge->type = PolytopeElementType::Edge;
  edge->vertices_of_edge[0] = v1;
  edge->vertices_of_edge[1] = v2;
  edge->faces_of_edge[0] = edge->faces_of_edge[1] = nullptr;

  // The min distance
  const Vector3<T>& a = v1->vertex_location;
  const Vector3<T>& b = v2->vertex_location;
  pointToSegmentSquaredDistance<T>(Vector3<T>::Zero(), a, b,
                                   edge->min_distance);

  // add to edge
  addElementToList(edge_list_, edge);
  return edge;
}

template <typename T>
PolytopeFace<T>* Polytope<T>::AddNewFace(PolytopeEdge<T>* e1,
                                         PolytopeEdge<T>* e2,
                                         PolytopeEdge<T>* e3) {
  // Do not continue if input invalid
  if (e1 == nullptr || e2 == nullptr || e3 == nullptr) return nullptr;

  // Do malloc
  if (face_malloc_list_.malloc_next == nullptr) return nullptr;
  auto* face = static_cast<PolytopeFace<T>*>(face_malloc_list_.malloc_next);
  face_malloc_list_.malloc_next = face->malloc_next;

  // Init the face
  face->id = next_face_id_;
  next_face_id_ += 1;
  face->type = PolytopeElementType::Face;
  face->edges_of_face[0] = e1;
  face->edges_of_face[1] = e2;
  face->edges_of_face[2] = e3;

  // The min distance
  const Vector3<T> a = e1->vertices_of_edge[0]->vertex_location;
  const Vector3<T> b = e1->vertices_of_edge[1]->vertex_location;
  Vector3<T> c;
  auto e = face->edges_of_face[1];
  if (e->vertices_of_edge[0] != face->edges_of_face[0]->vertices_of_edge[0] &&
      e->vertices_of_edge[0] != face->edges_of_face[0]->vertices_of_edge[1]) {
    c = e->vertices_of_edge[0]->vertex_location;
  } else {
    c = e->vertices_of_edge[1]->vertex_location;
  }
  pointToTriangleSquaredDistance<T>(Vector3<T>::Zero(), a, b, c,
                                    face->min_distance);

  // add to face list
  addElementToList(face_list_, face);

  // update the element in edges
  for (auto i = 0; i < 3; i++) {
    PolytopeEdge<T>* edge_i = face->edges_of_face[i];
    if (edge_i->faces_of_edge[0] == nullptr) {
      edge_i->faces_of_edge[0] = face;
    } else {
      if (edge_i->faces_of_edge[1] != nullptr) {
        // This might due to wrong connection at degeneration
        // std::cout << "Wrong edge" << std::endl;
        return nullptr;
      }

      assert(edge_i->faces_of_edge[1] == nullptr);
      edge_i->faces_of_edge[1] = face;
    }
  }

  return face;
}

template <typename T>
PolytopeElementBase* Polytope<T>::ComputeMinDistanceToOrigin(
    bool exclude_vertex) {
  // Empty polytope, do nothing
  if (isEmptyList(vertex_list_)) return nullptr;

  // Update the min-distance if necessary
  T min_distance_squared = std::numeric_limits<T>::infinity();
  PolytopeElementBase* min_distance_element{nullptr};
  auto visit_element_functor =
      [&exclude_vertex, &min_distance_element,
       &min_distance_squared](PolytopeElementBase* element) -> void {
    // Now it holds a distance
    switch (element->type) {
      case PolytopeElementType::Vertex: {
        auto* vertex = static_cast<PolytopeVertex<T>*>(element);
        assert(vertex->min_distance.witness_in_simplex);
        if (vertex->min_distance.min_distance_square < min_distance_squared) {
          min_distance_element = vertex;
          min_distance_squared = vertex->min_distance.min_distance_square;
        }
        return;
      }
      case PolytopeElementType::Edge: {
        auto* edge = static_cast<PolytopeEdge<T>*>(element);
        if ((!exclude_vertex) && (!edge->min_distance.witness_in_simplex))
          return;
        if (edge->min_distance.min_distance_square < min_distance_squared) {
          min_distance_element = edge;
          min_distance_squared = edge->min_distance.min_distance_square;
        }
        return;
      }
      default: {
        assert(element->type == PolytopeElementType::Face);
        auto* face = static_cast<PolytopeFace<T>*>(element);
        if (!face->min_distance.witness_in_simplex) return;
        if (face->min_distance.min_distance_square < min_distance_squared) {
          min_distance_element = face;
          min_distance_squared = face->min_distance.min_distance_square;
        }
        return;
      }
    }
  };

  // Go
  if (!exclude_vertex) visitList(vertex_list_, visit_element_functor);
  visitList(edge_list_, visit_element_functor);
  visitList(face_list_, visit_element_functor);
  return min_distance_element;
}

template <typename T>
bool Polytope<T>::ComputeFaceNormalPointingOutward(const PolytopeFace<T>* face,
                                                   Vector3<T>& normal,
                                                   T* area) const {
  Vector3<T> a, b, c;
  auto vertex_ok = face->GetVertices(a, b, c);
  // This one should never fail
  assert(vertex_ok);
  (void)(vertex_ok);  // disable warning

  const Vector3<T> e1 = a - b;
  const Vector3<T> e2 = b - c;
  const Vector3<T> e1_cross_e2 = e1.cross(e2);
  const T norm = e1_cross_e2.norm();
  if (area) *area = T(0.5) * norm;
  Vector3<T> original_direction;
  if (norm <= 0.0) {
    // Degenerate case, this may fail.
    // Use one if its neighbor's face?
    // Or just inform external
    return false;
  } else {
    original_direction = e1_cross_e2 / norm;
  }

  // Case 1: we can directly determine from this face
  const Vector3<T>& o_to_a = a;
  const T o_to_a_dot_n = o_to_a.dot(original_direction);
  constexpr T normal_tolerance = 1e-4;
  if (std::abs(o_to_a_dot_n) >= normal_tolerance) {
    if (o_to_a_dot_n > 0)
      normal = original_direction;
    else
      normal = -original_direction;
    return true;
  }

  // Case 2: we need to rely on other vertices as o is on plane abc
  T max_positive_dot_product = 0.0;
  T min_negative_dot_product = 0.0;
  auto update_dot_product_functor =
      [&original_direction, &max_positive_dot_product,
       &min_negative_dot_product](const PolytopeElementBase* element) -> void {
    assert(element->type == PolytopeElementType::Vertex);
    const auto* vertex = static_cast<const PolytopeVertex<T>*>(element);
    const auto dot_with_original =
        vertex->vertex_location.dot(original_direction);
    if (dot_with_original > max_positive_dot_product)
      max_positive_dot_product = dot_with_original;
    if (dot_with_original < min_negative_dot_product)
      min_negative_dot_product = dot_with_original;
  };
  auto visit_ok =
      VisitElementList(PolytopeElementType::Vertex, update_dot_product_functor);
  assert(visit_ok);  // Again, this never fails
  (void)(visit_ok);  // disable warning in Release

  // Choose the larger one
  if (std::abs(max_positive_dot_product) > std::abs(min_negative_dot_product)) {
    normal = -original_direction;
  } else {
    normal = original_direction;
  }
  return true;
}

}  // namespace cvx_collide
}  // namespace fcl