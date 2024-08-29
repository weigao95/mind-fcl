## Several Design Decisions

##### 1. Separation distance is no longer supported

Separation distance is the minimum distance between two *non-overlapping* shapes. In robot motion planning, separation distance can be considered as a "safety margin" or an implicit "padding" of the original geometry. For a pair of convex shapes, separation distance can be computed using the GJK algorithm. This is an convex optimization problem and can be solved rather efficently.

On the other hand, separation distance between non-convex shapes (bvh/octree/heightmap) can be computationally expensive. During distance computation, the node pair in bvh traversal tree can not be pruned using binary overlapping predicate. A good distance upperbound is required to prune the node pair, which is not easy for tight bounding volumes (such as OBB). As a result, we can only use AABB pair to compute the distance upper-bound.

Thus, we decide to drop the support of separation distance in the `fcl::distance` interface, while the penetration distance can still be computed using the `fcl::collide` interface. In practice, we explicitly apply the "padding" to the geometries (e.g., enlarge the shape with the given "safety margin").

##### 2. Non-convex shapes

Robot links, environements and gripper tools can be offline convexified by computing the convex hull or convex decomposition. However, point cloud or objects without CAD models necessitate non-convex shapes such as octree, heightmap and/or general BVH.

Additionally, due to the complexity of these non-convex shapes (usually we need to handle point clouds with 10^6-10^8 points), their collision detection can account for a majority of the computation time (e.g., more than 90%) during motion planning. To satisfied this requirement, we re-write the traversal collision detection algorithms for bvh/octree/heightmap for conduct a variety of performance optimization.

##### 3. Continuous collision detection

The library implements translation continous collision detection between two convex and non-convex shapes. While most movement in robotic manipulation involves rotation (and is not translational), there are two importance special cases:

1. The approaching movement of robot end-effector immediately before picking up an object
2. The retraction movement of robot end-effector immediately after picking up an object

These two types of movement are typically short linear movement of the end-effector. Thus, we can use translational collision detection during these movements for gripper tool and picked objects (vs. point cloud and perceived objects). As mentioned above, these collision detection typcially involves very complex non-convex shapes and account for a majority of computation time. Thus, integrating the continuous collision detection can subsentailly improve the performance of the overall manipulation planning.
