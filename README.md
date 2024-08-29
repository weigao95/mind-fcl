## mind-fcl -- An extension of the Flexible Collision Library

[FCL](https://github.com/flexible-collision-library/fcl) was forked and a large part of the code has been rewritten. This library is the backbone for the manipulation planner of [Mech-Mind Robotics](https://www.mech-mind.com/) and has been used in thousands of robot workstations that we have deployed.

### New Features

Compared to the original [FCL](https://github.com/flexible-collision-library/fcl) library, the main new features are:
- a new incremental penetration distance algorithm described in [this paper](https://arxiv.org/abs/2304.07357)
- a new implementation of the convex collision algorithms (GJK/EPA/MPR, in [this subdir](./include/fcl/cvx_collide/)). The new EPA is 2x faster than [libccd](https://github.com/danfis/libccd) in the original [FCL](https://github.com/flexible-collision-library/fcl) library
- optional header-only integration, `Eigen` is the only dependency ([libccd](https://github.com/danfis/libccd) is no-longer required)
- the support of new geometry shapes such as height map and tetrahedron mesh
- re-write and performance optimization of traversal collision detection algorithms for bvh/octree/heightmap
- better multi-threading support with a new broadphase AABB-tree and read-only narrowphase interface
- efficient translational continuous collision detection
- various bug-fix and performance optimization

### Installation

Please make sure Eigen is installed on your system and can be found in the include path. After that, you can simply copy the `/include` directory (which contains `fcl/fcl.h` and others) to your own project. In this way, there is no binary library to link to.

Alternatively if you want to build a library, CMakeLists.txt can be used to generate makefiles in Linux or Visual studio projects in windows. In command line, run

``` cmake
mkdir build
cd build
cmake ..
```

Next, in linux, use make to compile the code. In windows, there will generate a visual studio project and then you can compile the code.

### Usage

The narrowphase interface `fcl::collide` is almost the same as the original [FCL](https://github.com/flexible-collision-library/fcl) library. In some cases, it can be used as a drop-in replacement. Please refer to the document in the origial [FCL](https://github.com/flexible-collision-library/fcl).

For the broadphase, a new AABB tree inspired by [JoltPhysics](https://github.com/jrouwe/JoltPhysics) is implemented for better concurrency support. The [test code](./test/broadphase) might be a good starting point.

### Design Decisions

We have learned a lot regarding collision detection from our experience of deploying of robot manipulators. A summary of several design decisions behind this library is [here](./Design.md)

### Contact

If you find a bug or have a new feature request, please send me an email.

