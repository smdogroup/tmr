[![Build and unit tests](https://github.com/smdogroup/tmr/actions/workflows/unit_tests.yml/badge.svg)](https://github.com/smdogroup/tmr/actions/workflows/unit_tests.yml)

# TMR #

TMR is a mesh generation and refinement tool for creating large-scale quadrilateral and octree meshes in parallel. TMR creates a coarse quadrilateral or hexahedral mesh and builds refined quadtree or octree meshes on each coarse element. To ensure compatibility, TMR computes the necessary compatibility relationships between adjacent quadrilateral or hexahedral elements of different refinement levels.

Online documentation and examples is located at [https://smdogroup.github.io/tmr/](https://smdogroup.github.io/tmr/)

# Please cite us #

The key contributions and development of TMR are described in the following paper: 
Ting Wei Chin, Mark K. Leader, Graeme J. Kennedy, A scalable framework for large-scale 3D multimaterial topology optimization with octree-based mesh adaptation, Advances in Engineering Software, Volume 135, 2019.

```
@article{Chin:2019,
         title = {A scalable framework for large-scale 3D multimaterial topology optimization with octree-based mesh adaptation},
         journal = {Advances in Engineering Software},
         volume = {135},
         year = {2019},
         doi = {10.1016/j.advengsoft.2019.05.004},
         author = {Ting Wei Chin and Mark K. Leader and Graeme J. Kennedy}}
```

TMR is open source with the Apache License.
