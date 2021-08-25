# image-meshing

*Contour-based meshing of the image domain developed for CINEMAX summer school.*

The module `meshing.py` contains the functionality for contour-based meshing of the image. Contour detection uses marching squares implementation from the [scikit image module](https://scikit-image.org/). For computing the conforming constrained Delaunay triangulation we rely on the functionality from [triangle package](https://rufat.be/triangle/), which in turn wraps around [Jonathan Richard Shewchuk’s mesh generator](http://www.cs.cmu.edu/~quake/triangle.html). 

For 3-phase segmentation check [bone meshing notebook](bone_meshing.ipynb), which is also a most documented notebook.

For 2-phase segmentation check [bundles meshing notebook](bundles_meshing.ipynb), or any of the other provided notebooks




 
