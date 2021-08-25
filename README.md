# image-meshing

*Contour-based meshing of the image domain developed for CINEMAX summer school.*

The module `meshing.py` contains the functionality for contour-based meshing of the image. Contour detection uses marching squares implementation from the [scikit image module](https://scikit-image.org/). For computing the conforming constrained Delaunay triangulation we rely on the functionality from [triangle package](https://rufat.be/triangle/), which in turn wraps around [Jonathan Richard Shewchuk’s mesh generator](http://www.cs.cmu.edu/~quake/triangle.html). 

For 2-phase segmentation check [bundles meshing motebook](bundles_meshing.ipynb).



 
