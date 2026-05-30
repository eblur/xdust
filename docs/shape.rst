xdust.graindist.shape
=====================

At this time, only spherical dust grains are supported. 
For non-spherical dust grains, the grain size distribution 
can be specified using the effective grain radius $a_{\rm eff}$, defined as
$$ V_g = \frac{4}{3} \pi a_{\rm eff}^3 $$
where $V_g$ is the volume of the grain.

Abstract Class: *Shape*
--------------------------

Abstract class *Shape* describes the shape of the dust grains.
It must contain the following attributes:

- `shape` (a string)

- `vol` (a function giving the volume of the grain [cm^3])

- `cgeo` (a function giving the geometric area of the grain [cm^2])

Sphere
------
.. autoclass:: xdust.graindist.shape.Sphere
