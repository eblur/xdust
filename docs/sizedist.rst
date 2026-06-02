xdust.graindist.sizedist
============================

Abstract Class: *Sizedist*
--------------------------

Abstract class *Sizedist* must contain the following:

**Atributes**

- ``dtype`` : a string description

- ``a`` : float or np.array [grain size, micron] or astropy Quantity

**Methods**

- ``ndens(md, rho, shape=xdust.graindist.shape.Sphere)`` 
   Returns number density of dust grains [e.g. cm^-2 um^-1]

   ``md`` -- float [dust mass column density in g cm^-2]

   ``rho`` -- float [dust mass density in g cm^-3]

   ``shape`` -- xdust.graindist.shape object (abstract class)

- ``mdens(md, rho, shape=xdust.graindist.shape.Sphere)`` 
   Returns the mass density of dust grains [e.g. in g cm^-2 um^-1]

   ``md`` -- float [dust mass column density in g cm^-2]

   ``rho`` -- float [dust mass density in g cm^-3]

   ``shape`` -- xdust.graindist.shape object (abstract class)

Grain
-----
.. autoclass:: xdust.graindist.sizedist.Grain

Powerlaw
--------
.. autoclass:: xdust.graindist.sizedist.Powerlaw

ExpCutoff
---------
.. autoclass:: xdust.graindist.sizedist.ExpCutoff

WD01
----
.. autoclass:: xdust.graindist.sizedist.WD01

Astrodust
---------
.. autoclass:: xdust.graindist.sizedist.Astrodust


