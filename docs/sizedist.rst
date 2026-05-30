xdust.graindist.sizedist
============================

Abstract Class: *Sizedist*
--------------------------

Abstract class *Sizedist* must contain attributes:

- `dtype` : a string description

- `a` : an array

- `ndens` (md, rho, shape) returns number density of dust grains [e.g. cm^-2 um^-1]

- `mdens` (md, rho, shape) returns mass density of dust grains [e.g. g cm^-2 um^-1]

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


