xdust.grainpop
==============

The ``SingleGrainPop`` object is used to model a population of dust grains 
based on three distinct properties:

* the grain size distribution (see ``xdust.graindist.sizedist``)
* composition type specifying the optical constants (see ``xdust.graindist.composition``)
* the scattering physics model (see ``xdust.scatteringmodel``)

Furthermore, the **GrainPop** object contains a list of **SingleGrainPop** objects, 
which can be called upon using keys, like a dictionary.  Keys have to be specified 
when the **GrainPop** object is initialized.

SingleGrainPop
--------------

.. autoclass:: xdust.grainpop.SingleGrainPop

GrainPop
--------

.. autoclass:: xdust.grainpop.GrainPop

Helper functions
----------------

.. autofunction:: xdust.grainpop.make_MRN
.. autofunction:: xdust.grainpop.make_MRN_RGDrude
