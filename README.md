What is MUNICH
==============

The Model of Urban Network of Intersecting Canyons and Highways (MUNICH) is
used to simulate subgrid concentrations in the urban canopy represented by the
street network. MUNICH is designed for the coupling to Polair3D
chemical-transport model to form a Street-in-Grid model (SinG) on the
Polyphemus air-quality modeling platform (see http://cerea.enpc.fr/polyphemus/).
MUNICH is distributed under the GNU General Public License.

For more information on MUNICH, see http://cerea.enpc.fr/munich/

Download
========

From Zenodo

https://doi.org/10.5281/zenodo.4168985

or from the git public repository

https://github.com/cerea-lab/munich


Installation
============

ssh-aerosol is used to model gas-phase and particle-phase species.

photochemistry is used to model only gas-phase species.

Building ssh-aerosol
--------------------

```
$ cd processing/ssh-aerosol
$ compile
```

Run ssh-aerosol
---------------

```
$ munich-ssh munich.cfg
```


Building photochemistry
-----------------------

```
$ cd processing/photochemistry
$ scons
```

Run photochemistry
------------------

```
$ munich munich.cfg
```


Input files
-----------

Input data for a test case are available at http://cerea.enpc.fr/munich/


Users's Guide
=============

http://cerea.enpc.fr/munich/doc/munich-guide-v2.pdf


Help
====

Send e-mail to munich-help@liste.enpc.fr
