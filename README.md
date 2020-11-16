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

https://gitlab.enpc.fr/cerea/munich


Installation
============

Buliding
--------

```
$ cd processing/photochemistry
$ scons
```

Run
---

```
$ munich munich.cfg
```

Input files
-----------

A Python preprocessing tool is available to generate the input data and
visualize then by contact with Youngseob Kim (youngseob.kim@enpc.fr).


## --- intersection.dat

This file contains the informations concerning intersections.

The format of the input data is as follows:

<intersection id>;<longitude>;<latitude>;<number of streets> (which are connected to the intersection);<a series of street id>;


## --- street.dat

This file contains the input data of each street segment.

The format of the input data is as follows:

<street id>;<intersection 1>;<intersection 2>;<street length>;<averaged builing height>;<street width>
