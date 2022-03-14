# About

A parallel approach, based on a set of successful algorithms, for real-time simulation of the interaction of fluids with rigid bodies. The method is based on the Smoothed Particle Hydrodynamics (SPH) and developed within a heterogeneous multi-core CPU and GPU architecture. 

Our discretisation of polygonal model into a set of particles is made by using a modified version of depth peeling, which was originally used for rendering transparent polygons. The technique has also been applied to various other operations related to collision detection and polygonal discretisation.

In order to compile the project, it is necessary to use the [jrfxgl](https://github.com/josericardojr/jrfxgl) engine.

# Challenges
The main challenge of this project was distribute tasks among CPU and GPU and each of these processors have its own architecture.
Also, in order to minimize the data communication between GPU and CPU, which tends to reduce the overall speed, a new data structure was created.
Interaction between fluid and rigid body is done by covering the 3D model surface with a collection of spheres using a technique called Depth Peeling.

# Features
* CPU and GPU usage.
* Two way interaction between fluids and rigid body.
* Depth Peeling technique for discretizing rigid bodies.


# Documentation
* [A heterogeneous system based on GPU and multi-core CPU for real-time fluid and rigid body simulation](http://dx.doi.org/10.1080/10618562.2012.683789)
* [Neighborhood grid: A novel data structure for fluids animation with GPU computing](http://dx.doi.org/10.1016/j.jpdc.2014.10.009)
* [An architecture for real time fluid simulation using multiple GPUs](http://www.sbgames.org/sbgames2012/proceedings/papers/computacao/comp-full_12.pdf)
* [Two-Way Real Time Fluid Simulation Using a Heterogeneous Multicore CPU and GPU Architecture](http://dx.doi.org/10.1109/pads.2011.5936750)
* [Fluid simulation with two-way interaction rigid body using a heterogeneous GPU and CPU environment](http://sbgames.org/papers/sbgames10/computing/full/full19.pdf)


# Videos
<a href="http://www.youtube.com/watch?feature=player_embedded&v=cTA0uPKqPpY
" target="_blank"><img src="http://img.youtube.com/vi/cTA0uPKqPpY/0.jpg" 
alt="IMAGE ALT TEXT HERE" width="480" height="360" border="10" /></a>

<a href="http://www.youtube.com/watch?feature=player_embedded&v=iyM1bkVSc6U
" target="_blank"><img src="http://img.youtube.com/vi/iyM1bkVSc6U/0.jpg" 
alt="IMAGE ALT TEXT HERE" width="480" height="360" border="10" /></a>
