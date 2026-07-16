# Introduction

Notes on creating this python package, the difficulties, and things learned along the way.

## Overview
### goal
Model the dissolution of star clusters and the consequent production of stellar streams to facilitate research in galactic archeology. 

### Constraints
This problem has a large parameter space to be explored and has a large data volume. We use the restricted three body problem instead of full n-body for faster computations, at the cost of simplifying the internal dynamics of the system. The valididty of this assumption depends on the density of the globular cluster system. 

The parameter space of the problem is constrained by the number of streams we wish to model, and the various parameters requied to model a given stream. For starters, the Milky Way has about ~170 globular clusters and about ~120 streams (as of 2026), with an overlap of about ~20. Each globular cluster system has observational uncertainties in both its kinematics and structural properies. 

- **Observational Uncertainties**
    - kinematics
        - proper motions
        - parallaxes
        - radial velocities 
    - structural 
        - mass
        - size (?)
            - Note. The Baumgardt catalog does not provide uncertainties for the sizes. Briefly, they take observables of each cluster, such as the light function and velocity disperson and find best fitting n-body models to match them. They are described with king models, which is a 3 parameter model: the total mass, the core radius, and the tidal radius. Only the mass has uncertainty. 
- **Milky Way Potential**
    - the mass distribution of the Milky Way, and it's time evolution, is not known. Improving our knowledge of the Milky Way's gravitational field is a goal of the field. So, need to investigate various models. We can decompose the Milky Way as: 
        - __Static models__: Instead of considering the Milky Way as a collection of stars, it can be approcimated as a continuous density distribution and it's corresponding gravitational potential, which are linked through Poisson's equation. Also, we can consider this to be time-independent. Although the code should be extendable to analyze varying said components. There are many Milky Ways
            - Pouliasis et al. (2017)
            - Ibata et al. (2024)
            - McMilian et al. (2017)
            - Bovy et al. (2014)
        - __time varying components__: we may wish to include, as a higher order effect, time-vary components
            - Milky Way satellites: the LMC, SMC, Sagittarius, etc. 
            - Dark Matter Subhalos
            - globular clusters
            - spiral arms
            - the galactic bar


The number of simulations is thus vast. We can have:

$$\rm{Data~Volume} \approx N_{\rm{MWmodels}}\times N_{\rm{GCs}}\times N_{\rm{MC}}\times N_{\rm{p}}\times 6\times 4~\rm{bytes}$$

In our first publication, we explored 3 MW models, for 150 globular clusters, with 50 Monte-Carlo samplings of the uncertainties, and modelled each with 100,000 star particles. Is we save single precision data (4 bytes), we get a volume of approximately 430 gigabytes. 


### How
Solve hamilton's equations numerically using the restricted three body problem. 

Write a code in Fortran to do so. Build a python package with `f2py` so the computations can be handled with fortran and the i/o and analyses can be handled in python. 

I wanted as little dependencies as possible. This required writing the poisson solver for the gravitational field by hand. This is okay when the potentials were given analytically, or with density distributions that used the spherical harmonics only. However, as I extended to very flat systems, I needed to find the potential from the density distribution only. I implemented a solution of Poisson's equations in cylindrical coordinates using a hankle transform. This gives periodic boundary conditions and wrong solutions when you are far from the systems center. Then, I changed gears to use `agama` as a back end for the gravitational solver. 
