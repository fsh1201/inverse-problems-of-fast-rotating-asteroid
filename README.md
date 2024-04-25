# Description
The program (written in Fortran by Mikko Kaasalainen and converted to C by Josef Ďurech)
computes the shape+spin+scattering model that gives the best fit to the input lightcurves
(calibrated, uncalibrated, or sparse). Scattering law+shape representation are simple and
robustly converging. The shape representation this procedure obtains is the Gaussian image
of a convex polyhedron, i.e., the areas of the facets (with fixed outward normals). The vertices
are directly solved by the Minkowski procedure minkowski.

This program is a **modified version** that takes the **exposure time** into account for those **high 
exposure-period-ratio** situations.

# Cite this work
    @ARTICLE{2024MNRAS.528.3523F,
           author = {{Feng}, Shuai and {Hu}, Shaoming and {Chen}, Xu and {Li}, Yang and {Du}, Junju and {Yang}, Zhitao and {Cao}, Hai and {Gan}, Qingbo and {Liu}, Shuqi and {Jiang}, Yuchen},
            title = "{Spin state and convex shape inversion from light curves of fast-rotating asteroids}",
          journal = {\mnras},
         keywords = {methods: data analysis, methods: numerical, minor planets, asteroids: general},
             year = 2024,
            month = feb,
           volume = {528},
           number = {2},
            pages = {3523-3530},
              doi = {10.1093/mnras/stae250},
           adsurl = {https://ui.adsabs.harvard.edu/abs/2024MNRAS.528.3523F},
          adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }



# Compile
    $ make  
g++, gcc, and gfortran are required

# Usage
## Inversion of spin state
    $ ./spin_axis lcs input_period_scan out_periods
    
## Inversion of shapes
    $ cat lcs | ./convexinv -s -o out_areas -p out_par input_convexinv out_lcs
    $ cat out_areas | ./minkowski | ./standardtri | ./shape2obj > model.obj
    
## Input lightcurves (lcs)
The input file contains lightcurve data and the corresponding geometry, it is read from the
standard input. 

The first line gives the total number of lightcurves, then the individual
lightcurves follow in ‘blocks’. 

Each lightcurve starts with the number of points and 0/1 code
for a relative (0) or calibrated (1) lightcurve. 

Then there are lines with the epoch in JD
(light-time corrected!), the brightness in intensity units (reduced to unit distances from the
Earth and the Sun when calibrated), the ecliptic astrocentric cartesian coordinates $x, y, z$ of
the Sun and of the Earth in $AU$, and the exposure time in $day$.  

For a slowly moving main belt asteroid, the coordinate vectors can be approximated as
to be constant for a single-night lightcurve.

An [example file](https://github.com/fsh1201/inverse-problems-of-fast-rotating-asteroid/blob/main/lc.txt)

## input_period_scan
**Initial period** - The first line of the input period scan gives the initial period, the final period, and the
coefficient $p$ of the period step. The interval is scanned with the step $p\Delta P$ (always set $p < 1$,
a recommended value is $p = 0.8$). 

The initial values are followed by 0 or 1 depending on whether those are **fixed** (0) or
**free** (1). 

**Convexity regularization weight** – this is often needed to keep the shape formally convex.
A typical value is $0.1$, but it may often have to be orders of magnitude larger or smaller.

There is a dark facet area that makes the whole set of facets convex. Always try to put
it below 1% of the total area by increasing the convexity regularization parameter. If you
cannot get a good fit with a small dark area and if the fit significantly improves when the
dark facet is large, it means that there is likely to be an albedo variation over the surface.

**Laplace series expansion** – degree l, and order m. It affects the number of shape parameters –
if l = m, the number of shape parameters is $(l + 1)^2$
(and $(l + 1)^2 − 1$ for relative lightcurves).
A good choice is $l = m = 6$.

**Resolution** – number n of triangulation rows per octant (typically $8–10$). The number of
surface areas of the Gaussian image is $8n^2$.

**Light scattering parameters** – amplitude a, width d, slope k, and Lambert’s coefficient c. Initial
values are followed by 0/1 code for fixed (0) or free (1) parameters. If you use only relative lightcurves, 
set them fixed as $a = 0.5, d = 0.1, k = −0.5, c=0.1$.

**Iteration stop condition** – if it is an integral number higher than one, then it is the number
of iteration steps in the Levenberg-Marquardt loop. If it is lower than one, then it is the
smallest difference in rms deviation between two subsequent steps – when the steps have
smaller difference, the iteration loop is stopped.

An [example](https://github.com/fsh1201/inverse-problems-of-fast-rotating-asteroid/blob/main/input_period_scan)

## out_periods
From the 1st column to 7th column: it's $P$ (h), $\lambda$ (deg), $\beta$ (deg), RMS residual, 
residual, iteration times, dark area %

## out_areas
The first line gives the number of facets, then follow facet areas and outward unit normal
$x, y, z$ coordinates. This is used as an input for minkowski procedure that creates a 3D model
(see below). Note that the program writes the size and normal coordinates of a dark facet
at the end, so the number of facets is $8n^2 + 1$. The dark facet should be very small (control
with convexity regularization parameter). It is needed to make the collection convex, but it
does not contribute photometrically.

## out_par
This file contains the solution for the spin vector direction, period, and scattering parameters
in the format:  
$λ$ $β$ $P$  
$t_0$ $φ_0$  
$a$ $d$ $k$  
$c$  
If some of the parameters are fixed ($t_0$ and $φ_0$ are always fixed), they have the same values
as in the **input_convexinv** file.

## input_convexinv
**Initial spin** – asteroid’s initial ecliptic pole coordinates $λ, β$ (deg), and the rotation period $P$
(hours). The initial values are followed by 0 or 1 depending on whether those are fixed (0) or
free (1). The initial value could be taken from the best result of spin state inversion.

**Zero time** $t_0$ (JD), **initial rotation angle** $φ_0$ (deg) – they are needed for the transformation
between vectors rast in the asteroid co-rotating coordinate frame and vectors recl in the
ecliptic coordinate frame.
If $t_0 ≤ 0$ in the input file, then it is set to the lowest JD epoch of the data set. So the
basic choice is $t_0 = 0$ and $φ_0 = 0$.

Other input parameters have the same meaning as in the
**input_period_scan** file.

An [example](https://github.com/fsh1201/inverse-problems-of-fast-rotating-asteroid/blob/main/input_convexinv)

## out_lcs
This file has the same data map of input lightcurves (**lcs**), but without exposure time.

## a test
    $ ./spin_axis ./lc.txt ./input_period_scan ./out_periods
    $ cat ./lc.txt | ./convexinv -s -o ./out_areas.txt -p ./out_par.txt ./input_convexinv ./out_lcs.txt
    $ cat ./out_areas.txt | ./minkowski | ./standardtri | ./shape2obj > ./model.obj

# More
This README is pretty much a copy of [doc.pdf](https://github.com/fsh1201/inverse-problems-of-fast-rotating-asteroid/blob/main/doc.pdf). Read it for more details. 
Also see [DAMIT website](https://astro.troja.mff.cuni.cz/projects/damit/)
