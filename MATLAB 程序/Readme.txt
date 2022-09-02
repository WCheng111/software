THANK YOU FOR USING MPUI - A GRAPHICAL USER INTERFACE TO EXPLORE, STUDY AND ANALYSE THE GEOMETRIC PROPERTIES OF MOIRÉ PATTERNS IN THE CONTEXT OF 2D MATERIALS AND ATOMIC LATTICES.

*******************************************************************************
MPUI version 1.0.1 -- 14 November 2019
*******************************************************************************

AUTHORS:
Maxime Le Ster, Tobias Märkl, Simon Brown
The MacDiarmid Institute for Advanced Materials and Nanotechnology &
School of Physical and Chemical Sciences
University of Canterbury, Private Bag 4800, Christchurch 8140, New Zealand

DISCLAIMER:
While you are welcome to change the code according to your needs and redistribute it under the terms of the license, this is our final version of the program. We appreciate your interest but we cannot offer any support, and we will not develop the program any further.


********************************************************************************
By using this program you agree to cite the following research article if this
program contributes to a publication which you are (co-)authoring:

M. Le Ster, T. Märkl & S. Brown, "Moiré Patterns: A Simple Analytical Model"
2019, 2D Materials 7, 011005
DOI: https://dx.doi.org/10.1088/2053-1583/ab5470
********************************************************************************


CONTENTS OF THIS DOCUMENT:
1) A list of files
2) An introduction to the program
3) Instructions for using the UI
4) Known issues

********************************************************************************

1) LIST OF FILES
The MPUI zip archive should contain:
* this Readme document
* a matlab function file MPUI.m
* a matlab function file twistAngleDependence.m
* a License.txt file



2) INTRODUCTION
When two regular, periodic lattices are overlaid usually a third, new spatial periodicity arises as a consequence of the mismatch between the two layers. This is a manifestation of so-called moiré patterns (MPs), and in particular it is an interesting phenomenon in the nanoscale world of 2-dimensional materials on surfaces. 

The program allows you to do a couple of things regarding these MPs:
* visualise the MPs emerging from the mismatch of a pair of atomic lattices with arbitrary lattice symmetry and lattice constants
* examine the reciprocal lattices of both layers and their interference terms
* find the best candidates (in terms of a moiré period and orientation) that match the visible patterns
* study the effect of the relative twist angle between the two layers on the individual interference terms (i.e. their period and angle) 

For the last aspect in particular we implement a set of formulas (equations 4 & 7 from the research article cited above) that allows us to calculate the characteristic wavelengths and angles of moiré fringes for a pair of atomic lattices, called the under- and overlayer (UL/OL). 

	
	
3) INSTRUCTIONS
Extract the .m files to a convenient location and navigate there in Matlab, then run it by calling "MPUI", or open the MPUI.m in Matlab's editor and run it from there. The twistAngleDependence.m file must be in the same folder as MPUI.m, or its location must be added to the Matlab "path" variable.

After starting, a figure will be created that contains several plots and UI elements, grouped together in certain blocks:

	A) Visualization (big blue box)
	This is the biggest block of the UI and mainly contains two plots that show the two lattices and the emerging MPs on the left hand side, as well as the reciprocal space lattices on the right hand side. In both plots, the underlayer is shown in black, the overlayer in red. You can zoom and pan in all three panels and there are two reset buttons return the main views to a default. The "Highlight ..." button is explained below in section E).
	
	A small inset in the "Real Space" plot shows a close-up of the origin, independently of the zoom level of the large panel. It will return to a default view automatically when you make changes to the geometry.
	
	The blue dots in the reciprocal space plot are difference vectors of pairs of reciprocal lattice vectors of the under- and overlayer (one each). MPs can be represented as a superposition of all these difference vectors but it is often sufficient to consider only the ones that originate from the few lowest order reciprocal lattice points. For that reason a cut-off in k-space can be imposed which is represented in the plot by the shaded, pale green circle.
	
	
	B) Input Geometry
	This is the main input block. Here you can enter the geometrical details (R1/2 for the lattice parameters, omega for the angle between the unit cell basis vectors) of the two lattices and their relative twist angle theta. In order to more closely resemble experimental data, a global rotation angle gamma can be set as well which is merely a rotation of the coordinate system and does not affect the characteristic MP quantities.
	
	A slider is used to change the twist angle in small increments by clicking on the arrows/in the space between bar and arrow (0.1°/1°). This allows to access the range between 0° and 360° but a desired, specific value can be directly set in the input field above the slider.
	
	Furthermore there is a dropdown menu from which pre-defined examples can be loaded. These illustrate either widely-studied cases (such as twisted bilayer graphene) or cases with different symmetries such as a rectangular 2 monolayer-thin bismuth phase on a substrate with hexagonal symmetry (alpha-Bi on MoS2).
	
	
	C) Display options
	You can alter the appearance of the Real Space plot in this block in a couple of ways. The number of visible unit cells for each layer can be specified, this is useful when the two lattices have strongly different lattice constants. Note that (2n+1)^2 lattice points will be shown if the entered number is n. We limit these numbers to values between 4 and 100; if you really want to display more or less unit cells, this can be changed in the function "inputNumberOfUnitCells".
	Furthermore, you can choose whether and where to show a second atom of the lattices which can be desired in the case of a two-atom basis such as graphene. The coordinates here are relative values with respect to the size of a unit cell, i.e. allowed values are between 0 and 1.
	
	
	D) Found moiré patterns
	This block contains information about those MP contributions with the largest MP period (= the shortest k-vectors) for the specific geometry and twist angle defined in the block above. The period and angle of up to three sets of MP fringes (M1-M3) will be displayed in the table and the corresponding points are shown as dark blue dots in the reciprocal space plot.
	The input field and slider control the above-mentioned cut-off in reciprocal space. Per default the cut-off cannot be less than the either of the shorter reciprocal lattice vectors of under- and overlayer, so that there is always at least one difference vector. Its upper limit is set to ~4 times the largest of the reciprocal lattice basis vectors. This can be changed in the code if desired. (change the default values of ksp.ULnUnitCells & ksp.OLnUnitCells.)
	
	
	E) Twist angle dependence
	While the reciprocal space plot in A) show all difference vectors compatible with the cut-off, these last two plots show the twist angle dependence of only ONE specific difference vector. This vector is determined by the four indices m,n,p,q which are used to select reciprocal lattice vectors of the underlayer (m,n) and overlayer (p,q) - one each. All three points can be highlighted in the reciprocal lattice plot using the "Highlight M,N,P,Q" button. The range of twist angles theta is between +/-180° with steps of 0.1°.
	These two plots use the external function twistAngleDependence.m which can be used independently of this UI. Refer to its code for usage instructions.



4) KNOWN ISSUES
We are currently aware of a few (performance) issues:

* When using the twist angle slider for an extended time the UI can become unresponsive for a couple of seconds (usually 5-10s but we have experienced up to ~20 seconds delay in extreme cases). It will return to normal operation mode after that.
* Particularly on slower machines, laptops in energy-saving mode, or when hardware-accelerated rendering is not available, the output may be slow to update after changing the twist angle theta.
* The size and position of all UI elements is hard-coded as a number of pixels. If you use a high-DPI screen and have changed the operating system setting, the UI might appear improperly scaled and/or some elements may not be visible. In that case, please change the commented lines in the code at the beginning of the "addUIcomponents" function (h.fig = ...). That will allow resizing of the figure window and you should then be able to adjust it to a size that you can use.