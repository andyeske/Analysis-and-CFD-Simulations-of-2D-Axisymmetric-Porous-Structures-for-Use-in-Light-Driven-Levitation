# Analysis and CFD Simulations of 2D Axisymmetric Porous Structures for Use in Light Driven Levitation

In this repository, we present the various MATLAB functions and scripts utilized in our work, alongside an explanation of their use, inputs and outputs.

## MATLAB Code

```
calc_F
```

**Type:** Function

**Description:** The purpose of this code is to provide a framework for calculating the flow-through velocity (vft) and outlet velocity (vout) of an arbitrary 3D geometry.

**Inputs:** altitude, geom_param, chan_param

* <ins>altitude:</ins> this is a vector containing the specific altitudes (in km) at which the vft will be calculated at (given that properties such as pressure, temperature and viscosity are altitude-dependent). An example vector would be 0:5:80 

* <ins>geom_param:</ins> this is a vector containing important information regarding the geometrical properties of the 3D structure. This vector has length 4, with the first entry indicating the shape (2 – cone, 3 – sphere, 4 – rocket), the second entry the characteristic radius of the structure, the third entry the second characteristic length of the structure (such as the cone’s and rocket’s length or the outlet of the sphere), and then the fourth the number of suns the geometry will be subjected to. An example vector would be (3, 0.01, 0.005, 2), for a sphere of radius 1cm, outlet radius 5mm and 2 suns intensity.
  * geom_param(1): option, geometry chosen
  * geom_param(2): Ra, characteristic radius
  * geom_param(3): l, length of cone (option 2), outlet radius of the sphere (option 3), or length of the rocket (option 4)
  * geom_param(4): N, number of suns

* <ins>chan_param:</ins> this is a vector containing additional information regarding the channel properties of the 3D structure. This vector has length 6, with the first entry dictating the channel width A, the second the channel length B, the third the channel thickness L, the fourth the number of channels X, the fifth the channel spacing S, and lastly the sixth the ALD thickness t. An example vector would be (50x10^-6, 500x10^-6, 100x10^-6, 5, 50x10^-6, 50x10^-9). 
  * chan_param(1): A, channel width
  * chan_param(2): B, channel length
  * chan_param(3): L, channel thickness
  * chan_param(4): X, number of channels in cell
  * chan_param(5): S, channel spacing
  * chan_param(6): t, ALD thickness

**Outputs:** net_lift, fit, vft, deltaP, deltaT, vft2, vft3

* <ins>net_lift:</ins> this is the overall lift produced by the structure after accounting for its own weight

* <ins>fit:</ins> this is the estimated total force produced by the structure, based-off from the computed vft and vout

* <ins>vft:</ins> this is the computed flow-through velocity across the structure’s channels

* <ins>deltaP:</ins> this is the computed pressure differential across the structure’s channels

* <ins>deltaT:</ins> this is the computed temperature differential across the structure’s walls

* <ins>vft2:</ins> this is the computed flow-through velocity using an approximation for the deltaT term and neglecting the deltaP term

* <ins>vft3:</ins> this is the computed flow-through velocity neglecting the deltaP term but using the usual deltaT term

---
```
testing_Isun
```

**Type:** Script

**Description:** The purpose of this script is to calculate the dependency of the maximum flow-through velocity as a function of the sun intensity

**Inputs:** The user has to specify:

* geom_param(1) - option, geometry chosen (2, 3 or 4)
* geom_param(2) - Ra, characteristic radius(in m)
* geom_param(3) - Outlet radius l (in m)
* chan_param(6) – t, ALD thickness (in m)
* altitude – Altitude Vector (in km)
 
   The user has to also narrow the optimization range, which is set as a logarithmic spacing. This is given by: <br />
   
   L_vec = logspace(log10(1x10^-6),log10(10^-2),100) - for channel thickness L <br />
   A_vec = logspace(log10(1x10^-8),log10(500x10^-5),100) - for channel width A <br />
   I_vec = logspace(log10(1),log10(100),50) - for the number of suns N <br />
   
   As can be seen, the code will go through a triple for loop.

**Outputs:** N/A



## Authors

Thomas Celenza, Andy Eskenazi and Igor Bargatin <br />
University of Pennsylvania, 2022 <br />
