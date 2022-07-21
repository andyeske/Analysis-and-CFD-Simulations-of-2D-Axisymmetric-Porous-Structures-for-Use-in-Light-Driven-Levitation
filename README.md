<a name="back_to_top"></a>
# 3D photophoretic microflyers require ultrathin, ultralight porous materials

In this repository, we present the various MATLAB functions and scripts utilized in our work, alongside an explanation of their use, inputs and outputs.

## Summary

The MATLAB code presented in this repository accompanies the paper titled "3D photophoretic microflyers require ultrathin, ultralight porous materials", a work that proposes that photophoretic microflyers would greatly benefit from a three-dimensional (3D) geometry that uses Knudsen pump to create a high-speed jet. This code is based on the theory developed on this paper's supplementary information section, and essentially seeks to find the optimal set of both geometrical and porous parameters that would return the highest payload capabilities for various types of 3D geometries, such as cones, spheres and rockets. This repository contains the code to calculate the lift forces produced by each geometry (calc_F), as well as several scripts that perform parametric studies where two or three of the main variables are modified.

## MATLAB Code

### Table of Contents

1. [ calc_F ](#calc_F)
2. [ testing_Isun ](#testing_Isun)
3. [ testing_L_A ](#testing_L_A)
4. [ testing_t_A ](#testing_t_A)
5. [ testing_L ](#testing_L)
6. [ new_param_sweep ](#new_param_sweep)
7. [ plotting_forces ](#plotting_forces)

---
<a name="calc_F"></a>
### calc_F

**Type:** Function

**Description:** The purpose of this code is to provide a framework for calculating the flow-through velocity (vft) and outlet velocity (vout) of an arbitrary 3D geometry.

**Inputs:** 
```
altitude, geom_param, chan_param
```

* <ins>altitude:</ins> this is a vector containing the specific altitudes (in km) at which the vft will be calculated at (given that properties such as pressure, temperature and viscosity are altitude-dependent). An example vector would be 
   ```
   0:5:80 
   ```
   for an array of altitudes ranging from 0 to 80 km in multiples of 5.

* <ins>geom_param:</ins> this is a vector containing important information regarding the geometrical properties of the 3D structure. This vector has length 4, with the first entry indicating the shape (2 – cone, 3 – sphere, 4 – rocket), the second entry the characteristic radius of the structure, the third entry the second characteristic length of the structure (such as the cone’s and rocket’s length or the outlet of the sphere), and then the fourth the number of suns the geometry will be subjected to. 
   ```
   geom_param(1): option, geometry chosen
   geom_param(2): Ra, characteristic radius
   geom_param(3): l, length of cone (option 2), outlet radius of the sphere (option 3), or length of the rocket (option 4)
   geom_param(4): N, number of suns
   ```
   An example vector would be 
   ```
   [3, 0.01, 0.005, 2]
   ```
   for a sphere of radius 1cm, outlet radius 5mm and 2 suns intensity.
   

* <ins>chan_param:</ins> this is a vector containing additional information regarding the channel properties of the 3D structure. This vector has length 6, with the first entry dictating the channel width A, the second the channel length B, the third the channel thickness L, the fourth the number of channels X, the fifth the channel spacing S, and lastly the sixth the ALD thickness t. 
    ```
   chan_param(1): A, channel width
   chan_param(2): B, channel length
   chan_param(3): L, channel thickness
   chan_param(4): X, number of channels in cell
   chan_param(5): S, channel spacing
   chan_param(6): t, ALD thickness
   ```
   An example vector would be 
   ```
   [50x10^-6, 500x10^-6, 100x10^-6, 5, 50x10^-6, 50x10^-9]
   ``` 
   for a channel of width radius 50 microns, length 500 microns, thickness 100 microns, with 5 channels in total and spacing of 50 microns. 
   The ALD thickness in this case would be 50 nm.
   
  
**Outputs:** 
   ```
   net_lift, fit, vft, deltaP, deltaT, vft2, vft3, aerial
   ```

* <ins>net_lift:</ins> this is the overall lift produced by the structure after accounting for its own weight
* <ins>fit:</ins> this is the estimated total force produced by the structure, based-off from the computed vft and vout
* <ins>vft:</ins> this is the computed flow-through velocity across the structure’s channels
* <ins>deltaP:</ins> this is the computed pressure differential across the structure’s channels
* <ins>deltaT:</ins> this is the computed temperature differential across the structure’s walls
* <ins>vft2:</ins> this is the computed flow-through velocity using an approximation for the deltaT term and neglecting the deltaP term
* <ins>vft3:</ins> this is the computed flow-through velocity neglecting the deltaP term but using the usual deltaT term
* <ins>aerial:</ins> this is the aerial densities that can be sustained by the structure's material at each altitude

[ back to top ](#back_to_top)

---
<a name="testing_Isun"></a>
### testing_Isun

**Type:** Script

**Description:** The purpose of this script is to calculate the dependency of the maximum flow-through velocity as a function of the sun intensity

**Inputs:** The user has to specify:

   ```
   geom_param(1) - option, geometry chosen (2, 3 or 4)
   geom_param(2) - Ra, characteristic radius(in m)
   geom_param(3) - Outlet radius l (in m)
   chan_param(6) – t, ALD thickness (in m)
   altitude – Altitude Vector (in km)
   ```
 
   The user has to also narrow the optimization range, which is set as a logarithmic spacing. This is given by: 
   ```
   L_vec = logspace(log10(1x10^-6),log10(10^-2),100) - for channel thickness L 
   A_vec = logspace(log10(1x10^-8),log10(500x10^-5),100) - for channel width A 
   I_vec = logspace(log10(1),log10(100),50) - for the number of suns N 
   ```
   As can be seen, the code will go through a triple for loop.

**Outputs:** N/A

**Plots:** 
<p float="left">
<img src="https://github.com/andyeske/MATLAB-fluidflow-parametric-studies/blob/main/Sample%20Plots/testing_Isun1.png" width="400"> 
<img src="https://github.com/andyeske/MATLAB-fluidflow-parametric-studies/blob/main/Sample%20Plots/testing_Isun2.png" width="400">
</p>

---
<a name="testing_L_A"></a>
### testing_L_A

**Type:** Script

**Description:** The purpose of this script is to calculate the dependency of the maximum flow-through velocity as a function of L, the channel thickness, and A, the channel width

**Inputs:** The user has to specify:

   ```
   geom_param(1) - option, geometry chosen (2, 3 or 4)
   geom_param(2) - Ra, characteristic radius(in m)
   geom_param(3) - Outlet radius l (in m)
   geom_param(4) - N, the number of suns
   chan_param(6) – t, ALD thickness (in m)
   altitude – Altitude Vector (in km)
   ```

   The user has to also narrow the optimization range, which is set as a logarithmic spacing. This is given by 
   ```
   L_vec = logspace(log10(1x10^-6),log10(10^-2),100) - for channel thickness L 
   A_vec = logspace(log10(1x10^-8),log10(500x10^-5),100) - for channel width A 
   ```

   As can be seen, the code will go through a double for loop. However, notice that this code assumes that B = 10A and S = A. <br />

**Outputs:** Right now, the code has no specific outputs, although it can very easily be modified to return back the combination of A and L that yielded the max net lift, max force, or max vft

**Plots:** 

![testing_L_A](https://github.com/andyeske/MATLAB-fluidflow-parametric-studies/blob/main/Sample%20Plots/testing_L_A.png) 

---
<a name="testing_t_A"></a>
### testing_t_A

**Type:** Script

**Description:** The purpose of this script is to calculate the dependency of the maximum flow-through velocity as a function of t, the ALD thickness, and A, the channel width

**Inputs:** The user has to specify:
   
   ```
   geom_param(1) - option, geometry chosen (2, 3 or 4)
   geom_param(2) - Ra, characteristic radius(in m)
   geom_param(3) - Outlet radius l (in m)
   geom_param(4) - N, the number of suns
   chan_param(3) – L, the cannel thickness (in m)
   altitude – Altitude Vector (in km)
   ```

   The user has to also narrow the optimization range, which is set as a logarithmic spacing. This is given by 
   ```
   vec_t = logspace(log10(1x10^-9),log10(100x10^-9),1000) – for ALD thickness t 
   vec_A = logspace(log10(1x10^-7),log10(50x10^-5),1000) – for channel width A 
   ```

   As can be seen, the code will go through a double for loop. Notice that the code assumes that B = 10A and S = A. <br />

**Outputs:** Right now, the code has no specific outputs, although it can very easily be modified to return back the combination of A and t that yielded the max net lift, max force, or max vft

**Plots:** 

![testing_t_A](https://github.com/andyeske/MATLAB-fluidflow-parametric-studies/blob/main/Sample%20Plots/testing_t_A.png) 

---
<a name="testing_L"></a>
### testing_L

**Type:** Script

**Description:** The purpose of this script is to calculate the dependency of the maximum flow-through velocity as a function of L, the channel thickness

**Inputs:** The user has to specify:

   ```
   geom_param(1) - option, geometry chosen (2, 3 or 4)
   geom_param(2) - Ra, characteristic radius(in m)
   geom_param(3) - Outlet radius l (in m)
   geom_param(4) - N, the number of suns
   chan_param(1) – A, channel width (in m)
   chan_param(2) – B, channel length (in m)
   chan_param(5) – S, channel spacing (in m)
   chan_param(6) – t, ALD thickness (in m)
   altitude – Altitude Vector (in km)
   ```

   The user has to also narrow the optimization range, which is set as a logarithmic spacing. This is given by
   ```
   vec = (10x50x10^-9):(1x10^-8):(10^-2) - for channel thickness L
   ```

   As can be seen, the code will go through a double for loop

**Outputs:** Right now, the code has no specific outputs, although it can very easily be modified to return back the values of L that yielded the max net lift, max force, or max vft

**Plots:** 

![testing_L](https://github.com/andyeske/MATLAB-fluidflow-parametric-studies/blob/main/Sample%20Plots/testing_L.png) 

---
<a name="new_param_sweep"></a>
### new_param_sweep

**Type:** Script

**Description:** The purpose of this script is to calculate the dependency of the maximum flow-through velocity as a function of L, the channel thickness, A, the channel width, and l, a geometrical parameter controlling Ain

**Inputs:** The user has to specify:
   
   ```
   geom_param(1) - option, geometry chosen (2, 3 or 4)
   geom_param(2) - Ra, characteristic radius(in m)
   geom_param(4) - N, the number of suns
   chan_param(1) – A, channel width (in m)
   chan_param(2) – B, channel length (in m)
   chan_param(4) – X, number of channels (in m)
   chan_param(6) – t, ALD thickness (in m)
   altitude – Altitude Vector (in km)
   ```

   The user has to also narrow the optimization range, which is set as a logarithmic spacing. This is given by  
   ```
  min_A = 1*10^-8;
  max_A = 500*10^-5;
  num_A = 100; % number of points in the vector
  A_vector = logspace(log10(min_A),log10(max_A),num_A);
  
  min_L = 1*10^-6;
  max_L = 10^-2;
  num_L = 100; % number of points in the vector
  L = logspace(log10(min_L),log10(max_L),num_L);

  min_ra = 10^-4;
  max_ra = 10^-2;
  num_ra = 100; % number of points in the vector
  ra = logspace(log10(min_ra),log10(max_ra),num_ra);
   ```

   As can be seen, the code will go through a double for loop. Notice that the code assumes that B = 10A and S = A. 

**Outputs:** The code returns as its main output the minimum altitude at which a positive payload was carried (alongside that payload) and the maximum payload carried (alongside the altitude at which that is possible). In both cases, the combination of A, L and l that was able to produce such structures is returned. For example:

> For the given constraints, the minimum altitude at which the net payload becomes positive is 30 km 
> and has the potential of carrying 0.65669 mg as its maximum payload 
> The corresponding parameters that yielded these results are A = 6.2996e-05 m, L = 1e-06 m, and Ra_{out} = 0.0079248 m 
> Additionally, the maximum payload found was 368.3004 mg and corresponded to an altitude of 70 km 
> The corresponding parameters that yielded these results are A = 0.001019 m, L = 3.1257e-05 m, and Ra_{out} = 0.01 m 
> Finally, the aerial density corresponding to this maximum payload capability was 120.7 g/m2 

**Plots:** 

<p float="left">
<img src="https://github.com/andyeske/MATLAB-fluidflow-parametric-studies/blob/main/Sample%20Plots/Sphere_Pay3.jpg" width="300">
<img src="https://github.com/andyeske/MATLAB-fluidflow-parametric-studies/blob/main/Sample%20Plots/Sphere_Pay4.jpg" width="300">
<img src="https://github.com/andyeske/MATLAB-fluidflow-parametric-studies/blob/main/Sample%20Plots/Aereal_Sphere.jpg" width="300">
</p>

---
<a name="plotting_forces"></a>
### plotting_forces

**Type:** Script

**Description:** The purpose of this script is to calculate the dependency of the maximum flow-through velocity as a function of L, the channel thickness, A, the channel width, and l, a geometrical parameter controlling Ain

**Inputs:** The user has to specify:

   ```
   geom_param(1) - option, geometry chosen (2, 3 or 4)
   geom_param(2) - Ra, characteristic radius(in m)
   geom_param(3) - Outlet radius l (in m)
   geom_param(4) - N, the number of suns
   chan_param(1) – A, channel width (in m)
   chan_param(2) – B, channel length (in m)
   chan_param(3) – L, channel thickness (in m)
   chan_param(5) – S, channel spacing (in m)
   chan_param(6) – t, ALD thickness (in m)
   altitude – Altitude Vector (in km)
   ```

**Outputs:** Right now, the code has no specific outputs, although it can very easily be modified to return back the values of altitude that yielded the max net lift, max force, or max vft

**Plots:** 

![plotting_forces](https://github.com/andyeske/MATLAB-fluidflow-parametric-studies/blob/main/Sample%20Plots/plotting_forces.png) 

## Authors

Tom Celenza, Andy Eskenazi and Igor Bargatin <br />
Department of Mechanical Engineering and Applied Mechanics <br />
University of Pennsylvania, 2022 <br />
