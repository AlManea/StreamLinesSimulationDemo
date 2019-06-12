# StreamLinesSimulationDemo

## Approach and Assumptions 

We assumed that the system is immiscible and of a unit mobility ratio, such that there is no need to calculate the BL solution along each streamline to find the water front, as it would be simply found by tracing the initial particles travelling from the injection well along each streamline. 

We used the infinity (∞-norm) to calculate the errors.
We launch streamlines from a uniformly distributed points on a circle of unit-radius around the wellbore. 
We used the following formula to convert the calculated velocity into volumetric flux:

<h5 align=center> q = ∅Av = ∅[2πh/N<sub>SL</sub> ]v,

where: 

q is the volumetric flow rate,

∅ is the porosity,

A is the cross-sectional area,

v is the interstitial velocity,

h is the reservoir thickness, and 

N<sub>SL</sub>  is the number of streamlines.
		
For reasonable computational times, we use a maximum value on the number of points that every streamline can have. Thus, when a streamline exceeds that value, it is chopped. However, this does NOT mean that we could not trace that streamline. It only means that we will need extremely long processing time to trace it, as the particles move very slowly in that streamline direction.

The code uses a simple API, and it is very generic and easy to run from MATLAB’s prompt, as follwos:

SL_Demo(x, y, q, sl);

where:

x is the x coordinates of the wells, starting with the injector at location (0,0)
  
y is the y coordinate of the wells, starting with the injector at location (0,0)

q is the well rates, starting with the injector

sl is the number of desired streamlines/well.
