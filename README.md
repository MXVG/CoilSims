# Coilgun-Simulator
Simulation using FEMM to simulate a iron slug accelerated by a magnetic field induced by current in inductors part of a RLC circuit.

Simulator is intended for personal use to help with circuit component selection for building a [coilgun](https://en.wikipedia.org/wiki/Coilgun)

Just run python file in folder with StageSetup.xlsx (your coil configuration) to start simulation

# Coil Configuration (StageSetup.xslx)

![Screenshot from 2024-01-19 18-34-31](https://github.com/MXVG/Coilgun-Simulator/assets/91643168/0ca3895f-cffe-4d91-8c0a-7499452bd479)

## W_mm
Diameter of wire with copper enamel

## CMAT
Gauge of wire used in inductor

## D_ic
Inner diameter of inductor

## D_oc
Outer diamter of inductor

## L_c
Length of inductor

## N_Turns
Number of turns in the inductor

## R
Resistance of inductor

## L
Inductance

## C
Total capacitance of capacitors (An RLC circuit is assumed to give current to inductor)

## V
Voltage level of capacitors at time of discharge through RLC circuit

## Z_proj
How far back metal rod starts from simulated inductor in FEMM (mainly used to adjust for different coil lengths between RLC stages)

## ESR
Not used.


# Examples
Folders contain example output

# Dependencies
FEMM: https://www.femm.info/wiki/HomePage 
`pip install pyfemm`

Matplotlib: `pip install matplotlib`

Numpy: `pip install numpy`
