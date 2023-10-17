# VO2 Thermal Sim

The program was used  to simulate the thermal perofrmance of Vanadium Dioxide. Results from the simulation have been published in a paper as seen here: https://pubs.acs.org/doi/full/10.1021/acs.nanolett.3c02251



# Program Details

The program uses Finite Difference Methods (FDM) With discretized timesteps. The program takes into account 3 main factors:
  1) Thermal conduction of heat within the material
  2) Heat added from an external laser
  3) Internal Joule heating from electricity flowing through the material

The resistance of Vanadium Dioxide (VO2) is extremely nonlinear (We modelled it using a 5th order polynomial) so we were unable to use alternative methods of thermal simulation which are faster.

Thermal conduction is modelled using matrix addition
Electrical simulation models discretized resistor cubes and does power calculations uding pySpice (A python extension for SPICE software)

## Program Speed

Because of the simulation method used, The program inherently takes a long time to run. As such, we have two methods of running the simulation

### 1) Thermal and Electrical simulation
   - Here, we use both simulation methods. As such, we are left with standard CPU processing for thermal conduction.
   - The main bottleneck for this method is actually moving data between different arrays as matrices can get extremely large
     
### 2) Only Thermal Simulation
   - If we assume little impact from electrica heating (For example, towards the beginning of the simulation) we can turn off the electrical simulation segments.
   - Here, we can take advantage of GPU acceleration. Using the pytorch library, we can use GPU acceleration to drastically increase the it/s of the FDM loop. Below are the options for simulation modes, which can be adjusted in the code:
   
#### Setup device for PyTorch calculation
   1) cpu: Uses the computer's CPU.
        - This is required if you want to also have electrostatic sim
        - Can support 64-bit floating point
        - Generally slower than the other two ptions
  
   1) cuda: Uses Nvidia's CUDA GPU processors
        - This is only available on GPUs with Nvidia graphics (unfortunately, I cannot afford one)
        - This supports 64-bit floating point
        - Probably the fastest option
  
   3) mps: Uses apple's metal GPU acceleration
        - This is what is used to test GPU acceleration
        - Only supports 32-bit floating point (sadge) (Only used for testing, not for simulation)
        - Similar to cuda but 32-bit is not accurate enough for our puroses
      
   For testing purposes, I am usually using gpu. For running the actual trial, its probably best to use CUDA if available for ~7.5x speed increase from my testing
