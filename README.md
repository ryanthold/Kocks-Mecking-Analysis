# Kocks-Mecking Analysis
This project is an open-source MATLAB code with functions for extracting tensile testing properties and producing Kocks-Mecking (KM) curves for work hardening analysis.

### __For further understanding of this code and to cite its use, view "CITE PAPER".__

Its main functions are to:
- Smooth the stress-strain curve to limit noise for KM analysis
- Estimate tensile properties such as the proportional limit, instability point, 0.2% offset yield stress, Ultimate tensile strength, and elongation
- Estimate work hardening coefficient values for stage 3 and stage 4 of plastic deformation using KM analysis

All functions are provided with explanation headers.  Additionally a GUI is provided for specific smoothing operations.

Functions:
- **SS_KMsmoothapp.m** - GUI application
- **executioncode.m** - an example that executes all functions listed below for many datasets
- **cutend.m** - removes datapoints on the stress-strain curve that occur after fracture
- **spline_fit.m** - creates a smoothing spline fit of stress-strain data
- **ESS_to_TSS** - converts engineering stress-strain data to true-strain data
- **PL02** - determines the proportional limit, 0.2% offset yield stress, and elastic modulus
- **instability.m** - determines the instability point through Considere's criterion
- **wrkhard.m** - calculates work hardening rate, which is the instantaneous slope of true stress-strain
- **cb3_calc** - determines the work hardening coefficient for stage 3
- **cb4_calc** - determines the work hardening coefficient for stage 4

*MATLAB's Curve Fitting Toolbox, Signal Processing Toolbox, and Deep Learning Toolbox are utilized in this work.*

This work was in coordination with:
1. Sandia National Laboratories, Livermore, CA 94550
2. Texas A&M, Department of Materials Science and Engineering, Department of Mechanical Engineering, College Station, TX 77843
3. University of California, Irvine, Department of Materials Science and Engineering, Henry Samueli School of Engineering, Irvine, CA 92697
