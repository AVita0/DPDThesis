# Adaptive Digital Predistortion
This is the source code and manuscript for my masters thesis on adaptive digital predistortion (DPD).

This work was completed while at The University of New Haven.

## Background 

Radio Frequency (RF) Power Amplifiers (PA) are fundamentally non-linear due to their design 
and environmental effects, this non-linearity reduces performance and introduces distortion to 
the system. To improve the linear region of operation an adaptive Digital Predistortion Process 
(DPD) is used consisting of an adaptive Finite Impulse Response (FIR) filter utilizing the Least 
Means Square (LMS) algorithm.

This process pre distorts the PA input to be the inverse of the amplifier distortion. As the pre distorted input 
passes though the PA and the distortion is added the output becomes increasingly linear as the 
distortion has been accounted for.

![DPDHighLevel](https://github.com/AVita0/DPDThesis/assets/80707753/e2bc5cbf-0233-4129-b160-a27f5aadb6e1)

## Block Diagram 

![DPD_BlockDiagram](https://github.com/AVita0/DPDThesis/assets/80707753/eb12eea6-e223-4091-9e4d-8cc5fd5ce0be)

## Simulations 

The two source codes in this repository are DPD_Heavy and DPD_Slim. The source code was developed in MatLAB R2022b.

DPD_Heavy contains all of the test cases and results provided in the manuscript, it is also not very organized. 

DPD_Slim contains one test case, is well organized, and easier to understand. 
