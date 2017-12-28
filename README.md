# CompositeSandwichOptimizer
Routines for optimization of composite sandwich based on classical lamination theory

This project consists of a Matlab script (sandwich_optimizer.m) in which it is possible to define a composite material sandwich geometry (supposed a rectangle) and material characteristics.

The other input to the file is the maximum force (acting at the middle of the panel).

The script automatically calculates the ply failures, for every combination of face and core thickness (using the energetic approach). At the end the script plots the minimum face and core thicknesses that complies with the maximum stress calculated. At the same time it plots the total weight for the solution and suggests the thickness combination that results in minimum panel weight.

An example of the usage of it is in the "Carbon_skate.pdf" file
