# Code for generating the data of "Polarization-induced stress in the noisy voter model"

The scripts and the code needed to reproduce the results of the article  "Polarization-induced stress in the noisy voter model" (arXiv:2402.00516). 

A program in C lenguaje can be found in the directory code/.

To compile the code one have to run the comand "make".
The number of agents of the simulation (N) is a compilation parameter and can be modified in the Makefile.  
Once the program is compiled, you can run it with the usual "./stressed_mf_gille_N??".
The input parameters of the program are: 
1. Number of calls to the Gillespie rutine.
2. Random numer generator seed (if seed=0, the seed is taken from the system clock).
3. Initial fraction of agents in the state 1.
4. The discretization step of the grid.

The program produce two kind of output files:
a) "histogram_N??_eMCS?_a???_b??.dat"--> The histogram for a given value of the parameters a/h and b/h.
b) "fase_N50_eMCS6_delta0.1000.dat"--> The phase (W,U,B or M) at wich the system is for each pair of values (a/h, b/h).

The phases are codified using the following dictionary:
//Fase Dictionary:
//W: fase=0
//B: fase=1
//U: fase=2
//M: fase=3
