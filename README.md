### warning :: Work in progress 

# kondo-nrg
This is a Numerical Renormalization Group(NRG) code for the single site Anderson 

 impurity model. I have followed Rev. Mod. Phys. 80, 395(2008)

(https://doi.org/10.1103/RevModPhys.80.395) . 
 
## Dependence 
It uses LAPACKE (http://www.netlib.org/lapack/lapacke.html) for
 
 matrix operations such as diagonalization, and multiplication.  It also uses OpenMP
 for parallel computing.  
 
 
## RG flow
 
 Energy flow Renormalization group: This figure should be compared with Fig.4 of the RMP paper.
![alt text](src/c++/rgflow/out.png)
