## NRGCore
This is the main repository of the NRG library. 

Use this repository as a `git submodule` to build a NRG model.


 
## Dependence 
It uses LAPACKE (http://www.netlib.org/lapack/lapacke.html) for
 
 matrix operations such as diagonalization, and multiplication.  It also uses OpenMP
 for parallel computing.  
 
 
## RG flow
 
 Energy flow Renormalization group: This figure should be compared with Fig.4 of the RMP paper.
![alt text](src/c++/rgflow/out.png)

### How to build a Model

`examples1` directory contains examples of different Models. 
1. Copy a folder here from `examples` to start woking on a specific model.
1.  Edit `CMakeLists.txt` to add the model folder.
1. `build` the project by invoking `make`.

