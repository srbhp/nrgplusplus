## ⚠ This is a work in progress repository


##  Documentation
![srbhp.github.io/nrgplusplus](https://srbhp.github.io/nrgplusplus/)



 
## Dependence 

- `intel-mkl` : Lapack  
- `hdf5` : We use hdf5 file for the I/O
- to genarate documentation we these packages 
    - `sphinx-rtd-theme, breathe,sphinx-sitemap,sphinx,exhale,alabaster `
 
## RG flow
 
 Energy flow Renormalization group: This figure should be compared with Fig.4 of the RMP paper.
![alt text](src/c++/rgflow/out.png)

### How to build a Model

The `examples` directory contains examples of different Models. 

- Copy a folder from `examples` to start working on a specific model.
- Edit `CMakeLists.txt` to add the model folder.
- `build` the project by invoking `make`.

