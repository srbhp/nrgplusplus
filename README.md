# nrgplusplus

[![.github/workflows/cmake.yml](https://github.com/srbhp/nrgplusplus/actions/workflows/cmake.yml/badge.svg?branch=main&event=push)](https://github.com/srbhp/nrgplusplus/actions/workflows/cmake.yml)

This is a modern C++ implementation of the 
Numerical Renormalization Group (NRG) Method 
to  solve Quantum Impurity problems embeded in a Bath. 
We tried to follow the convention of 
[Bulla, 2008](https://doi.org/10.1103/RevModPhys.80.395). 

##  Documentation
https://srbhp.github.io/nrgplusplus/


### Quick Start Guide
https://srbhp.github.io/nrgplusplus/build.html




## Example : Single Impurity Anderson Impurity (SIAM)
`(See : examples/rgflowSIAM/main.cpp)`

Define the impurity Model wth onsite energy `eps` and Coulomb energy `U_int`. 

```cpp
  spinhalf impurity(eps, U_int);
```

The bath for the SIAM is also constructed in the same way. 

```cpp
  spinhalf bathModel(0, 0); // set parameters
```

Once we have created the bath and the impurity we can 
construct a `nrgcore` object which will take care of 
many things that we need to do for the `NRG iterations`.
This includes calculating static and dynamic quantities
of the Impurity.

```cpp
  nrgcore<spinhalf, spinhalf> siam(impurity, bathModel);
  siam.set_parameters(1024);       // set max number of states to be kept
  siam.add_bath_site({V, V}, 1.0); // V is the coupling og the impurity and first bath site. 
  siam.update_internal_state();
```
Next we iteratively add the bath sites in the same way. We 
also create HDF5 file object to save the Eigenvalues.


```cpp
  // file where outputput will be wriiten
  h5stream::h5stream rfile("resultSIAM.h5");
  // Iterative add bath sites
  for (int in = 0; in < nMax; in++) {
    double rescale = 1.0;
    if (in > 0) {
      rescale = std::sqrt(LAMBDA);
    }
    siam.add_bath_site({hopping(in, LAMBDA), hopping(in, LAMBDA)}, rescale);
    // Update System Operators now here if we need to.
    // This has to be done before updating the systems internal state.
    siam.update_internal_state();
    // Save the eigenvalue of the current iteration
    rfile.write(siam.all_eigenvalue, "Eigenvalues" + std::to_string(in));

  }
  rfile.close();
```
Plot the RG flow `(See : examples/rgflowSIAM/plot.py)`.
![SIAM RG Flow](docs/image/rgflowSIAM.png)

## Contributing

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.

Please make sure to update tests as appropriate.

## Similar Open-source Projects 

  1. http://www.phy.bme.hu/~dmnrg/
  2. https://github.com/rokzitko/nrgljubljana

## Acknowledgment

During the development of this project I was employed by
  - Prof. Jonas Fansson, Uppsala University and 
  - Prof. Frithjof Anders, Technical University Dortmund.

This project would have been not possible without their support. 
