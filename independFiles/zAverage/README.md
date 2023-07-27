### z-average for a bath with energy dependent DOS
 
The bath Hamiltonian is

$$ H = \sum_k \epsilon_k c_{k}^\dag c_{k} $$

$$  = \sum_{j=1}^{N} E_j^z(b^\dag_{j0}b_{j0} - b^\dag_{-j0}b_{-j0}  ) $$

Here 2N is number of discretization points of the conduction bath. 

$$E_j = \frac{\int_{I_j}  d\epsilon}{\int_{I_j}\frac{ d\epsilon}{\epsilon}}$$

$I_J = [\epsilon_{j-1} , \epsilon_{j}]$ 

The coupling hamiltonian for the first Wilson's site 
is 

$$ A (f_0^\dag c_{imp} + h.c.) $$

Here A ensures the $\{f_0^\dag, f_0\} = 1$ and $f_0$ is defined as 


$$  f_0 = \sum_{j}\sqrt{\frac{|\epsilon_j - \epsilon_{j+1} |}{2D}} (b_{j0} + b_{-j0})$$


Starting with these definitions, We start the Lanczos algorithm

### Lanczos algorithm

![Lanczos](./lanczos.png)



---
[1.] Campo, V. L. ; Oliveira, L. N.: Alternative discretization in the numerical renormalization-group
method. In: Phys. Rev. B 72 (2005), p. 104432. doi: 10.1103/PhysRevB.72.104432
