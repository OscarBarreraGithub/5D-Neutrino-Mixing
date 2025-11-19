# Bessel Functions

This solver finds the Kaluza–Klein (KK) masses in Randall–Sundrum models by solving the Bessel-function boundary conditions that quantize 5D bulk fields. The bulk equations reduce to Bessel’s equation in conformal (z) coordinates; the UV/IR boundary conditions give transcendental equations whose roots are the KK masses.

### Inputs

* **Geometry / scale**

  * $k$ and $\Lambda \equiv ke^{-\pi k r_c}$ in the format of [`warpConfig`](warpConfig)
 
 
* **Field & boundary data**

  * **Species**: `gauge` or `fermion`.
  * **BC choice**:
    * `gauge`: Must use `NN`.
    * `fermion`: Must use `++` (LH zero mode) or `--` (RH zero mode).
  * **Fermions only**: bulk mass parameter (c) (dimensionless), with $\alpha \equiv |c+\tfrac12|$.
* **Numerics**

  * Number of roots to find (N).
  * `exact=True|False` (exact ratio equation vs. IR-only approximation).
  * Root-finding tolerances and bracketing options.


### Outputs

* **Mass eigenvalues** $m_n$
  Expect $m_n = x_n \Lambda$, where $x_n$ are dimensionless roots.
* **(Optional) Mode data**

  * $b_n$: the Bessel mix $J_\nu + b_n Y_\nu$ fixed by BCs.
  * $N_n$: normalization constants for profiles.


## Where the Bessel functions enter

Bulk EOMs in the RS background reduce to

$$
\left[z^2 \partial_z^2 + z \partial_z + (m^2 z^2 - \nu^2)\right]\Phi(z)=0
$$

$$
\quad\Rightarrow\quad \Phi(z) = A J_\nu(m z)+B Y_\nu(m z).
$$

Boundary conditions tie $A$ & $B$ together and quantize $m$. It is numerically cleaner to work with dimensionless

$$
x \equiv m z_\nu \quad\Rightarrow\quad m=\frac{x}{z_\nu}=x\Lambda.
$$

------
### Algorithm (high level)

1. Set the scale $\Lambda$
2. **Choose equation:** build $F(x)$ from the exact ratio condition above (or $F(x)=J_\nu(x)$ in IR-only mode).
3. **Seed intervals:** bracket roots near the tabulated zeros of $J_\nu$.
4. **Solve:** use a bracketed method (e.g., Brent) to find $x_n$ with tolerance $\delta x$.
5. **Map to masses:** $m_n = x_n / z_\nu = x_n \Lambda$.


----

### 1. Numerical Stability

The standard condition for KK modes often looks like a ratio:

$$
 \frac{J_\nu(x)}{Y_\nu(x)} = \frac{J_\nu(\epsilon x)}{Y_\nu(\epsilon x)} 
$$

Directly implementing this is numerically dangerous for two reasons:
1.  **Singularities**: $Y_\nu(z)$ has zeros, causing the ratio to blow up to infinity.
2.  **Small Argument Divergence**: As $z \to 0$, $Y_\nu(z) \to -\infty$.

The solver instead finds roots of the "cross-product" equation:

$$
F(x) = J_\nu(x) Y_\nu(\epsilon x) - J_\nu(\epsilon x) Y_\nu(x) = 0 
$$

This function behaves well when $\epsilon$ is extremely small, which is typical for RS models.

### 2. Intelligent Seeding

Finding roots of oscillating functions is tricky; if you guess wrong, you might find the 2nd root when you wanted the 1st.
The solver uses a robust seeding strategy:
*   **Integer Order**: If $\nu$ is an integer, we use `scipy.special.jn_zeros` to get the exact locations of the zeros of $J_\nu(x)$.
*   **Non-Integer Order**: We use the Bessel asymptotic expansion for large arguments:

    $$
     x_n \approx \left(n + \frac{\nu}{2} - \frac{1}{4}\right)\pi 
     $$
    

These seeds are used to create **brackets** (intervals $[a, b]$ where the sign of $F(x)$ changes). We then pass these brackets to `scipy.optimize.brentq` (Brent's method), which is guaranteed to converge if a sign change exists.

### 3. The "IR-Only" Approximation

The code supports an `exact=False` mode. In RS, the hierarchy $\epsilon = \Lambda/k$ is tiny ($\sim 10^{-15}$).
*   For $x \sim O(1)$, the argument $\epsilon x$ is effectively zero.
*   $J_\nu(\epsilon x) \to 0$ (for $\nu > 0$) and $Y_\nu(\epsilon x) \to -\infty$.
*   The cross-product term dominated by $Y_\nu(\epsilon x)$ forces $J_\nu(x)$ to be very close to 0.

Thus, solving $J_\nu(x) = 0$ is a very nice approximation. The solver includes this mode for speed, though the exact mode is fast enough for most purposes.

