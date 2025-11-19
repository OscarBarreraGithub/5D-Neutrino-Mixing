
This directory provides utilities to compute derived 5D warp parameters and the fermion zero-mode overlap
functions used throughout the project.

-----

Start by choosing the two physical inputs you care about:

- $k$: the AdS curvature (mass dimension, in GeV). We usually take this to
  be of order the Planck mass, e.g. $k \sim M_{\text{Pl}}$.
- $\Lambda$ : the IR / KK mass scale you want (in GeV). Default is 3 TeV.

The code then computes the warp factor $\epsilon$ and the orbifold radius $r_c$ from the relation

$$
\Lambda = k e^{-\pi k r_c} = k\,\epsilon.
$$

With those quantities available you can immediately compute
$f_{IR}(c, \epsilon)$ and $f_UV(c, \epsilon)$ as defined in `wavefuncs.py`

The IR “overlap factor” for a canonically normalized fermion zero mode is

$$
f_{\mathrm{IR}}^2=\frac{\tfrac12 - c}{1 - \epsilon^{1-2c}}
$$

while for the UV we have

$$
\left(f_N^{\mathrm{UV}}\right)^2=\frac{\tfrac12 - c}{\epsilon^{2c-1}-1}
$$

Note: the `f_IR` function returns $f_{\mathrm{IR}}$, not $f_{\mathrm{IR}}^2$



What is in this directory
-------------------------

- ``baseParams.py``
  - The central helper: call ``get_warp_params(k, Lambda_IR)`` and it returns a
    dictionary containing the derived parameters: ``epsilon``, ``rc``,
    ``warp_log`` (= $\pi k r_c$), and the conformal brane positions. Input validation is performed and the implementation is numerically stable to small $\epsilon$.

- ``wavefuncs.py``
  - Implements the zero-mode overlap functions ``f_IR(c, epsilon)`` and
    ``f_UV(c, epsilon)``. Both accept scalars or NumPy arrays for ``c`` and
    treat the special ``c = 1/2`` limit analytically to avoid division by
    zero.

- ``wavefuncsTest.ipynb``
  - Small interactive notebook that demonstrates reloading the local modules,
    computing overlaps for a range of ``c``, and plotting the results for a
    few choices of ``Lambda``.


Updating $\Lambda$ and recomputing
-----------------------------------
It is intentionally simple to change the IR scale and recompute the derived
quantities. For a Monte-Carlo sweep or a script, just call ``get_warp_params``
for each desired value of ``Lambda_IR`` and pass the returned ``epsilon`` into
``f_IR``/``f_UV``. The functions are lightweight and vectorized, so a loop or
vectorized mapping will be efficient for typical scans.


