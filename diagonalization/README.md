The 5D model contains three kinds of lepton fields living in the bulk:

* L: Left-handed SU(2) doublets (contain $\nu_L$ and $e_L$),
* E: right-handed charged singlets ($e_R$),
* N: right-handed neutral singlets (sterile neutrinos).

Each has a bulk mass and couples on the IR brane to a “Higgs-like” field through Yukawa matrices $Y_E$ and $Y_N$.  In addition, N gets a large Majorana mass on the UV brane because lepton number is gauged and broken there.

When the theory is written in 4D KK language, two different mass matrices appear:

**(a) Charged-lepton (Dirac) mass matrix**
$$
m_E \simeq 2 v k, f_L Y_E f_E ,
$$
relating left-handed doublets to right-handed singlets.
This is an ordinary **Dirac-type** matrix: it couples L↔E, and has one block per KK level.

**(b) Neutral-lepton (Majorana + Dirac) mass matrix**
In the basis ($\nu_L$, N, $N^(1)$, $N^c(1)$) it looks schematically like
$$
M_N =
\begin{pmatrix}
0      & Y_N v f_L f_N & Y_N v f_L & 0\\
Y_N v f_L f_N & M_{Maj}^{(0)} & M_{Maj}^{(01)} & 0\\
Y_N v f_L     & M_{Maj}^{(01)T} & M_{Maj}^{(1)} & M_{KK}\\
0              & 0               & M_{KK}        & 0
\end{pmatrix},
$$
where the $M_{Maj}$ entries come from UV-brane Majorana terms, and the $M_{KK}$ entries from the extra-dimension excitations.

This is the full object we want to diagonalize mixes light neutrinos, heavy sterile states, and their first KK partners.



1. We need the physical eigenstates. The particles that propagate in the $\mu \rightarrow e \gamma$ loop are the *mass eigenstates* of this combined system, not the “interaction-basis” fields.
2. We want to do it exactly. Rather than use a *mass-insertion* approximation (treating the small Dirac terms as perturbations), we build the entire matrix numerically and perform the full diagonalization to capture all mixings.
3. The neutral block is **symmetric, not Hermitian**, so we must use a *Takagi factorization* routine. The charged block is ordinary Dirac and needs SVD.


### Matrix Diagonalization: Hermitian, SVD, and Takagi

Different kinds of matrices can be “diagonalized” in different ways, depending on their symmetry and properties. 

Hermitian matrices represent self-adjoint transformations. They can be “rotated” by a unitary $U$ to reveal real eigenvalues $\Lambda$, corresponding to intrinsic directions and magnitudes of stretching. On such matrices we can use Eigendecomposition:

- **Eigendecomposition**  
  - **Applies to:** Hermitian (or real symmetric) matrices, where $A = A^\dagger$.  
  - **Form:** $A = U \Lambda U^\dagger$.  

For more general complex matrices, we use the fact that every linear transformation can be decomposed into a rotation $V^\dagger$, a scaling $\Sigma$, and another rotation $U$. The *singular values* in $\Sigma$ quantify how $A$ stretches space along orthogonal directions.

- **Singular Value Decomposition (SVD)**  
  - **Applies to:** Complex matrices (square or rectangular).  
  - **Form:** $A = U \Sigma V^\dagger$.  

However, SVD assumes you have two independent vector spaces (e.g. left and right handed components in a Dirac mass term). For Majorana terms, the field equals its charge conjugate, so those two spaces are identified.   So, SVD would overcount degrees of freedom and destroy the symmetry $M_N = M_N^T$ that the Majorana mass matrix must satisfy.

Instead, we use the fact that a symmetric complex matrix can be diagonalized by a single unitary matrix on both sides (with a transpose instead of a dagger). The diagonal entries of $\Sigma$ are the singular values of $A$.

- **Takagi Factorization**  
  - **Applies to:** Complex symmetric matrices, where $A = A^T$.  
  - **Form:** $A = U \Sigma U^T$.  
