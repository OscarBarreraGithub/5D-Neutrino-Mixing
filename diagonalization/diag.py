from typing import Tuple
import numpy as np

__all__ = ["SVD", "Takagi"]


def SVD(M: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Biunitary (SVD) diagonalization for Dirac-type mass matrices.

    Returns
    -------
    U_L, s, U_R  such that  U_L.conj().T @ M @ U_R = diag(s)  with s >= 0
    """
    M = np.asarray(M, dtype=np.complex128)
    U, s, Vh = np.linalg.svd(M, full_matrices=True)
    return U, s, Vh.conj().T

def _sqrtm_via_eig(M: np.ndarray) -> np.ndarray:
    """
    Matrix square root for *normal* matrices (our use: unitary/symmetric blocks).
    Falls back to eigendecomposition: V diag(sqrt(w)) V^{-1}.
    """
    w, V = np.linalg.eig(M)
    S = np.diag(np.sqrt(w))  # principal branch
    return V @ S @ np.linalg.inv(V)

def Takagi(
    A: np.ndarray,
    *,
    svd_order: bool = True,
    tol: float = 1e-12,
    symmetrize: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    r"""
    Autonne–Takagi decomposition of a complex symmetric matrix.

     See `Houde et al. Matrix decompositions in Quantum Optics:
     Takagi/Autonne, Bloch-Messiah/Euler, Iwasawa, and Williamson 


    Computes U, r such that A ≈ U @ diag(r) @ U.T with r ≥ 0 and U unitary.

    Parameters
    ----------
    A : (n, n) complex ndarray
        Input matrix (should satisfy A.T = A up to `tol`).
    svd_order : bool, default True
        If True, return singular values in descending order (NumPy SVD order).
        If False, return ascending order.
    tol : float, default 1e-12
        Tolerance for symmetry check and zero tests.
    symmetrize : bool, default False
        If True, project to the symmetric subspace via (A + A.T)/2 before factoring.
        If False, raise if A is not numerically symmetric.

    Returns
    -------
    r : (n,) float ndarray
        Nonnegative Takagi singular values (identical to the singular values of A).
    U : (n, n) complex ndarray
        Unitary Takagi factor so that A ≈ U @ diag(r) @ U.T.

    Notes
    -----
    - Real symmetric case uses an `eigh` shortcut with column phase adjustments.
    - Complex case uses the standard SVD congruence:
          A = U0 Σ V0* ,  Z := (V0* @ conj(U0)).T
          U = U0 sqrtm(Z)   (principal square root)
      so that A = U Σ U^T.
    - Requires SciPy for the most robust `sqrtm`; otherwise a stable eigen fallback
      is used (adequate here because Z is normal).
    """
    A = np.asarray(A)
    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        raise ValueError("A must be a square 2D array.")
    n = A.shape[0]

    # Optional least-squares projection onto the symmetric subspace
    if symmetrize:
        A = 0.5 * (A + A.T)

    # Enforce symmetry unless explicitly projected
    if np.linalg.norm(A - A.T) > tol:
        raise ValueError(
            "Input is not (numerically) symmetric; set symmetrize=True to project."
        )

    # Trivial zero case
    if np.allclose(A, 0, atol=tol, rtol=0):
        return np.zeros(n, dtype=float), np.eye(n, dtype=complex)

    # Fast real path via eigh
    if np.max(np.abs(A.imag)) <= tol:
        A = A.real  # be explicit for eigh
        evals, Q = np.linalg.eigh(A)  # ascending
        r = np.abs(evals)

        # Column phases: +1 for λ >= 0, +i for λ < 0
        signs = np.sign(evals)
        signs[signs == 0] = 1.0
        phases = np.sqrt(signs.astype(np.complex128))
        U = Q @ np.diag(phases)

        idx = np.argsort(r)[::-1] if svd_order else np.argsort(r)
        return r[idx], U[:, idx]

    # General complex-symmetric case
    U0, r, Vh = np.linalg.svd(A, full_matrices=True)  # r is descending by default
    Z = (Vh @ np.conjugate(U0)).T  # equals conj(U0.T @ V0); Z is symmetric & unitary

    # Try SciPy sqrtm; fall back to eig-based sqrt if SciPy is unavailable
    try:
        from scipy.linalg import sqrtm  # type: ignore
        S = sqrtm(Z)
    except Exception:
        S = _sqrtm_via_eig(Z)

    U = U0 @ S

    if not svd_order:  # return ascending instead of NumPy's default descending
        r = r[::-1]
        U = U[:, ::-1]

    return np.real_if_close(r), U
