import numpy as np

from yukawa.neutrino import compute_neutrino_yukawas


def test_compute_neutrino_yukawas_matches_eq6_scalar_inversion():
    """Scalar universal-limit inversion should match the Eq. (6) algebra."""
    f_L = 0.015976144656324655
    f_N = 0.4795831681596669
    f_N_UV = 0.00012321551579024553
    M_N = 1.22e18
    masses_eV = np.array([0.002, 0.009, 0.05])
    v = 174.0
    k = 1.2209e19

    result = compute_neutrino_yukawas(
        f_L=f_L,
        f_N=f_N,
        f_N_UV=f_N_UV,
        M_N=M_N,
        neutrino_masses_eV=masses_eV,
        k=k,
        v=v,
    )

    expected_y = np.sqrt(
        masses_eV * 1e-9 * (f_N_UV**2) * M_N / (2.0 * k**2 * v**2 * f_L**2 * f_N**2)
    )
    expected_y_bar = 2.0 * k * expected_y

    assert np.allclose(result["Y_N"], expected_y, rtol=1e-12, atol=0.0)
    assert np.allclose(result["Y_N_bar"], expected_y_bar, rtol=1e-12, atol=0.0)
