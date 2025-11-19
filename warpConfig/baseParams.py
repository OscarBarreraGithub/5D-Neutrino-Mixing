import math

# Electroweak vev (GeV)
V_EWSB = 174.0

MPL = 1.2209e19      # Planck mass (GeV)
MBARPL = 2.435e18    # reduced Planck mass (GeV)


# Default values. Can be adjusted after importing 
# later with e.g. params_local = get_warp_params(Lambda_IR=3000)
DEFAULT_K = MPL
DEFAULT_LAMBDA_IR = 3000.0  # 3 TeV

def get_warp_params(k=DEFAULT_K, Lambda_IR=DEFAULT_LAMBDA_IR):
    """
    Calculate and return warp geometry parameters.
    
    Parameters:
    -----------
    k : float
        AdS curvature (mass dim +1). Default is MPL.
    Lambda_IR : float
        IR scale (KK mode mass scale) = k * exp(-pi k rc). Default is 3000 GeV (3 TeV).
        
    Returns:
    --------
    dict
        Dictionary containing:
        - k: curvature
        - Lambda_IR: IR scale
        - epsilon: warp factor (Lambda_IR / k)
        - rc: radius of the orbifold
        - z_h: UV brane position (1/k)
        - z_v: IR brane position (1/Lambda_IR)
        - warp_log: pi * k * rc
    """
    # Derive epsilon from Lambda_IR and k
    # Lambda_IR = k * epsilon  =>  epsilon = Lambda_IR / k
    epsilon = Lambda_IR / k
    
    if epsilon <= 0:
        raise ValueError("epsilon must be > 0. Ensure Lambda_IR and k have the same sign.")
        
    # Derive rc from epsilon
    # epsilon = exp(-pi k rc)  =>  ln(epsilon) = -pi k rc  =>  rc = -ln(epsilon) / (pi k)
    warp_log = -math.log(epsilon)
    rc = warp_log / (math.pi * k)

    return {
        "k": k, 
        "Lambda_IR": Lambda_IR, # IR "KK scale" knob: Lambda = k * exp(-pi k rc)
        "epsilon": epsilon, # epsilon = exp(-pi k rc) = Lambda / k
        "rc": rc, # radius of the orbifold (mass dim -1)
        "z_h": 1.0 / k, # UV brane position in conformal coord (1/k)
        "z_v": 1.0 / Lambda_IR, # IR brane position in conformal coord (e^{pi k rc}/k) = 1/(k*eps)
        "warp_log": warp_log # pi*k*rc = ln(k/Lambda)
    }


