There are many copies of the muon event data. The different versions use different equations to calculate the quantity etastar.

VERSION ONE: v1

# Function calculates eta of the particle from its eta_star, the eta with which it leaves its vertex
def calc_etastar_from_eta_v1(eta, phi, x0, y0, z0, zstar=850.):
  # Propagate to station 2 (z = 850 cm)
  # Note: x0, y0, z0 in cm. Assume pT -> inf.
  if eta < 0:
    zstar *= -1
  cot = np.sinh(eta)
  delta_r = np.abs((zstar - z0)/cot)
  xstar = x0 + np.cos(phi) * delta_r
  ystar = y0 + np.sin(phi) * delta_r
  rstar = np.hypot(xstar, ystar)
  cotstar = zstar/rstar
  etastar = np.arcsinh(cotstar)
  return etastar
  
VERSION TWO: v2

def calc_etastar_from_eta_v2(invpt, eta, phi, x0, y0, z0, zstar=850.):
  # Propagate to station 2 (z = 850 cm)
  # Note: x0, y0, z0 in cm. Assume pT -> inf.
  if eta < 0:
    zstar *= -1
  # Assume a simplified magnetic field where it is 4T (or 3.811T)
  # inside the solenoid and 0T outside (boundary at z = 650 cm)
  zstar_4T = 650.
  if eta < 0:
    zstar_4T *= -1
  B = 3.811
  R = -1.0 / (0.003 * B * invpt)            # R = -pT/(0.003 q B)  [cm]
  cot = np.sinh(eta)
#   if np.abs(zstar_4T) < np.abs(zstar):
#     delta_r = np.abs(2 * R * np.sin((zstar_4T - z0)/(2 * R *cot)))  # with magfield
#     delta_r += np.abs((zstar - zstar_4T)/cot)                       # without magfield
#   else:
#     # Also need to check for the boundary at r, ignore for now
  delta_r = np.abs(2 * R * np.sin((zstar - z0)/(2 * R *cot)))     # with magfield -- indent this to return to previous
  xstar = x0 + np.cos(phi) * delta_r
  ystar = y0 + np.sin(phi) * delta_r
  rstar = np.hypot(xstar, ystar)
  cotstar = zstar/rstar
  etastar = np.arcsinh(cotstar)
  return etastar

VERSION THREE: v3

def calc_etastar_from_eta_v3(invpt, eta, phi, x0, y0, z0, zstar=850.):
  # Propagate to station 2 (z = 850 cm)
  # Note: x0, y0, z0 in cm. Assume pT -> inf.
  if eta < 0:
    zstar *= -1
  # Assume a simplified magnetic field where it is 4T (or 3.811T)
  # inside the solenoid and 0T outside (boundary at z = 650 cm)
  zstar_4T = 650.
  if eta < 0:
    zstar_4T *= -1
  B = 3.811
  R = -1.0 / (0.003 * B * invpt)            # R = -pT/(0.003 q B)  [cm]
  cot = np.sinh(eta)
  if np.abs(zstar_4T) < np.abs(zstar):
    arg_term = (zstar_4T - z0)/cot                      # with magfield
    sin_term = (2 * R) * np.sin(arg_term/(2 * R))       # with magfield
    cos_term = (2 * R) * (1 - np.cos(arg_term/(2 * R))) # with magfield
    arg_term = (zstar - zstar_4T)/cot                   # without magfield
    sin_term += arg_term                                # without magfield
    cos_term += 0.5 * np.square(arg_term)               # without magfield
  else:
    # Also need to check for the boundary at r, ignore for now
    arg_term = (zstar - z0)/cot                         # with magfield
    sin_term = (2 * R) * np.sin(arg_term/(2 * R))       # with magfield
    cos_term = (2 * R) * (1 - np.cos(arg_term/(2 * R))) # with magfield
  xstar = x0 + np.cos(phi) * sin_term - np.sin(phi) * cos_term
  ystar = y0 + np.sin(phi) * sin_term + np.cos(phi) * cos_term
  rstar = np.hypot(xstar, ystar)
  cotstar = zstar/rstar
  etastar = np.arcsinh(cotstar)
  return etastar

VERSION FOUR: v4

def calc_etastar_from_eta_v4(invpt, eta, phi, x0, y0, z0, zstar=850.):
  # Propagate to station 2 (z = 850 cm)
  # Note: x0, y0, z0 in cm. Assume pT -> inf.
  if eta < 0:
    zstar *= -1
  # Assume a simplified magnetic field where it is 4T (or 3.811T)
  # inside the solenoid and 0T outside (boundary at z = 650 cm)
  zstar_4T = 650.
  if eta < 0:
    zstar_4T *= -1
  B = 3.811
  R = -1.0 / (0.003 * B * invpt)            # R = -pT/(0.003 q B)  [cm]
  cot = np.sinh(eta)
  if np.abs(zstar_4T) < np.abs(zstar):
    arg_term_4T = np.abs((zstar_4T - z0)/cot)                 # with magfield
    sin_term_4T = (2 * R) * np.sin(arg_term_4T/(2 * R))       # with magfield
    cos_term_4T = (2 * R) * (1 - np.cos(arg_term_4T/(2 * R))) # with magfield
    arg_term_0T = np.abs((zstar - zstar_4T)/cot)              # without magfield
    sin_term_0T = arg_term_0T                                 # without magfield
    cos_term_0T = 0                                           # without magfield
  else:
    # Also need to check for the boundary at r, ignore for now
    arg_term_4T = np.abs((zstar - z0)/cot)                    # with magfield
    sin_term_4T = (2 * R) * np.sin(arg_term_4T/(2 * R))       # with magfield
    cos_term_4T = (2 * R) * (1 - np.cos(arg_term_4T/(2 * R))) # with magfield
    arg_term_0T = 0                                           # without magfield
    sin_term_0T = 0                                           # without magfield
    cos_term_0T = 0                                           # without magfield
  phistar_4T = phi + arg_term_4T/(2 * R)  # phi at the boundary where 4T -> 0T
  xstar = x0 + np.cos(phi) * sin_term_4T - np.sin(phi) * cos_term_4T + \
          np.cos(phistar_4T) * sin_term_0T - np.sin(phistar_4T) * cos_term_0T
  ystar = y0 + np.sin(phi) * sin_term_4T + np.cos(phi) * cos_term_4T + \
          np.sin(phistar_4T) * sin_term_0T + np.cos(phistar_4T) * cos_term_0T
  rstar = np.hypot(xstar, ystar)
  cotstar = zstar/rstar
  etastar = np.arcsinh(cotstar)
  return etastar

