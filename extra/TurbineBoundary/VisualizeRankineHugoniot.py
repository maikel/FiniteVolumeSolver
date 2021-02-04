import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math
import scipy.optimize
import numpy

def set_size(width, fraction=1, subplots=(1, 1)):
    """Set figure dimensions to avoid scaling in LaTeX.

    Parameters
    ----------
    width: float or string
            Document width in points, or string of predined document type
    fraction: float, optional
            Fraction of the width which you wish the figure to occupy
    subplots: array-like, optional
            The number of rows and columns of subplots.
    Returns
    -------
    fig_dim: tuple
            Dimensions of figure in inches
    """
    if width == 'thesis':
        width_pt = 426.79135
    elif width == 'beamer':
        width_pt = 307.28987
    else:
        width_pt = width

    # Width of figure (in pts)
    fig_width_pt = width_pt * fraction
    # Convert from pt to inches
    inches_per_pt = 1 / 72.27

    # Golden ratio to set aesthetic figure height
    # https://disq.us/p/2940ij3
    golden_ratio = (5**.5 - 1) / 2

    # Figure width in inches
    fig_width_in = fig_width_pt * inches_per_pt
    # Figure height in inches
    fig_height_in = fig_width_in * golden_ratio * (subplots[0] / subplots[1])

    return (fig_width_in, fig_height_in)

#plt.style.use('seaborn')
tex_fonts = {
    # Use LaTeX to write all text
    "text.usetex": True,
    "font.family": "serif",
    # Use 10pt font in plots, to match 10pt font in document
    "axes.labelsize": 9,
    "axes.titlesize": 9,
    "axes.labelsize": 9,
    "font.size": 9,
    # Make the legend/label fonts a little smaller
    "legend.fontsize": 9,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7
}

plt.style.use('seaborn')
plt.rcParams.update(tex_fonts)
#plt.rcParams.update({'axes.grid' : False})

class PerfectGas:
  pass

equation = PerfectGas()
equation.gamma = 1.24
equation.gm1 = equation.gamma - 1.0
equation.gp1 = equation.gamma + 1.0
equation.oogm1 = 1.0 / equation.gm1
equation.gm1_over_gp1 = equation.gm1 / equation.gp1
equation.gm3 = equation.gamma - 3.0
equation.R = 1.0

class State:
  pass

left = State()
left.rho = 1.0
left.p = 1.0
left.a = math.sqrt(left.p / left.rho * equation.gamma)
left.u = 0.0
left.rhou = left.rho * left.u
left.rhoE = left.p / equation.gm1 + 0.5 * left.rhou * left.rho

rhos = numpy.linspace(0.001, 4.0, 2000)
mdot = 5.0
u_mdots = mdot / rhos

mdot2 = 0.025
u_mdots2 = mdot2 / rhos

# def u(rho):
#   return uL + 2.0 * math.sqrt(gamma) / (gm1 * math.sqrt(rhoL)) * (math.sqrt(pL) - numpy.sqrt(rho / rhoL)**gm1)

# def solver_star_state(mdot):
#   f = lambda rho: rho * (uL + 2.0 * math.sqrt(gamma) / (gm1 * math.sqrt(rhoL)) * (math.sqrt(pL) - math.sqrt(rho / rhoL)**gm1)) - mdot
#   result = scipy.optimize.root(f, rhoL)
#   if not result.success:
#     print(result)
#     raise "Failed"
#   rho_star = result.x[0]
#   u_star = mdot / rho_star
#   p_star = (rho_star / rhoL)**gamma
#   return (rho_star, u_star, p_star)

# (rho_star, u_star, p_star) = solver_star_state(0.6)
# print((rho_star, u_star, p_star))
# print(rho_star * u_star)

def f_shock(p, w, eq):
  assert(p > 0.0)
  A = 2.0 / (eq.gp1 * w.rho)
  B = w.p * eq.gm1_over_gp1
  ooQ = math.sqrt(A / (p + B))
  return (p - w.p) * ooQ

def df_shock(p, w, eq):
  assert(p > 0.0)
  A = 2.0 / (eq.gp1 * w.rho)
  B = w.p * eq.gm1_over_gp1
  Q = A / (p + B)
  DQ = -A / (p + B)**2
  sqrQ = math.sqrt(A / (p + B))
  DsqrQ = 0.5 / sqrQ * DQ
  return sqrQ + (p - w.p) * DsqrQ

def u_shock(p, w, eq):
  return w.u - f_shock(p, w, eq)

def Du_shock(p, w, eq):
  return df_shock(p, w, eq)

def rho_shock(p, w, eq):
   return w.rho * (equation.gm1_over_gp1 + p / w.p) / (equation.gm1_over_gp1 * p / w.p + 1.0)

def delta_rhou_shock(p, rhou_star, w, eq):
  u_star = u_shock(p, w, eq)
  rho_star = rho_shock(p, w, eq)
  return rho_star * u_star - rhou_star

def Ddelta_rhou_shock(p, rhou_star, w, eq):
  u_star = u_shock(p, w, eq)
  du_star = Du_shock(p, w, eq)
  g = (equation.gm1_over_gp1 + p / w.p)
  Dg = 1.0 / w.p
  h = (equation.gm1_over_gp1 * p / w.p + 1.0)
  Dh = equation.gm1_over_gp1 / w.p
  rho_star = w.rho * g / h
  drho_star = w.rho * (Dg * h - g * Dh) / h**2
  return rho_star * du_star + drho_star * u_star

def f_rarefaction(p, w, eq):
  return 2.0 * w.a * eq.oogm1 * ((p / w.p)**(eq.gm1 / 2.0 / eq.gamma) - 1.0)

def Df_rarefaction(p, w, eq):
  exponent = eq.gm1 / 2.0 / eq.gamma
  return 2.0 * w.a * eq.oogm1 * exponent * (p / w.p)**(exponent - 1.0) 

def u_rarefaction(p, w, eq):
  return w.u - f_rarefaction(p, w, eq)

def Du_rarefaction(p, w, eq):
  return Df_rarefaction(p, w, eq)

def rho_rarefaction(p, w, eq):
  return w.rho * (p / w.p)**(1.0 / eq.gamma)

def delta_rhou_rarefaction(p, rhou_star, w, eq):
  u_star = u_rarefaction(p, w, eq)
  rho_star = rho_rarefaction(p, w, eq)
  return rho_star * u_star - rhou_star

def Ddelta_rhou_rarefaction(p, rhou_star, w, eq):
  u_star = u_rarefaction(p, w, eq)
  du_star =  Du_rarefaction(p, w, eq)
  rho_star = w.rho * (p / w.p)**(1.0 / eq.gamma)
  drho_star = w.rho * (p / w.p)**(1.0 / eq.gamma - 1.0) / eq.gamma
  return rho_star * du_star + drho_star * u_star

def f(p, w, eq):
  if p <= w.p:
    return f_rarefaction(p, w, eq)
  else:
    return f_shock(p, w, eq)

def delta_rhou(p, rhou_star, w, eq):
  if p <= w.p:
    return delta_rhou_rarefaction(p, rhou_star, w, eq)
  else:
    return delta_rhou_shock(p, rhou_star, w, eq)

def Ddelta_rhou(p, rhou_star, w, eq):
  assert(p > 0)
  if p <= w.p:
    return Ddelta_rhou_rarefaction(p, rhou_star, w, eq)
  else:
    return Ddelta_rhou_shock(p, rhou_star, w, eq)


upper = 4.0
lower = 0.375
upper = 0.5 * (lower + upper)
upper = 0.5 * (lower + upper)
upper = 0.5 * (lower + upper)
upper = 0.5 * (lower + upper)
upper = 0.5 * (lower + upper)
upper = 0.5 * (lower + upper)
upper = 0.5 * (lower + upper)
upper = 0.5 * (lower + upper)
upper = 0.5 * (lower + upper)
lower = 0.5 * (lower + upper)
lower = 0.5 * (lower + upper)
upper = 0.5 * (lower + upper)
lower = 0.5 * (lower + upper)
current = 0.5 * (lower + upper)
rhou = current

f_opt = lambda p: delta_rhou(p, rhou, left, equation)

#p = 1.0598 * left.p

# Df_opt = lambda p: Ddelta_rhou(p, rhou, left, equation)
p0 = left.p

print(f_opt(p0))
if f_opt(p0) < 0.0:
  upper = p0
  lower = 0.5 * p0
  while (f_opt(lower) < 0.0):
    lower = 0.5 * lower
else:  
  lower = p0
  upper = 2.0 * p0
  while (f_opt(upper) > 0.0):
    upper = 2.0 * upper

print(f_opt(lower))
print(f_opt(upper))

result = scipy.optimize.root_scalar(f_opt, bracket=(lower, upper), method='bisect')
print(result)
p = result.root
print('p: {}, delta(rhou): {}'.format(p, f_opt(p)))
if p <= left.p:
  u = u_rarefaction(p, left, equation)
  rho = rho_rarefaction(p, left, equation)
else:
  u = u_shock(p, left, equation) 
  rho = rho_shock(p, left, equation)
print(p, u, rho)
# print(p / rho)
print(rho * u)
a = math.sqrt(equation.gamma * p / rho)
print(u - a)

right = State()
right.rho = rho
right.u = u
right.p = p
right.a = math.sqrt(equation.gamma * right.p / right.rho)
right.rhou = right.rho * right.u
right.rhoE = right.p / equation.gm1 + 0.5 * right.rhou * right.rho

result = scipy.optimize.root(lambda p: f(p, left, equation) + f(p, right, equation) + right.u - left.u, 0.5 * (left.p + right.p))
print(result)
u = 0.5 * (left.u + right.u) + 0.5 * (f(p, right, equation) - f(p, left, equation))
print(u)
# f, ax = plt.subplots(nrows=1, ncols=1, figsize=set_size('thesis'))
# ax.plot(rhos, u(rhos), label='Locus Rarefaction')
# ax.plot(rhos, u_mdots, label='$\\rho u = {}$'.format(mdot))
# ax.plot(rhos, u_mdots2, label='$\\rho u = {}$'.format(mdot2))
# ax.scatter(rhoL, uL)
# ax.annotate('({}, {})'.format(rhoL, uL), (rhoL, uL), textcoords='offset points', xytext=(0, 10), ha='center')
# ax.annotate('$(\\rho_L, u_L)$', (rhoL, uL), textcoords='offset points', xytext=(0, 10), ha='center')
# ax.set(xlabel='$\\rho$', ylabel='$u$')
# ax.legend()
# f.savefig('locus_rarefaction.png')
# f.clear()