The equivalent system of ordinary differential equations is:

$$
N = S + I_1 + I_r + R + P,\quad\dot{N} = 0 \\
\dot{S} = \mu N - \beta S \frac{I_1 + \kappa I_r}{N} - \mu S \\
\dot{I_1} = \beta S \frac{I_1 + \kappa I_r}{N} - \left(\gamma+\mu\right) I_1 \\
\dot{I_r} = \kappa \beta \frac{I_1 + \kappa I_r}{N}P - \left(\frac{\gamma}{\kappa}  + \mu\right) I_r \\
\dot{R} = \gamma\left(I_1 + \frac{I_r}{\kappa} \right) - \rho R - \mu R \\
\dot{P} = \rho R - \kappa \beta \frac{I_1 + \kappa I_r}{N}P - \mu P
$$

Since population is constant, we can simply express the state variables as fractions of the total population:

$$
\dot{s} = \mu - \beta s \left(i_1 + \kappa i_r\right) - \mu s \\
\dot{i_1} = \beta s \left(i_1 + \kappa i_r\right) - \left(\gamma+\mu\right) i_1 \\
\dot{i_r} = \kappa \beta \left(i_1 + \kappa i_r\right)p - \left(\frac{\gamma}{\kappa}  + \mu\right) i_r \\
\dot{r} = \gamma\left(i_1 + \frac{i_r}{\kappa} \right) - \rho r - \mu r \\
\dot{p} = \rho r - \kappa \beta \left(i_1 + \kappa i_r\right)p - \mu p
$$

In the equilibrium condition, all $\dot{x}=0$, so:

$$
0 = \mu - \beta s \left(i_1 + \kappa i_r\right) - \mu s \\
0 = \beta s \left(i_1 + \kappa i_r\right) - \left(\gamma+\mu\right) i_1 \\
0 = \kappa \beta \left(i_1 + \kappa i_r\right)p - \left(\frac{\gamma}{\kappa}  + \mu\right) i_r \\
0 = \gamma\left(i_1 + \frac{i_r}{\kappa} \right) - \rho r - \mu r \\
0 = \rho r - \kappa \beta \left(i_1 + \kappa i_r\right)p - \mu p
$$

Adding the first two equations together:

$$
0 = \beta s \left(i_1 + \kappa i_r\right) - \left(\gamma+\mu\right) i_1 + \mu - \beta s \left(i_1 + \kappa i_r\right) - \mu s \\
0 =  - \left(\gamma+\mu\right) i_1 + \mu - \mu s \\
s = 1 - \frac{\gamma+\mu}{\mu} i_1 = s(i_1)
$$

and substituting that into second equation for $s$:

$$
0 = \beta s(i_1) \left(i_1 + \kappa i_r\right) - \left(\gamma+\mu\right) i_1 \\
\left(\gamma+\mu\right) i_1 - \beta s(i_1) i_1= \beta s(i_1) \kappa i_r \\
\left(\frac{\left(\gamma+\mu\right)}{\beta s(i_1)} - 1\right)\frac{i_1}{\kappa} = i_r(i_1)
$$

using that and re-arranging the fourth equation:

$$
0 = \gamma\left(i_1 + \frac{i_r}{\kappa} \right) - \rho r - \mu r \\
r(i_1) = \frac{\gamma}{\rho + \mu}\left(i_1 + \frac{i_r(i_1)}{\kappa} \right)
$$

now substituting into the third equation:

$$
0 = \kappa \beta \left(i_1 + \kappa i_r(i_1)\right)p - \left(\frac{\gamma}{\kappa}  + \mu\right) i_r(i_1) \\
\left(\frac{\gamma}{\kappa}  + \mu\right) \frac{i_r(i_1)}{\kappa \beta\left(i_1 + \kappa i_r(i_1)\right)} = p(i_1) \\
$$

substituting these into the constant population equation:

$$
1 = s + i_1 + i_r + r + p\\
1 = s(i_1) + i_1 + i_r(i_1) + r(i_1) + p(i_1)
$$

There are some options for simplification, but reduction still apparently leaves an equation requiring numerical solution.