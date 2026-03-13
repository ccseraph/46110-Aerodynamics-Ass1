# Sources used for expected test results

1. **Flow around a circular cylinder (Cambridge MDP)**  
   https://www-mdp.eng.cam.ac.uk/web/library/enginfo/aerothermal_dvd_only/aero/fprops/poten/node37.html  
   Used for: cylinder potential/stream function forms, polar velocity components, and zero-circulation pressure coefficient relation $C_p = 1 - 4\sin^2\theta$.

2. **Intermediate Fluid Mechanics (Oregon State) – Potential Flows**  
   https://open.oregonstate.education/intermediate-fluid-mechanics/chapter/potential-flows/  
   Used for: circulation-augmented cylinder formulas and stagnation-point relation on the surface.

3. **Kutta–Joukowski theorem (Wikipedia)**  
   https://en.wikipedia.org/wiki/Kutta%E2%80%93Joukowski_theorem  
   Used for: lift per unit span relation $L' = \rho U_\infty \Gamma$.

4. **Panel method note and MATLAB reference (course material)**  
   Literature/Panel-methods.pdf and matlab-files/PanelAirfoil.m  
   Used for: source-panel local influence formulas, self-induction values, panel geometry notation, and circulation extension workflow.

5. **Classic airfoil theory (ERAU open textbook)**  
   https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/classic-airfoil-theory/  
   Used for: Joukowsky transform and derivative forms (equivalent notation), trailing-edge preimage/cusp context, and conformal-mapping velocity relation.

6. **Complex Analysis textbook: Joukowsky map examples**  
   https://complexanalysis.org/web/sec_joukowski-airfoil.html  
   Used for: circle-to-ellipse mapping relations that support semi-axis ratio identities used to connect `a` and `c`.

7. **Potential Flow Theory course PDF + MATLAB snippet (course material)**  
   Literature/Potential-flow-theory_v2.pdf and Literature/Potential-flow-theory_v2.txt  
   Used for: symmetric displaced-circle Joukowsky setup (`a = c + s1`), Kutta circulation form `Gamma = 4*pi*(c+s1)*U_inf*sin(alpha)`, and MATLAB relation `c = a*sqrt((1-r)/(1+r))` with `r` semi-axis ratio.

8. **Joukowsky transform (Wikipedia, backup reference only)**  
   https://en.wikipedia.org/wiki/Joukowsky_transform  
   Used for: corroborative context on transform conventions and Kutta-condition narrative; kept as secondary source.