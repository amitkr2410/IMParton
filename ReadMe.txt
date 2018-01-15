IMParton v1.0
-- IMParton = I'm Parton

1. Description
This package gives parton distribution functions (PDFs) of the proton
starting from low Q^2 ~ 0.07 GeV^2, which are based on the analysis to
deep inelastic scattering data applying DGLAP equations with nonlinear
corrections. Refs. Rong Wang, Xu-Rong Chen, Chinese Physics C Vol. 41,
No. 5 (2017) 053103 [arXiv:1609.01831v4]; Xurong Chen et al, Int. J. 
Mod. Phys. E 23 (2014) 1450057 [arXiv:1306.1872v2]; Wei Zhu et al, 
Eur. Phys. J. Plus 131 (2016) 6 [arXiv:1404.0759v2].

2. Usage
The library consists of a C++ class named IMParton (read ./IMParton.h
and ./IMParton.cpp for details). IMParton has a method
IMParton::getPDF(Iparton, X, Q2), which is suggested to be called
in users' programs. Iparton set as -4, -3, -2, -1, 0, 1, 2, 3, 4
corresponds to getting cbar, sbar, dbar, ubar, gluon, u, d, s, c
quark/gluon number density distributions respectively.
The other important method of IMParton is setDataSet(int).
setDataSet(1) corresponds to use the data set A, and setDataSet(2)
corresponds to use the data set B.

We also provide methods getXUV(X, Q2), getXDV(X, Q2), getXUSea(X, Q2),
getXDSea(X, Q2), getXSSea(X, Q2), getXCSea(X, Q2) and getXGluon(X, Q2)
to get the up valence quark, down valence quark, up sea quark,
down sea quark, strange sea quark, charm sea quark and gluon distribution
functions respectively. It should be noted that the returned values
are the distributions multiplied by x, for these methods.
(up valence quark = x(u - ubar)
down valence quark = x(d - dbar)
up sea quark = x*2*ubar
down sea quark = x*2*dbar
strange sea quark = x*2*sbar
charm sea quark = x*2*cbar)

./test.cpp gives an example to get up and down valence quark distributions
at low Q^2, and sbar, gluon, dbar/ubar distributions at high Q^2.
./test.cpp can be modified as users' wants.
To run the example,
>tar -zvxf IMParton.tar.gz
>cd IMParton
>make
>./testIMParton

3. Method
We provide two data sets of PDFs obtained from the global analysis to
experimental data. Set A is from the three valence quarks non-perturbative
input, and set B is from the non-perturbative input of three valence quarks
adding flavor-asymmetric sea quarks components. This package contains two
table grids, which are ./grid_1_1_SetA.dat and ./grid_1_1_SetB.dat.
The table grids are generated in the kinematic range of 10^-6 < x < 1
and 0.125 < Q^2 < 2.68435e+08 (GeV^2). The number of the grid points is 32
for the variable Q^2 and 61 for the variable x. For PDF values in small x
region (x<10^-4), we use function form A*x^B to do interpolation. In large
x region (x>0.5), we use function form A*(1-x)^B to do the interpolation.
In middle x region (10^-4<x<0.5), we use quadratic interpolation. The PDF
values outside of the grid range are given using extrapolation method.
The sophisticated extrapolation method is expected to be effective.

4. Questions
If you have detailed questions concerning these distributions,
or if you find problems/bugs using this package, direct inquires to
wangrong11@mails.ucas.ac.cn (rwangcn8@gmail.com), xchen@impcas.ac.cn,
or wzhu@phy.ecnu.edu.cn.
