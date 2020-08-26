This collection of codes is based on "Estimating overidentified, non-recursive, time varying coefficients structural VARs"
by Fabio Canova and Fernando J. Pérez Forero (Quantitative Economics, Volume 6, Issue 2 (July 2015)).

Modified to run on GNU Octave 5.2.0 by Huachen Li

File list:
carter_kohn1beta.m
data.mat
draw_beta_conditional.m
MAIN_EST.m
PLOTIRF.m
PLOTSV.m
PRIOR_EST.m
Readme.txt
rows.m

Estimation instructions:
1. Unzip and copy all files to the same folder as those by Fabio Canova and Fernando J. Pérez Forero.
2. Run PRIOR_EST.
3. Run MAIN_EST.
4. Run PLOTSV and PLOTIRF for SV and IRF results. 
Note: the program is optimized to run on GUN Octave 5.2.0. 
If error "maximum_recursion_depth_exceeded" occurs, delete quantile.m from the root folder and reload workspace.