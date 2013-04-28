Msreject uses the tbs feature of ms (see ms documentation for details) which allows you to simulate across a range of parameter values. Basically the program calculates summary statistics for simulations and only keeps simulations that approximate the observed data. It outputs the parameters from those accepted simulations, allowing the user to calculate the acceptance rate for a given model ( = # acceptances/ # sims run ) and estimate the posterior distribution of the parameters of interest (essentially just a histogram of the parameters returned by the program).

Currently the code assumes you are simulating two populations, but it can be easily modified to allow for more. The version here only calculates four summary statistics of the simulations: the number of shared (Ss) and fixed (Sf) sites between populations and the number of unique sites in each population (Sx1 and Sx2). Note that because msreject only reads ms output, this version can be used for ANY demographic model of two populations. One can use the number of acceptances for each of several potential models to calculate a Bayes Factor to choose among models (see ref [2] below for an example of this).

Please cite [1] or [2] when using this code.

Documentation

msreject is normally run by piping ms output directly to the program, like:

ms 10 1000 0t 5 -I 2 5 5 .5 | msreject -f output.file -d data.file -c cutoff -l 100


The required arguments above and additional optional arguments are listed below:

-d name of observed data file. must be a single line with the observed summary statistics in the same order as produced by msreject (in the version here, this should be as follows: mean Ss, variance Ss, mean Sf, variance Sf, mean Sx1, variance Sx1, mean Sx2, variance Sx2, all separated by whitespace)


-c cutoff. this is the percentage cutoff used for rejection. for example, if cutoff = 0.3, then only simulations in which EVERY summary statistic is within 30% of its observed value


-f output filename


-l the number of loci being simulated


-A when this argument is present msreject accepts all simulations regardless of the cutoff value or the summary statistics produced. useful for checking the range of summaries produced, or for using msreject to do the regression-ABC approach of Beaumont [3].

-m model number. this is useful for ABC model selection. the model number will be printed out with the other parameters

-V prints summary statistics and distances for every simulation run to
STDOUT. this will generate a LOT of text, so be ready for some big
files if you save this. useful for evaluating the behavior of certain
summary statistics and for doing regression ABC. this output will look
something like:


fabs: 2.24675 8.13341 0.82906 0.86788 0.186813 0.483558 0.434159 1.79404 1
simdata: 10 77.36 0.230769 1.38462 17.0769 130.474 10.1538 141.015 


where
fabs is the absolute % difference between observed and simulated data
for each simulation, and simdata are the simulated summary statistics
themselves (mean and variance across loci)


-M rejection done only on the mean of the summary statistics. msreject
normally rejects based on both the mean and variance of each summary
statistic.


-q removes singletons from the data before calculating summary statistics. 

not currently available


Download and Installation

You will need to have ms, msnsam, or some other simulation program that can produce ms-like output installed. You will also need to have the libsequence C++ library installed.

rejection.cc, the main file

SimDataTBS.hpp and SimDataTBS.cc ( found here ) are additional files needed to read in the tbs line from ms/msnsam

to compile the code, run something like:

g++ -I. -I.. -I/usr/local/include -Wall -W -ansi -pedantic -DNDEBUG -O2 -c rejection.cc SimDataTBS.cc 
g++ -DNDEBUG -O2 -o msreject rejection.o SimDataTBS.o -lsequence 


Example files

data.palu.txt, a file of observed data of mean and variance of SS, Sf, Sx1, and Sx2

temp.sh, a bash script to run rejection sampling on a cluster running sge, using msnsam, msreject, and priorgen


References

[1] Ross-Ibarra J, Wright SI, Foxe JP, Kawabe A, DeRose-Wilson L, et al. 2008. Patterns of Polymorphism and Demographic History in Natural Populations of Arabidopsis lyrata . PLoS ONE 3(6): e2411

[2] Ross-Ibarra, J., M. Tenaillon, J. Doebley, and B.S. Gaut. 2008. Historical divergence in the genus Zea. In Review.

[3] Beaumont, M.A., Zhang, W., and Balding, D.J. (2002) Approximate Bayesian Computation in Population Genetics. Genetics, 162: 2025-2035.

