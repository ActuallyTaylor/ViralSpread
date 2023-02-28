# TESTING THE ACCURACY OF A COMPUTER SIMULATION OF A VIRAL SPREAD
## A resarch paper submitted to the Owen J Roberts High School Science Department
This paper was wirtten for my Biotechnology class in High School. It is focused on the idea of simulating viral spread throughout a population. 

You can read the entire paper by looking at my [Final Paper](Final%20Paper.pdf)

## Summary of Trials
The following results were collected from a series of simulations all with different parameters. The biggest and most prominent simulation was completed over 18 hours with a population size of 500,000 units, 100 chunks, and a starting chunk density of 5,000. Other simulations included; a second 500,000 unit simulation, a 100,000 unit simulation, a 1,000 unit simulation and a 100 unit simulation. Prior to the completion of those simulations, there is a population size 100 simulation with health degradation removed.
Simulations were performed on two computers. Smaller, quicker simulations (100 - 1,000) were performed on a mac laptop with an Intel i9-9880H processor. These simulations were compiled using the Clang C++ compiler. Larger simulations (10,000 - 500,000) were run on an Intel i7-7700k processor. These simulations were compiled using the Visual C++ compiler. Differences in how the program performed between the two compilers were not noticed. The only difference between the two computers was the speed of the simulation. The i7-7700k was able to perform significantly better, resulting in faster simulations. This is because it was not thermally limited as it was in a desktop computer instead of a laptop.
