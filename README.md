# Coalescent_Internal_Branches
Software for calculating the distribution of internal branch lengths in the Kingman Coalescent, as well as related summary statistics.

## Recommended Installation instructions

###Option 1
If the program cmake is available on the computer where you are installing, both the internal branches program and the modified coalescent simulator can be installed in the following way. First clone the repository or download the tarball with the necessary files.

In cloned repository:
   
   cmake .
   make   	     

Using tarball:
      
      tar -xvf Coalescent_Internal_Branches.tar.gz
      cd Coalescent_Internal_Branches_
      cmake .
      make

###Option 2
If cmake is unavailable, the programs can be individually compiled using the following commands

   g++ -o FTEC_branchlength -I . FTEC_branchlength.cpp
   g++ -o internal_branches -I . internal_branches.cpp
	  

### Commands 

#### Start Probabilties
Outputs a matrix with the probability that a coalescent branch with a given size arises at a particular coalescent event. Rows represent coalescent events, with the first row giving the initial state before any coalescent events have occurred (100% probability that all branches have size 1). Each column represents the corresponding branch size. For example, column 2 has a 1 in row 2 and zeros everywhere else, meaning there is a 100% probability of a branch with size 2 arising at the first coalescent event. 

Example command to output matrix for a sample size of 10, for branches up to size 7:

`./internal_branches -o start_probabilities -n 10 -s 7`

#### End Probabilties
Outputs a matrix with the probability that a coalescent branch ends at a particular coalescent event, conditional on the event at which it arose. Rows correspond to starting event, columns correspond to ending event. First row corresponds to external branches existing before first coalescent event, and first column is all zeros to represent the probability of a branch ending before first coalescent event. These matrices are independent of branch size. 

Example command to output matrix for a sample of size 10:
 
`./internal_branches -o end_probabilities -n 10 -s 2`

#### Branch Numbers
Outputs a table with the probabilities of observing a given number of branches with a specific size in a genealogy. 

Example command to calculate probabilities for branches with sizes 4 and 5 in a sample of size 25:
 
`./internal_branches -o branch_numbers -s 4,5 -n 25`
  
##### Output from the proceeding commands can be stored and used as input for the following commands
  
#### Length Distribution
Given a series of "breaks", this will output the cumulative probability of branches being shorter than the given values. Can use a file to enter breaks, where there is a single line with comma separated values. Users input the number of random branches to base the distribution on, with more branches increasing the accuracy but also the run time. 

Example command to output probability branches of size 5 are less than lengths 0.1,0.5, and 1 in a sample size of 25
 
`./internal_branches -o length_distribution -s 5 --breaks 0.1,0.5,1 -n 25 -x 2500000`
  
Example command to output probability branches of size 5 are less than lengths 0.1,0.5, and 1 in a sample size of 1000 using precomputed start and end probability files
  
`./internal_branches -o length_distribution -s 5 --breaks 0.1,0.5,1 -n 1000 --starts starts_n1000.txt --ends ends_n1000.txt`
  
#### Genealogy Summary Statistics
Output total number of branches, their summed length, and the inter-branch variance for a specified number of genealogies. Columns after these statistics contain the the individual branch lengths.

Example command to output statistics for branches with size 5 in 50 genealogies from a sample of size 100
 
`./internal_branches -o genealogy_portions -s 5 -n 100 -x 50`

####  R2 Probabilities
Will give you the expected probability that two observed variants of the same size arose on the same ancestral branch of a genealogy (in the absence of recombination). For this command -x inputs the number of genealogies to sample rather than the number of branches. The number of branches in each genealogy is determined by the calculations for branch_numbers.   
  
Example command to output trend in probability between sample sizes of 50 and 100 for branches of size 5
 
`./internal_branches -o r2_prob -s 5 -n 50,60,70,80,90,100 -x 50000`


Packages required for compiling the raw source code

Tclap - http://tclap.sourceforge.net/
Boost Random - http://www.boost.org/doc/libs/1_40_0/libs/random/
Eigen - http://eigen.tuxfamily.org/index.php?title=Main_Page

####Helpful Files
StartProbabilities_n1000.txt.gz includes the starting probabilities for branches up to size 15 for all sample sizes up to 1000, and can be used to speed up other functions when uncompressed and included with --starts command.