Instructions to run the function dtree in MATLAB.
The name of the function is dtree.
In the MATLAB command window enter the following commands.


1.type  " dtree pendigits_training.txt pendigits_testing.txt optimized '50' " in the command window in Matlab
2.Press enter



The first argument is the training file name which can be changed as and when needed.
The second argument is the testing file name which can be changed as and when needed.
The third argument is the option to choose from optimized trees, randomized trees or forests.
Th fourth argument is the pruning threshold.

The other possible running combinatins are:

dtree pendigits_training.txt pendigits_testing.txt randomized '50'
dtree pendigits_training.txt pendigits_testing.txt forest3 '50'