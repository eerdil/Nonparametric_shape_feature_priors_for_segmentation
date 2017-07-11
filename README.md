The implementation was used in:

[1] Ertunc Erdil, M. Usman Ghani, Lavdie Rada, A. Ozgur Argunsah, Devrim Unay, Tolga Tasdizen, Mujdat Cetin "Nonparametric Joint Shape and Feature Priors for Image Segmentation", IEEE Transactions on Image Processing 2017

Any papers using this code should cite [1] accordingly.

The software has been tested under Matlab R2015a and Microsoft Visual C++ 2010.

After unpacking the file and installing the required libraries,
start Matlab and run the following int the root directory:

>> compile_mex_files

If no errors are reported, you can then run "nonparametric_shape_feature_priors.m" in the root directory.

If errors are reported, you probably have problem with your mex compiler. Please make
sure that your mex compiler works properly. 

If you still have problems, you can email me at ertuncerdil@sabanciuniv.edu
I will do my best to help.

As is, the code produces the results given as the first experimental setting with MNIST data set in the paper.
You can create a new folder similar to them to test the algorithm on various data sets. 
Note that, you may also need to change some parameters in "nonparametric_shape_feature_priors.m" where I tried to comment heavily.

Please also report any bug to ertuncerdil@sabanciuniv.edu
