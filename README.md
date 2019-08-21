Matlab library for flexible normalization that can be applied to
activation maps from deep Convolutional Neural Networks. It is based
on a mixture of Gaussian Scale Mixture Model described in:

L. Sanchez Giraldo and O. Schwartz. Integrating Flexible Normalization
into Mid-Level Representations of Deep Convolutional Neural Networks.
To Appear in Neural Computation, 2019.

The model differs from the flexible normalization introduced in: 
R. Coen-Cagli, P. Dayan, O. Schwartz. Statistical models of linear and
nonlinear contextual interactions in early visual processing. NIPS, 2009.

The model introduced here is designed to work with an arbitrary filter set, 
such as a learned layer of a deep convolutional neural network. In technical 
terms, in this model, the surround independent components have independent
mixers, which avoids the need for any symmetry assumption in the
parameters.

To use the code in this library, which is mainly contained in the "src"
directory, run the following command to add the relevant paths before
running any of your own scripts:  
>> run load_paths.m


The examples directory contains two scripts that illustrate the processes
of training and evaluation of a flexible normalization model on an
activation map from a Convolutional Neural Network or any related models.

The input array to be processed is a 4-dimensional array of size:
    	  	      N x C x H x W
N: number of instances
C: number of filters
H: spatial height
W: spatial width

Other important parameters are the size of the center neighborhod "n_cross_neigh,"
the distance between center and the spatial surround units "p_dist," and the spatial
configuration "spatial_conf."

For a demo script that trains a normalization model run:
"./examples/train_example_from_tensor_filelist.m"  


For a demo script that takes a trained model and gives back the normalized responses run:
"./examples/normalize_example_tensor_filelist.m"  


Both scripts take a list of .mat files. Each file should contain a 4-dimensional
array as the one described above, and C, H, W should be the same for all files.
See scripts for additional parameters.

Luis Gonzalo Sanchez Giraldo
Odelia Schwartz 
August, 2019

Last Rev: Aug-16-2019
