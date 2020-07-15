### Dataset and code released for
# Xun Xu, Loong-Fah Cheong and Zhuwen Li, Motion segmentation by exploiting complementary geometric models, CVPR 2018
#
# By: Xun Xu, Sep 2018
# Contact: alex.xun.xu@gmail.com

Dataset:

The KT3DMoSeg dataset [1] was created upon the KITTI benchmark [2] by manually selecting 22 sequences and labelling each individual foreground object. We select sequence with more significant camera translation so camera mounted on moving cars are preferred. We are interested in the interplay of multiple motions, so clips with more than 3 motions are also chosen, as long as these moving objects contain enough features for forming motion hypotheses. 22 short clips, each with 10-20 frames, are chosen for evaluation. We extract dense trajectories from each sequence using [3] and prune out trajectories shorter than 5 frames.

Raw Sequence:

We provide the raw sequence under ./OriginalSequence/

Data Format:

All meta data are saved as Matlab mat files under ./Data/ . All information are fields of struct varaible "Data". All dense trajectories are indexed by "Data.yAll" and all sparse trajectories are indexed by "Data.ySparse". The ground-truth label for all sparse trajectories is indexed by "Data.GtLabel". Motion segmentation is only evaluated on sparse trajectories. The visible mask is indexed by "Data.visibleSparse/visibleAll". The visibleSparse/visibleAll is a boolean matrix with NxL dimension where N is the number of trajectories and L is the total number of frames. Additional original video sequences can be accessed at https://www.dropbox.com/sh/4u1p3xwe6v48kww/AACLNzNSWg_YYZ6dsIXtLRoTa?dl=0

Demo Code:

A demo code to sample hypotheses (Affine, Homography and Fundamental Matrix), Compute Ordered Residual Kernel (ORK) [4] and single-view or multi-view spectral clustering is provided. The code is based on sampling hypotheses from all sparse trajectories and test/evaluate on sparsely labelled trajectories.
To run this demo, follow the steps:
(1) Run hypotheses generation code from each model under ./script/Hypo/
(2) Run ORK kernel computing code under ./script/Kernel/
(3) Run single-view (script_X_MoSeg_RandSamp.m, where X is Affine/Homography/Fundamental) or multi-view (script_AHF_X_RandSamp.m, where X is KerAdd/CoReg/Subset) motion segmentation under ./script/MoSeg/
(4) Check motion segmentation results by running the code under ./script/CheckPerf/
Due to randomness in hypotheses sampling, the performance may vary in each run. An optimal alpha is around 10.


Please cite [1] if you wish to use this dataset and demo code.

The libraries and code to support this demo include:

MATLAB Functions for Multiple View Geometry
David Capel, Andrew Fitzgibbon, Peter Kovesi, Tomas Werner, Yoni Wexler, and Andrew Zisserman
http://www.robots.ox.ac.uk/~vgg/hzbook/code/

Accelerated Hypothesis Generation for Multi-Structure Robust Fitting
T.-J. Chin
https://cs.adelaide.edu.au/~tjchin/doku.php?id=publications_code

Sparse Subspace Clustering (SSC)
Ehsan Elhamifar
http://www.ccs.neu.edu/home/eelhami/codes.htm

[1] X. Xu, L.F. Cheong, and Z. Li. Motion segmentation by exploiting 
complementary geometric models, In CVPR 2018.
[2] A. Geiger, P. Lenz, C. Stiller, and R. Urtasun. Vision meets robotics: The kitti dataset. International Journal of Robotics Research, 2013.
[3] N. Sundaram, T. Brox, and K. Keutzer. Dense point trajectories by GPU-accelerated large displacement optical flow. In ECCV, 2010.
[4] T. Chin, H. Wang and D. Suter. The ordered residual kernel for robust motion subspace clustering. In NIPS, 2009.