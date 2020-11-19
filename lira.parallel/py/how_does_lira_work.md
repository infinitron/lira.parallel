# The basic skeleton of LIRA workflow is as follows:

1. Take an X-ray image and create a psf of the core 
2. Create a baseline model from the deconvolved image of the observation i.e., just a core and a background. The total number of counts remain the same in both the observed and the baseline images. This will be the *null model*.
3. Simulate a few realizations of the null model.
4. Now for each image i.e., the observed image and the relplicates, LIRA fits an additional component to the null model to reproduce the image.
5. For each region, take the fraction of the additional component: zeta = t_add/(t_null+t_add)
6. If the counts in the image were just because of a random fluctuation the distrubution of zeta for the replicated images and the observed images would be the same. A deviation would indicate an additional component.
7. Then compute an upper bound on the p value to  obtain the significance.