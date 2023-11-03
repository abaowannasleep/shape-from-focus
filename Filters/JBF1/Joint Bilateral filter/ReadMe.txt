bfilter2 function perfroms two dimensional bilateral gaussian filtering.
The standard deviations of the bilateral filter are given by sigma1 and
sigma2, where the standard deviation of spatial-domain is given by sigma1
and the standard deviation intensity-domain is given by sigma2.
This function presents both bilateral filter and joint-bilateral filter.
If you use the same image as image1 and image2, it is the normal bilateral
filter; however, if you use different images in image1 and image2, you can
use it as joint-bilateral filter, where the intensity-domain (range weight)
calculations are performed using image2 and the spatial-domain (space weight) calculations are performed using image1.

Author: Mahmoud Afifi, York University. 

Usage:
  %Example1: normal bilateral filter using 5x5 kernel, spatial-sigma=6, and
  %intensity-sigma= 0.25:
  image=bfilter2(I1,I1,5,1.2,0.25);
  
  %Example2: joint-bilateral filter using 5x5 kernel, spatial-sigma=1.2,
  %and range-sigma= 0.25, the spatial-domain calculations are performed
  %using image (I1) and the intensity-domain calulcations (range weight)
  %are performed using image (I2):
  image=bfilter2(I1,I2,5,1.2,0.25);
  
  %Example3: use the default values for n, sigma1, and sigma2
  image=bfilter2(I1);

Input:
  -image1: the spatial-domain image
  -image2: the intensity-domain (range weight) image (use the same image
  for the normal bilateral filter. Use different images for joint-bilateral
  filter.
  (default, use the same image; i.e. image2=image1)
  -n: kernel (window) size [nxn], should be odd number (default=5)
  -sigma1: the standard deviation of spatial-domain (default=1.2)
  sigma2: the standard deviation of intensity-domain (default=0.25)

