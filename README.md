# EnKF_EnOI_ES_EnKS
   
  A toy DA system that uses a 1D linear heat diffusion model to compare
  the following ensemble DA schemes: 
  - Ensemble Kalman Filter: EnKF
  - Ensemble Optimal Interpolation: EnOI
  - Ensemble Smoother: ES
  - Ensemble Kalman Smoother: EnKS
  
 The update scheme considers the observations all at once (i.e., batch style)
 and uses the transformation matrix (X5; Evensen, 2003). I also provided
 an EnKS function that assimilates the observations serially and uses 
 DART's style (2-step update, Anderson, 2003). 
 
 This a mere educational package. The coding style is not the best. 
 It's intended to become familiar with the different ensemble schemes, 
 their implementation and performance.
 
<object data="EnKF_EnKS.pdf" type="application/pdf" width="100%"> 
</object>

![smoother vs filter](EnKF_EnKS.pdf)
