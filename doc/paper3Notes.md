# Dynamically Removing False Features in Pyramidal Lucas-Kanade Registration

## Abstract

The natural LK algorithm might include falsely selected features. Traditional methods deal with this after the flow computation is completed, which incorporates additional cost. This paper proposes a handy tool to remove false features without degrading the algorithm's efficiency. This method uses a confidence preditor to directly evaluate an LK system's ill-posedness from the underlying data.

## Introduction

1. Introduce KL algorithm 

2. To tackle large displacement, implement flow computation by multi-level($l$) multi-stage($S$) refinement. An **L-level Gaussian pyramid**, solving a series of LK system
   $$
   A^{l,s}\delta X^{l,s+1} = b^{l,s}
   $$
   , where $X^{L,0}$ is initialized by a zero vector. The flow is updated by 
   $$
   X^{l,s+1} = X^{l,s} + \delta X^{l,s+1} \\
   X^{l-1,0} = \alpha X^{l,S}
   $$

3. **Why LK is fast and what did they select?**

   + recovers each flow vector individually, instead of computing the flow vectors of all image points jointly
   + image points that have strong local variation are selected for registration

   + **Erroneous of LK**: mang image points around flow boundaries have very large local variation, they are also selected as qualified features to register. However, flow recovery around flow boundaries is generally severely erroneous. (**TODO: why?**)
     + besides, when $AX=b$ is ill-posed, it may have 
       + multiple solutions (iff $A$ is rank deficient), aperture problem
       + no solution (iff the system is inconsistent)
     + "In a real application, the LK algorithm is applied selectively to texture image points. **This bypasses the aperture problem but inevitably mistakes flow boundary points in textured areas for reliable features.**"

4. Introduce several solutions to the previous problem

## Predicting LK Flow Confidence

A confidence predictor **P**, 1) has low cost, 2)decreases if the likelihood of the aperture problem or flow discontinuity increases.
$$
P &= \frac{f(\lambda_2^{2D})}{g(AX^0-b)}\\
$$
Then they claim the requirement and approximate the function $f$ and $g$, get this **STAR**(Spatial-Temporal vAriation Ratio)
$$
P = \frac{det(A^TA)}{||A^TA||_1\sum_NE_t^2}
$$
, where the cost is dominated by $4m$ multiplications and $4(m-1)$ additions, indicating that a flow vector can only be reliably recovered if the temporal sampling rate is "slower" than the spatial sampling rate.

## Dynamically detecting and discarding false features by STAR

**Why we can apply STAR?**

The LK algorithm formulates each feature's flow recovery independently, and hence the number of refinement stages is allowed to vary with the features.

### Algorithm

If $P$ is smaller than the threshold $T, the computation breaks and the image point $(x,y,t) is immediately removed from the registration list. if $P$ is larger than $T$, the LK flow computation proceeds to the next stage or level.