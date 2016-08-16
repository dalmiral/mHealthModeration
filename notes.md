1. My understanding is that we should not refer to a "conditional mean model for E[Y_t+k | A_t , H_t ]" when discussing the weighted estimator...since the weighted estimator is not estimating this conditional mean model. 
Do I have this correct?  

​2. About routine regres. assumption R2: This assumption says that this true:
(R2)         E[Y_t+k | A_t , H_t ] = S_kt α_k + A_t S_kt β_k    

Can this assumption instead be written as having three parts: 
a) S_kt α_k is a correct working model for E[Y_t+k | A_t = 0, H_t ], and 
b) S_tk contains all moderators in H_t
c) Treatment effect model S_kt β_k  is correct

i.e., "a) and b) and c)" imply R2?  

4. For weighting or centering, there currently is no simulation showing variance reductions when "better" working models are used (e.g., putting in some predictive covariates in working model vs an intercept only model). 
I think we should have this simulation. Do you agree?

5. Proposition 3.2 says "...and at least one of C2 and C3."
This twisted me up. 
Should it be "...and either C2 or C3" ?

6. Can I change the text describing Simulation 6.2 so that only the 3 conditions being shown in the table are the only ones we discuss/present in the lead up?  I had a really hard time getting 6.2 straight .  And, could I be in touch with more questions about it once I start getting it straighter in my head?

7. Can you take a look at my changes to text in 6.1 and tell me what you don't agree with?

- run case where there are no trt effects

11/23

- mention up front that simulations consider only randomized studies for clarity
- don't bother implementing a working variance model
- check indexing of variables
- defer SD, ASE appendix - mention in caption std dev and avg std dev are equal up to two decimal places we relegate this to the appendix

- weight stabilization - trt prob stratified on binary variable
- P_n rho_t versus P_n A_t
- P_n sum_t rho_t(H_t) / T

- baseline candidate moderators for stabilization? don't bother

- Heagerty simulation - bias variance trade-off

- make generative, analysis models more similar

- AR(1) correction of SEs - known, true, not estimated for simplicity

- not working model - conditional mean model correct
- working model - 

- trt discontinuity study - implemented *a* dynamic trt regime

- random effects, not multilevel model

- currently we are involved in multiple MRTs, methods are useful for these studies

- in example code - do we need to fix the scale?

- check centered estimator

wgt stabilization scenario
- show problems when numerator is a function of time
- same generative model as omitted moderator?

- 4pm tues


sim(25, 25, 25, scenario = "omit")
sim(25, 25, 25, scenario = "stable", beta0 = c(-0.8, 0, 0, 0))
sim(25, 50, 25, scenario = "ar1", beta0 = c(-0.8, 0, 0, 0), beta1 = c(-0.8, 0, 0), eta = rep(0, 5))

challenges by domain top 3 or 4
- technical

outline

missingness - mar or not?

notification windows

non-response by age, gender

phone on/off? capture restart

tnecamp


- interpretation of beta_01 under xi = 0 or beta_11 = 0
- use of term conditional mean model

- in tables (i) and (ii)

- check in Mancl & DeRouen what is considered small (50?)