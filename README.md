Extinction and origination dynamics in Triassic-Jurassic
================
ALLuza, MG Bender, CS Dambros, L Kerber - Departamento de Ecologia e
Evolução, Universidade Federal de Santa Maria
2023-03-23

<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->
<!-- badges: end -->

# 1. Introduction

#### A bit about diversification and the fossil record

Extinctions are large depressions on Earth biodiversity history (Benton
1995). The diversity of all organisms on Earth, with an likely origin
during the Precambrian period (4560 to 570 Ma), has increased rapidly
and substantially after the Vendian and Early Cambrian (600 Ma) (Benton
1995). Since then, periods of diversification (i.e. the balance of
extinctions and speciations/originations) were marked by periods of high
rates of origination interrupted/interpassed by periods of biodiversity
crisis/depression, called ‘mass extinctions’ (Benton 1995). The fossil
record supports five mass extinction events, known as the ‘big five’:
end Permian, late Ordovician, late Devonian, Triassic-Jurassic boundary,
and Cretaceous-Paleogene boundary – that vanished from Earth a
considerable proportion of biodiversity, while speciation and radiation
to vacant niches occurred for remaining taxa from there on (Benton 1995,
May 1997, Reeves et al. 2020).

The Triassic-Jurassic extinction (Rhaetian-Hettangian intervals) played
an important role in the diversification rates of archosaurs, therapsid
reptiles, ammonites, among others (Allen et al. 2018), yet recent
research cast doubt on their existence (Bambach et al. 2004). The
statistical support for this extinction event (and its place as one of
the ‘big five’ extinctions) is limited due to the generally small sample
size at the Rhaetian-Hettangian intervals. Existing evidence suggests
that the lower number of taxa at this period can be a byproduct of high
extinction rates (e.g., Wignall & Atkinson 2020; Olsen et al. 2022)
–which seemed to exert an impressive effect on marine organisms (Benton
1995)– as well the low origination/speciation rates and nearly constant
extinction rates compared to previous time (Bamchach et al. 2004). Thus,
models that account for both uneven sampling size across periods and the
state of taxa occurrence at previous time could provide a more robust
inference on/shed light on the processes that caused the
Triassic-Jurassic extinction, whether low speciation/low extinction
rates or high extinction/low speciation rates.

#### The drivers of distribution

Understanding what drives shifts in temporal and spatial distribution of
taxa is one of the major goals in biology, disregarding whether the
studied taxa are extinct or extant. In the case of extinct taxa, the
abiotic environment and its variation over millennia years might affect
extinction and origination rates (e.g., Fritz et al. 2013), while
factors that hamper data collection and the taphonomic process (e.g.,
sedimentation, transport) might bias species detection in fossilized
samples (Liow 2013). The simultaneous accounting of ecological and
sampling processes yielding modeled data in hierarchical models has some
tradition in ecology (Williams et al. 2002, McKenzie et al. 2002, Kery
et al. 2013, Guillera-Arroita et al. 2016; reviewed by Kellner & Swihart
2014) but the use of these models to make inference in paleontology is
still in its infancy (Liow 2013, Lawing et al. 2021). Accounting for
these drivers of distribution might provide better maps, more robust
inferences and treatment of biological and observational uncertainty
producing paleontological data (Liow 2013 for a paleontological
perspective; Kery et al. 2013, Guillera-Arroita et al. 2016 show
ecological (neontological) perspectives).

#### Major aim

Here we used \~2 million fossil observations comprising 1139 genera and
1770 geologic formations, collected across the entire Triassic-Jurassic
periods, to estimate pairwise rates of extinction and origination across
17 time intervals of the Triassic and Jurassic periods (Olenekian (early
Triassic) to Tithonian (late Jurassic)). We emphasized i) temporal
variation on these rates, ii) the proportional change in number of genus
from the beggining of each interval to its end (Bambach et al., 2004),
and iii) the expected number of genera across these intervals. As the
limit between the Rhaetian (age) and Hettangian (age) marked the
Triassic-Jurassic extinction event and the split between these periods,
we expected to find i) higher rates of extinction and lower rates of
origination, ii) high change in the number of genus from the beggining
to the end of these periods, and iii) a lower expected number of taxa at
Rhaetian-Hettangian intervals when compared to other intervals in the
data.

# Methods

We used the Paleo Database (REF) as the main source of paleontological
data for the Triassic-Jurassic period. We used information at the genus
level, as genera give a closer view of species-level behavior through
time, despite comprising less complete data than one at a higher
taxonomic level because of the vagaries of preservation, collection, and
study (Benton 1995). Nonetheless, we accounted for such sampling biases
in our models.

## State-space models

State-space models are hierarchical models that enable the estimation of
biological parameters describing dynamic properties of ecological
systems (Kéry & Royle 2007, Ecology). They are called ‘first-order
Markovian models’ because they use information from previous time and/or
neighbor sites to project biological parameters in the next time. Thus,
they are useful to estimating parameters involved with inference in
paleontology, such as extinction (*epslon*), origination (*gamma*), and
turnover rates (Bambach et al. 2004, Benton 1995). In state-space
models, the true state of a system (e.g. the realized occurrence of the
fossil genera *g* at the time interval *t*) is denoted by the latent,
stochastic variable *z<sub>\[t,g\]*. The *z</sub>\[t,g\]* is the
realization of a Bernoulli process (‘dbern’) dependent on the variable
*muZ\~\[t,g\]*, which is the expected probability of occupancy of
interval *t* by the genus *g*. In its turn, the *muZ\~\[t,g\]* depends
on both the realized occurrence of genus *g* in the previous time *t-1*,
ie. *z\~\[t-1,g\]* (whose first time gives the initial conditions for
the dynamics develop) and the estimated probability of extinction
(*epslon\~t* = 1 - *phi<sub>t*) and origination (*gamma</sub>t*) in the
next interval if the previous interval *t-1* is either occupied or not
by the genus *g*. The model design looks like:


        model {
        
        #############################################################
        #                                                           #
        #                  Biological process                       #
        #                                                           #
        #############################################################
        
        
        # Priors
            
        for (t in 2:nint) {
        
           for (g in 1:ngenus) {
        
            
              ## colonization (origination)
              gamma [t,g] ~ dunif (0,1)
              ## extinction
              phi [t,g]~ dunif (0,1)
            
            }
           
           }
      
        
        ## set initial conditions
        ## priors for occupancy in time 1
         for (g in 1:ngenus) {
         
          psi1 [g] ~ dunif (0,1)
         
         }
        
        ############ Model ########################
        for (g in 1:ngenus) {
        
            z[1,g]~dbern(psi1[g]) # occupancy status initialization
        
                for (t in 2:nint){
                
                  # model likelihood
                  ### modeling dynamics conditional on previous time realized occurrence z
                  muZ[t,g] <- z[t-1,g] * (1-phi[t,g]) + ### if occupied, p of not getting extinct in the next time
                              (1-z[t-1,g]) * gamma[t,g] ###  if not occupied, p of getting colonized in the next time
                            
                z[t,g] ~ dbern(muZ[t,g])
        
            }#t
        }#i
        
        #############################################################
        #                                                           #
        #         Observation process across formations             #
        #                                                           #
        #############################################################
        
        # priors
        ## detection
        #for (t in 1:nint) {
        
          for (g in 1:ngenus) {
        
                p[g] ~ dunif (0,1)

          }
        #}

        ### model
        for (k in 1:nobs){
            
            y [k] ~ dbern (muY[form[k],int[k],genus[k]])
            muY [form[k],int[k],genus[k]] <- z [int[k],genus[k]] * p[genus[k]]
            
        }
        
        
        # -----------------------------------------------
        
        ## derived parameters
        # number of genus per interval
        for (t in 1:nint) {
            Ngen[t]<-sum(z[t,])
        }
        
        # average extinction and origination
        for (g in 1:ngenus) {
          avphi[g] <- mean(phi[2:nint,g])
          avgamma[g]<- mean(gamma[2:nint,g])
        }
          
        # turnover (proportional gain or loss)
        for (t in 2:nint) {  
          
            propcH [t] <-(sum (z[t-1,]) - sum(z[t,]))/sum(z[t-1,]) 
          
        }
        
        # equilibrium occupancy (which genus decline or increase over time)
        for (g in 1:ngenus) {
        
            psi.eq[g] <- mean(gamma[2:nint,g])/(mean(gamma[2:nint,g])+mean(1-phi[2:nint,g])) # Equilibrium occupancy
        
        }
        
        
        }## end of the model
        
        

Besides allowing the estimation of biological parameters depicting taxon
dynamics over time, state-space models enable the estimation of
parameters depicting the probability of detecting a genus *p\~\[g\]* in
any geological formation and interval conditional on the genus realized
occurrence *z\~\[t,g\]*. The new detection/detection dataset (or the
data that should be observed given the realization of the biological and
observational processes) of genus *g* in interval *t* and geological
formation *k*, denoted by the latent variable *y<sub>\[g,t,k\]*, is
another parameter estimated by the state-space model. The
*y</sub>\[g,t,k\]* is the realization of a Bernoulli process which
depends on *muY<sub>\[g,t,k\]*, i.e. the probability of detecting a
genus conditional on its realized occurrence (i.e., *z</sub>\[t,g\] x
p<sub>\[g\]*). In this model, *p</sub>\[g\]* only varies across *g* to
*G* genera. We estimated genus detection probability across geological
formations as each one of them represent a fieldwork effort to obtain
fossil data, and might represent differential physical processes causing
sedimentation/erosion, transport, deformation, and lost of
paleontological/taphonomic information/affect fosssil preservation. The
input data consist of detections and non-detections (coded as 1 and 0,
respectively) of each genus in each geological formation and interval. A
genus was ‘detected’ if it appeared at least once in a geological
formation and interval, and it was ‘not-detected’ if it not appeared in
formations and intervals with at least one detected genus. Missing
observations (coded as ‘NA’) comprise intervals and formations with no
detected genus, and were removed from the state-space model.
Nonetheless, the model estimated data for these missing observations
using the existing data and the prior information (see below “Modeling
approach”).

## Derived parameters

Using the state-space model structure we were able to estimate
parameters often used to make inference in paleontology, namely
‘extinction probability’, ‘origination probability’, ‘proportional
change’, and ‘equilibrium occupancy’. The first ones depict the
probability of genus extinction and origination from t to t+1 (i.e.,
between each pair of adjacent time intervals). There were 16 probability
estimates from Olenekian (early Triassic) to Tithonian (late Jurassic),
as there is no transition from Olenekian to previous time
(Changhsingiano, late Permian). We choose not to include Permian data in
our state-space model because it would include a large taxonomic
heterogeneity in our data, as the Permian-Triassic boundary (\~252 Ma)
was characterized by the most massive extinction on Earth’s history
(Benton 1995).

The proportional change (PCh) (Bambach et al. 2004) depicts the
difference in the number of taxa in the end relative to the beginning of
a given period, as follows:

*PCh\[t\]* = sum (z\[t-1,\]) - sum(z\[t,\]))/sum(z\[t-1,\])

Finally, the equilibrium occupancy showing the trend of a genus *g* to
either decline or increase over time was estimated as follows:

*psi.eq\[g\]* =
mean(gamma\[t:T,g\])/(mean(gamma\[t:T,g\])+mean(1-phi\[t:T,g\]))

## Modeling approach

We estimated the state-space model parameters using Bayesian inference.
Bayesian inference is a statistical method that applies the Bayes
Theorem to update prior information about parameter values with observed
data to then estimate the posterior probability of parameters. The prior
information, generally represented by statistical distributions, are
updated by their integration with data and likelihood estimation across
independent Monte Carlo Markov Chains (MCMC). Variation in parameter
estimates across posterior distribution draws represents an appropriate
measure of parameter uncertainty (Kruschke & Liddell, 2018). In this
sense, the Credible Intervals (CI), built using the quantiles of the
posterior distribution draws of each parameter, depict the interval
where most posterior distribution draws are and, therefore, delimit the
area in which we have large certainty of finding the true parameter
value (Kruschke & Liddell, 2018) \[COPY OF MY JBI PAPER\].

The priors for all model parameters were assumed to come from an Uniform
distribution ranging from 0 to 1 (*phi<sub>\[t,g\]*,
*gamma</sub>\[t,g\]*, *psi<sub>\[t=1,g\]*, *p</sub>\[g\]* \~ U(0,1)).

Our state-space model was run using three independent MCMC (using the
Gibbs Sampling algorithm implemented in JAGS; Plummer, 2003), with
15,000 iterations and a warm up period of 7,000 iterations each chain.
We retained 1000 post-warm up draws of each chain, resulting in 3000
posterior distribution draws of each model parameter which were used to
make statistical inference and test our hypotheses.

# Results

The analyzed dataset comprises 1139 genera observed at 17 time intervals
and 1769 geological formations. Thus, if all intervals and geological
formations were sampled we would have 30.073 possible observations per
taxon/genus with this dataset. However, the dataset is quite sparse,
with many missing observations, and had only 7.1 % (n=2.149
observations) out of the possible observations per taxon.
