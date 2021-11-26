<head>
    <script src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script>
    <script type="text/x-mathjax-config">
        MathJax.Hub.Config({
            tex2jax: {
            skipTags: ['script', 'noscript', 'style', 'textarea', 'pre'],
            inlineMath: [['$','$']]
            }
        });
    </script>
</head>

# HI_ProfileLine
Plot a HI profile line from TNG-100 simulation

# TODO
 - [ ] Luminosity distance fix: single fixed value -> value for each particle
 - [ ] LOSVector fix: similar to lumnosity distance fix
 - [ ] Calculate the rotational matrix of subhalo

# Logs

## 11/26/2021
Still dealing with previous problem. Some critical bug fixed, but result is still far from expectation. $frac{M_{simu}}{M_{real}}}$ seems like a sine wave and exceed 1, which is against prediction.

Trying to get rotate_star working. I'm not sure about how it works.

Remove class prototype 'coord' from subhalo.py (this class is probably not be used).

## 11/23/2021
Luminosity distance & LOSVector now is different for every gas particle.
However, the outcome is not ideal or, I would say totally opposite to anticipation.
Debugging...

## 11/22/2021
Profile line -> HI simulate mass function complete

## 11/19/2021
Generate 26 plot instead of 33 plot (better mapping)

## 11/18/2021
One-key generate 32 + 1 plot for a single subhalo.

## 11/17/2021
Main implentation of class subhalo completed.

## 11/16/2021
Reconstructing.
Looks like currently there is no need for coordinate class.

Create class structure for subhalo reading/ploting.
Utility for coordinate system change
## 11/3/2021
Reason for last error: no gas particle exist in some subhalos. Add check if gas exists.
## 11/1/2021
Known error for snapShot 99, subhalos with higher ordinal returns KeyError: 'GFM_Metals'. Need to check the shape of GFM_Metals for 43
Unexpected shape for subhalo: 
1. Probabaly too large?
2. Symmetry problem
7. Seems like not much HI present here 
42. Symmetry problem
43. Seems like not much HI present here 
