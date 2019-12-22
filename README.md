# fracdiff
fractional differencing in C with frequency domain discussion:  the properites of the fractional differencing and integration filters in the frequency domain

Plain C code, set up to build on Mac, but should be buildable by any C compiler.

See detailed notes in code and PDF slide show:  https://github.com/diffent/fracdiff/blob/master/freqrespfracdiff.pdf

For a slightly different take on fractional time series modeling methods, see the Fractional Empirical Motion
feature of our MCarloRisk3D iPad app, where the Hurst exponent can be specified for projecting
forward price / probability envelopes.  FEM is a generalization of FBM (Fractional Brownian Motion)
which does not assume a Gaussian for the underlying distribution, but rather allows
an empirical distribution to be used:

https://apps.apple.com/us/app/mcarlorisk3d/id641208540

In particular, see slide set 9 in this training guide for examples, and how to access the Hurst parameter:
https://diffent.com/mcrtrain/mcrtrain.html
