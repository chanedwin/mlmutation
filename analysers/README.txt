----------------------------------------------------------------------
This pipeline was written by Edwin Chan for his Final Year Project, 2017.

Problem Statement :

A lot of mutations called by variant callers are false positives.
We want a way to try to reduce the number of false positives using machine learning techniques.
We also want a way to rank the mutations in a manner that is understandable.

Solution :

Use a deep learning method to integrate the information from different variant callers.
Use a bayesian network to rank the genes.

Folder Structure :



It contains a set of tools which can be used to do two things -

1. use a machine learning method to 