Introduction
============

Molecular dynamics (MD) simulations are exceptionally important in the research field of polymers, soft matters, biomolecules, and etc. 
In general, all-atom or coarse-grained force fields are not easily ported between MD packages. It is quit difficult to realize 
a new method in a well-developed MD package by users. Thereby, we provide a MD simulation package named PYGAMD 
(Python GPU-Accelerated Molecular Dynamics Software) to solve these problems.

PYGAMD is a platform where users could build up their customized force fields including potential forms and parameters.
This is achieved by that PYGAMD is programmed based on Numba, a Just-In-Time Python Compiler. The potential forms and methods could be conveyed 
from user interface to the underlying computation. More important, PYGAMD provide a high performance on GPU computation up to traditional packages 
programmed by CUDA and C. 

This package is the version 1 of PYGAMD which includes MD engine **pygamd**, molecular configuration generator **molgen**, and data tackler **dataTackle** etc. 
The **pygamd** is purely written by Python language based on Python3 Numba compiler. The plugins **molgen** and **dataTackle** that are
written by C++ and CUDA C, need to be compiled. 


