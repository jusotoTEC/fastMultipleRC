# MATLAB Code for Numerical Experiments in the Paper "*Fast Multiple Rank-Constrained Matrix Approximation*"

## Authors

* Pablo Soto-Quiros (https://www.tec.ac.cr/juan-pablo-soto-quiros) - Email: jusoto@tec.ac.cr
* Jeffry Chavarría-Molina (https://www.tec.ac.cr/jeffrey-chavarria-molina) - Email: jchavarria@tec.ac.cr
* Juan José Fallas-Monge (https://www.tec.ac.cr/juan-jose-fallas-monge) - Email: jfallas@tec.ac.cr
* Anatoli Torokhti (https://people.unisa.edu.au/anatoli.torokhti) - Email: anatoli.torokhti@unisa.edu.au

Pablo Soto-Quiros; Jeffry Chavarría-Molina and Juan José Fallas-Monge are Associate Professors from the *Instituto Tecnológico de Costa Rica* (https://www.tec.ac.cr/) in Cartago, Costa Rica

Anatoli Torokhti is an Associate Professor from the *University of South Australia* (https://www.unisa.edu.au/) in Mawson Lakes, SA, Australia

## Description

* This repository contains the MATLAB code for numerical experiments presented in the paper "*Fast Multiple Rank-Constrained Matrix Approximation*". 
* This scientific paper has been published in the "SeMA Journal".
* The link of the paper is [https://doi.org/10.1007/s40324-023-00340-6](https://doi.org/10.1007/s40324-023-00340-6).
* Our work addresses  methods for a fast  multiterm  matrix approximation subject to multiple rank constraints. The problem arises in applications associated with data processing systems. For large matrices,  finding  acceptable matrix approximations  may require quite a long time. In practice, this issue may fail   associated with computation due to a conflict with available time and computer memory.
We provide  techniques to accelerate the associated computation and avoid the above  bottleneck. The proposed approach combines a fast pseudoinverse matrix computation, based on the use of a vector tensor product, with a fast low-rank matrix approximation, based on a new extension of a method of bilateral random projections. The provided theoretical and numerical studies  demonstrate the faster performance  of the proposed method  compared to  methods based on the SVD computation. It is achieved, in particular, at the cost of "a little bit" worse associated numerical error which, in many practical cases, might be acceptable.


<p align="center"><img width="1200" src="https://github.com/jusotoTEC/fastMultipleRC/blob/main/img/img.jpg"></p>
