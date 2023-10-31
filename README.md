# Takagi
Matlab solver for Riemannian optimization on the Stiefel manifold

## Problems
### This solver is to solve the following optimization problem
     min -1/2 tr(Jp \widetilde{U}' Jn \widetilde{A} \widetilde{U} \widetilde{N}) ,
     s.t.   \widetilde{U}' \widetilde{U} = I2p.
  
where \widetilde{U} is a 2n-by-2p matrix, Jn = [0 In; In 0], and I2p is the 2p-by-2p identity matrix.


--------------------------------------------------------------------------------
Contents:
--------------------------------------------------------------------------------


| name  | meaning |
| ---------- | -----------|
| Our\main_singular.m    | Script file that can run the experiments of truncated Takagi factorization   |
| Our\fun_singular.m    | Objective function and its Euclidean gradient of truncated Takagi factorization  |
| Our\wen_Stiefel_singular.m   | Algorithm 4.1 (the steepest descent method) of truncated Takagi factorization   |
| Our\zhu_StiQR_singular.m   | Algorithm 4.2 (the Riemannian nonmonotone conjugate gradient method) of truncated Takagi factorization   |
| Our\Newton.m   | Algorithm 4.3 (the Newton mathod) of truncated Takagi factorization  |
| Our\xishu.m   | The coefficient matrix M for solving Newton's equation   |
| Our\zhu_StiQR_singular2.m   | Sub-code for Algorithm 4.4 (Same as Algorithm 2 except for eps in stopping criterion) of truncated Takagi factorization   |
|Sato's\main_singular.m  |   Script file that can run the experiments |
|Sato's\fun_singular.m | Objective function and its Euclidean gradient of truncated svd|
|Sato's\Alg1.m    | Sato's steepest descent method of truncated svd|
|Sato's\Alg2.m  |Sato's conjugate gradient method of truncated svd   |
|Sato's\Alg3.m | Sato's Newton method  of truncated svd |
|Sato's\CR.m | Sub-code for Alg3（solve the system of linear equations）|
| README.txt   | this file   |

## Notes:
Table 5.1 and Figure 5.1 show the complex symmetric matrix A of Example1 by running "main_singular.m" in the "Our" folder and "Sato's" folder, respectively.

Table 5.2 and Figure 5.2 are calculated by running "main_singular.m" in the "Our" folder several times, and the reader can adjust the dimensionality by himself/herself.

## References
+ Kong L and Chen X. Riemannian optimization methods for the  truncated Takagi factorization. Numerical Algorithms, 2023
+ Sato H and Iwai T. A Riemannian optimization approach to the matrix singular value decomposition, SIAM J Optim, 23(2013): 188--212.
+ Stao H and Iwai T. A complex singular value decomposition algorithm based on the Riemannian Newton method, In: 52nd IEEE Conference on Decision and Control, (2013): 2972--2978.
+ Sato H. Riemannian conjugate gradient method for complex singular value decomposition problem, In: 53rd IEEE Conference on Decision and Control, (2014): 5849--5854.

## Authors
Lingchang Kong, Xiaoshan Chen




## Copyright
Copyright (C) 2023, Lingchang Kong, Xiaoshan Chen

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see [http://www.gnu.org/licenses/](http://www.gnu.org/licenses/)


