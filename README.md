# qpSWIFT
Light-weight sparse Quadratic Programming Solver


## Introduction
qpSWIFT is light-weight sparse Quadratic Programming solver targetted for embedded and robotic applications. It employs Primal-Dual Interioir Point method with Mehrotra Predictor corrector step and Nesterov Todd scaling. For solving the linear system of equations, sparse LDL' factorization is used along with approximate minimum degree heuristic to minimize fill-in of the factorizations

## Wiki
For more information, please check the repo wiki (in progress).

## Problem Structure
qpSWIFT is designed to solve Quadratic Programs of the following form \
`min. 0.5*x'Px + c'x`\
`s.t  Ax = b `\
`     Gx <= h `

## Features
 - Written in ANSI-C
 - Fully functional Quadratic Programming solver for embedded applications
 - Code Generation for target platform
 - Tested on multiple target architectures
    + x86
    + x86_64
    + ARM
  - Support for multiple interfaces
    + [C/C++](https://github.com/qpSWIFT/qpSWIFT/tree/main/src)
    + [Matlab](https://github.com/qpSWIFT/qpSWIFT/tree/main/matlab)
    + [Simulink](https://github.com/qpSWIFT/qpSWIFT/tree/main/simulink)


<!---
# Case Studies
  - BAMBY
  - Ghost Robotics Vision60
-->

## Future Updates
  - Quadratic Program with only equality constraints

## Note
**The project is still in active development. Feedback is highly appreciated. For any queries and suggestions please write to agp19@vt.edu, yanran@mit.edu or haewonpark@kaist.ac.kr**

## Citing qpSWIFT
If you like qpSWIFT and are using it in your work, please cite the following paper\
  `@article{pandala2019qpswift,`\
  `title     = {qpSWIFT: A Real-Time Sparse Quadratic Program Solver for Robotic Applications},`\
  `author    = {Pandala, Abhishek Goud and Ding, Yanran and Park, Hae-Won},`\
  `journal   = {IEEE Robotics and Automation Letters},`\
  `volume    = {4},`\
  `number    = {4},`\
  `pages     = {3355--3362},`\
  `year      = {2019},`\
  `publisher = {IEEE}`\
  `}`
