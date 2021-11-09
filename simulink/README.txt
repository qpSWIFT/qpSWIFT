-------------------------------------------------------------------------------
There are two ways of interacting with qpSWIFT simulink interface as
1) For general quadratic programs with no equality constraints
2) For general quadratic programs with equality constraints

Important Considerations
-> Change the value of the output NV to the desired number of output variables
-> Use qpdemo.slx and qpdemo_e.slx for model problems (1) and (2)
-> 




-------------------------------------------------------------------------------
Compilation -> Type Swift_make('qpSWIFT') in your matlab command window
            -> Add the mex file to your working directory to use qpSWIFT
-------------------------------------------------------------------------------
Usage       -> Instructions on how to use the mex-file are given in qpSWIFT.m
                                    or
            -> Type help qpSWIFT in your matlab command window
-------------------------------------------------------------------------------
Demo        -> Demo QP is given in demoqp.slx and demoqp_e.slxx
-------------------------------------------------------------------------------


-------------------------------------------------------------------------------
Note: Make sure you have compatible C compiler available for your matlab version
-------------------------------------------------------------------------------