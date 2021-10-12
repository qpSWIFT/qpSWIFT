function Swift_make(file)

curpath = cd;
curpath = curpath(1:end-6); % remove matlab from the curpath
if(isempty(file))
    disp('Enter the name of the s_function to be used');
else
% d = '-largeArrayDims -DDLONG -DLDL_LONG';
% d = '-DDLONG -DLDL_LONG';
%% ldl solver
    fprintf('Compiling LDL...\n');
    cmd = ['mex -g -c -largeArrayDims ','-I',curpath,'/include ','-I',curpath,'ldl/include ',curpath,'ldl/src/ldl.c'];
    eval(cmd);
    fprintf('LDL Compilation Done\n');
    
%% AMD Solver
    fprintf('Compiling AMD...\n');
    cmd = ['mex -c -g -largeArrayDims ','-I',curpath,'amd/include ',curpath,'amd/src/amd_order.c ',curpath,'amd/src/amd_dump.c ',curpath,'amd/src/amd_postorder.c ',...
            curpath,'amd/src/amd_post_tree.c ', curpath,'amd/src/amd_aat.c ',curpath,'amd/src/amd_2.c ',curpath,'amd/src/amd_1.c ',curpath,'amd/src/amd_defaults.c ',curpath,'amd/src/amd_control.c ', ...
             curpath,'amd/src/amd_info.c ', curpath,'amd/src/amd_valid.c ', curpath,'amd/src/amd_global.c ', curpath,'amd/src/amd_preprocess.c ' ];
    eval(cmd);
    fprintf('AMD Compilation Done\n');
%% QP (core IPM solver)
    fprintf('Compiling QP...\n');
    cmd = ['mex -c -g -largeArrayDims ','-I',curpath,'/include ','-I',curpath,'ldl/include ' ,'-I',curpath,'amd/include ',curpath,'src/Auxilary.c ',curpath,'src/Prime.c ',curpath,'src/timer.c '];
    eval(cmd);
    fprintf('QP Compilation Done !!!\n');
 
%% QP_mex
    fprintf('Compiling QP_sfunc...\n');
    cmd = ['mex -c -g -largeArrayDims ','-I',curpath,'/include ','-I',curpath,'ldl/include ', '-I',curpath,'amd/include ',file,'.c'];
    eval(cmd);
    fprintf('QP_sfunc Compilation Done\n');
    fprintf('Linking Object Files...\n');
    if (ispc) cmd = ['mex -g -largeArrayDims ','Prime.obj Auxilary.obj timer.obj ldl.obj amd_info.obj amd_1.obj amd_2.obj amd_control.obj amd_defaults.obj amd_aat.obj amd_global.obj amd_post_tree.obj amd_postorder.obj amd_preprocess.obj amd_dump.obj amd_order.obj amd_valid.obj ',file '.obj',' -output ', file]; end
    if (isunix) cmd = ['mex -g -largeArrayDims ','Prime.o Auxilary.o timer.o ldl.o amd_info.o amd_1.o amd_2.o amd_control.o amd_defaults.o amd_aat.o amd_global.o amd_post_tree.o amd_postorder.o amd_preprocess.o amd_dump.o amd_order.o amd_valid.o ',file '.o',' -output ', file]; end
    eval(cmd);
    fprintf('Successfully Compiled \n');

%% clean
    if( ispc ), delete('*.obj'); end
    if( isunix), delete('*.o');end
    fprintf('Object Files cleaned \n');
    clear;
end