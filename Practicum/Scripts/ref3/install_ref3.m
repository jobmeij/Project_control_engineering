disp('Ref3 - Third degree reference generator')
disp('---------------------------------------')

addpath(cd)
savepath
disp(['  Added ' cd ' to matlab path.'])

mex ref3.c
rehash toolbox

disp(' ')
disp('Installation of Ref3 completed!')
disp('  See ref3_example.mdl for an example of how to use ref3.')
disp('  Mail comments to M.J.G.v.d.Molengraft@tue.nl or G.Witvoet@tue.nl.')
disp('  =================================================================')