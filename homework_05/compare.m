
% GNU Octave code
function [] = compare(savefile)

   if nargin == 0
      savefile = 'data.h5';
   end

   load(savefile);
   Lambda = transpose(Lambda);

   [e] = eig(A);
   norm(sort(Lambda)-sort(e),inf) % same evals <=> small

end

