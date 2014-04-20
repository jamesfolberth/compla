
% GNU Octave code
function [] = compare(savefile)

   if nargin == 0
      savefile = 'data.h5';
   end

   load(savefile);
   Lambda = transpose(Lambda);

   [e] = eig(A);
   abs(sort(e)(2)-sort(Lambda)(2)) % difference between 0.7... eval
   norm(sort(Lambda)-sort(e),inf) % same evals <=> small

end

