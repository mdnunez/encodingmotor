function printiter(iter)
%PRINTITER- Overwrites previous iter-1 output and replaces with iter
%
%Usage: printiter(n)
%
%Input: n - Integer value
%
%-This function is useful within algorithms if a progress report is 
%required
%-Note that this will not work if any other output is required during your
%loop
%
%Usage example (within a script):
% fprintf('Testing ''printiter'' function at iteration -101');
% for n=(-100):100
%   pause(.1);
%   printiter(n);
% end
% fprintf('\n');

% Copyright (C) 2015 Michael D. Nunez, <mdnunez1@uci.edu>
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


%% Record of Revisions
%   Date           Programmers               Description of change
%   ====        =================            =====================
%  5/20/15        Michael Nunez                  Original Code

%% Code
if (iter == round(iter))
    fprintf([repmat('\b',1,numel(num2str(iter-1))),'%i'],iter);
else
    fprintf('?');
end
