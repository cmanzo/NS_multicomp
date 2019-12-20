% % Copyright (c) 2019-2020 Carlo Manzo & Tina Košuta 

% If you use this code, please cite: 
% (http://dx.doi.org/10.1039/C9CP05616E)
% (https://arxiv.org/abs/1909.13133)

% Permission is hereby granted, free of charge, to any person 
% obtaining a copy of this software and associated documentation files 
% (the "Software"), to deal in the Software without restriction, 
% including without limitation the rights to use, copy, modify, merge, 
% publish, distribute, sublicense, and/or sell copies of the Software, 
% and to permit persons to whom the Software is furnished to do so, 
% subject to the following conditions:

% The above copyright notice and this permission notice shall be 
% included in all copies or substantial portions of the Software.

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
% IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
% CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
% TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
% SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

function y = logdirichlet(x, delta)
if length(delta)==1 % simmetric Dirichlet
    delta=delta*ones(1,size(x,2));
end    
delta=repmat(delta,[size(x,1),1]);

y=sum((delta-1).*log(x),2)  -sum(gammaln (delta) ,2)+ gammaln(sum(delta,2));