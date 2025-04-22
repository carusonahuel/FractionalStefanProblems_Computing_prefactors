%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%         Author: Nahuel Caruso  (Cifasis-Conicet / Fceia-UNR)
%         mail: caruso.nahuel@gmail.com
%         Year: 2024
%
%         Especial Function: 
%
%         Creation:         2018-07-01
%         Modifications:    2019-01-07 
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function W = build_up_betainc_series_function(beta,varargin)
   %
   kmax = varargin{1};
   %if (nargin==3)
   %    kmax = varargin{2};
   %end
   k = [0:kmax]; 
   %
   W = @(x) sum(bsxfun(@power,repmat(x,[1,kmax+1]), k)*beta,2);
   %
end 