function p = polyfitweighted2(x,y,z,n,w)
%polyfitweighted2: Fit polynomial to data.
%   P = polyfitweighted2(X,Y,Z,N,W) finds the coefficients of a polynomial P(X,Y) of
%   degree N that fits the data Z best in a least-squares sense. P is a
%   row vector of length (N+1)*(N+2)/2 containing the polynomial coefficients in
%   ascending powers, p00 p10 p01 p20 p11 p02 p30 p21 p12 p03...
%
%   X,Y must be vectors
%   Z,W must be 2D arrays of size [length(X) length(Y)]
%
%   based on polyfit.m - see doc polyfit for more details
%
%   Class support for inputs X,Y,Z,W:
%      float: double, single
%
%   by SSR (2006) based on polyfit.m by The MathWorks, Inc.
%   
%
% The regression problem is formulated in matrix format as:
%
%    z = V*p    or
%
%                   2      2  3  2      2   3
%    z = [1  x  y  x  xy  y  x  x y  x y   y ]    [p00
%                                                  p10
%                                                  p01
%                                                  p20
%                                                  p11
%                                                  p02
%                                                  p30
%                                                  p21
%                                                  p12
%                                                  p03]
%
% where the vector p contains the coefficients to be found.  For 4th 
% order quadric surface, matrix V would be:
%
% V = [w wx wy wx2 wxy wy2 wx3 wx2y wxy2 wy3 wx4 wx3y wx2y2 wxy3 wy4];

x = x(:);
y = y(:);

lx=length(x);
ly=length(y);

if ~isequal(size(z),size(w),[ly lx])
    error('MATLAB:polyfitweighted2:XYSizeMismatch',...
          [' X,Y *must* be vectors' ...
          '  Z,W *must* be 2D arrays of size [length(X) length(Y)]'])
end

y=y*ones(1,lx);
x=ones(ly,1)*x';
x = x(:);
y = y(:);
z = z(:);
w = w(:);

pts=length(z);

% Construct weighted Vandermonde matrix.
V=zeros(pts,(n+1)*(n+2)/2);
V(:,1) = w;
%V(:,1) = ones(pts,1);
ordercolumn=1;
for order = 1:n
    for ordercolumn=ordercolumn+(1:order)
        V(:,ordercolumn) = x.*V(:,ordercolumn-order);
    end
    ordercolumn=ordercolumn+1;
    V(:,ordercolumn) = y.*V(:,ordercolumn-order-1);
end

% Solve least squares problem.
[Q,R] = qr(V,0);
ws = warning('off','all'); 
p = R\(Q'*(w.*z));    % Same as p = V\(w.*z);
warning(ws);
% if size(R,2) > size(R,1)
%    warning('MATLAB:polyfit:PolyNotUnique', ...
%        'Polynomial is not unique; degree >= number of data points.')
% elseif condest(R) > 1.0e10
%         warning('MATLAB:polyfit:RepeatedPointsOrRescale', ...
%             ['Polynomial is badly conditioned. Remove repeated data points\n' ...
%             '         or try centering and scaling as described in HELP POLYFIT.'])
% end
%r = z - (V*p)./w;
p = p.';          % Polynomial coefficients are row vectors by convention.
