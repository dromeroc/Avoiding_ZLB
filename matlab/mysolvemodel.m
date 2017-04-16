function [G1,impact,eu] = mysolvemodel(funcname,stst,varargin)
% This function gives the linear rational expectation solution for any
% model using Sims' package.
%
% Input
%
%   -funcname:  function name of the residual equations of the model
%   -stst:      vector of variables where the approximation must be taken
%               (typically the steady state of the model)
%   -varargin:  variable arguments (parameters vector, number of variables,
%               number of shocks and number of expectational errors)
%
% Output
%
%   -G1:        policy function with respect to variables
%   -impact:    policy function with respect to shocks
%   -eu:        2x1 Boolean vector for existence and uniqueness of solution

% Parameters
params  = varargin{1};
nvar    = varargin{2};
nshocks = varargin{3};
nEerr   = varargin{4};

% Check model at steady state
check = norm(funcname(stst,params),inf);
if check>=1e-5;error('Wrong steady state/equations');end

% Compute numerical derivatives
jac = fdjac(funcname,stst,params);
g0  = -jac(:,1:nvar);
g1  = jac(:,nvar+1:2*nvar);
c   = funcname(stst,params);
Psi = jac(:,2*nvar+1:2*nvar+nshocks);
Pi  = jac(:,2*nvar+nshocks+1:2*nvar+nshocks+nEerr);

% Check stationarity
div = 1;
eu = checkeu(g0,g1,Psi,Pi,div);
if sum(eu)~=2;warning('Problem with existence or solution');end

% Find policy functions
[G1,~,impact] = gensys(g0,g1,c,Psi,Pi,div);