function [B1,bound2] = preallocated(Grid,fspace,bN,fgN,QP)
% This function: computes pre-allocated variables of the NKM such as
% natural output (Yn) and next period shocks (A1 and B1)
%
% Inputs:
%     - Grid:     matrix of grid points
%     - fspace:   definition of interpolation scheme
%     - aN:       number of points in first state variable
%     - bN:       number of points in second state variable
%     - QP:       vector of quadrature nodes
% Output:
%     - Yn:       vector of natural output
%     - A1:       TFP in period t+1
%     - B1:       discount factor in period t+1
%     - bound1:   position of points in which TFP is binding
%     - bound2:   position of points in which discount factor is binding

% Read parameters
[~,~,~,~,BETA_LR,~,~,~,RHO_B] = parameters;

% b) Next period discount factor
B1 = kron((1-RHO_B)*log(BETA_LR)+RHO_B*Grid(:,1),ones(length(QP),1))+kron(ones(bN*fgN,1),QP);
bound2 = (B1>fspace.b(1)) + (B1<fspace.a(1));
B1 = min(max(B1,fspace.a(1)),fspace.b(1));