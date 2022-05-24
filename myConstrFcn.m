function [cneq, ceq, gradc, gradceq] = myConstrFcn(x)
 global nelx nely vol volfrac ang angle  penal rmin
  nele=nelx*nely;
    % Non-linear Constraints
    cneq  = sum(vol) - volfrac*nele;
    gradc = ones(nele,1);
    % Linear Constraints
    ceq     = [];
    gradceq = [];
end % mycon
