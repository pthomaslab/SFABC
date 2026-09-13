function feasibility = feas_test_tele_steady(us, tol, se, p, Ps, ...
                                        k_off, k_on, k_syn, k_deg, ...
                                        Q_off, Q_on, Q_syn, Q_deg, I)
%FEAS_TEST_TELE_STEADY Test steady-state CME feasibility for the telegraph model.
%
%   feasibility = feas_test_tele_steady(...) checks whether the truncated
%   state probability vector p can satisfy the CME and the moment bounds
%   around candidate probabilities Ps for the supplied rates and tolerance.
%
%   Output:
%     feasibility - Solver status code (0 = feasible, 1 = infeasible).
%                   For other errors, see
%                   https://yalmip.github.io/command/yalmiperror/.

%% Configure the feasibility problem
% d is the highest moment order; B is the number of unique observed states.
d = length(se) - 1;
B = length(us);

% Configure MOSEK with silent output and basis identification disabled.
ops = sdpsettings('solver','mosek','verbose',0,...
    'mosek.MSK_IPAR_INTPNT_BASIS','MSK_BI_NEVER');

%% Build constraints
% Require nonnegative probability mass with total mass at most 1.
Fr = [p>=0, sum(p)<=1];

% Enforce the truncated steady-state CME: Q*p = 0.
Fr = Fr + [k_off * (Q_off * p) + k_on * (Q_on * p) + ...
    k_syn * (Q_syn * p) + k_deg * (Q_deg * p) == 0];

% Aggregate joint gene-state/mRNA probabilities to marginal mRNA probabilities.
p_agg = I*p;

% Bound standardized moment differences for orders 0,...,d.
for i = 0:d
    Fr = Fr + [-tol <= (sum(p_agg(us+1)'.*us.^i) - ...
        sum(Ps(1:B).*us.^i))/se(i+1) <= tol];
end

%% Solve the feasibility problem
% Use a constant objective subject to constraints Fr.
sol = optimize(Fr,1,ops);
feasibility = sol.problem;
end
