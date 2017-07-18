% [x, info] = quadprogpp (H)
% [x, info] = quadprogpp (H, g)
% [x, info] = quadprogpp (H, g, A, b)
% [x, info] = quadprogpp (H, g, A, b, lb, ub)
% [x, info] = quadprogpp (H, g, A, b, lb, ub, Ain, lbin, ubin)
% [x, info] = quadprogpp (..., options)
%
% Solve the quadratic program
%
%      min 0.5 x'*H*x + x'*g
%       x
%
% subject to
%
%      A*x = b
%      lb <= x <= ub
%      lbin <= A_in*x <= ubin
%
% Input
%
%   Variables g, A, b, lb, ub, Ain, lbin, ubin can be empty matrices.
%
%   Options may include field 'regularization_factor', value of which is added
%   to the diagonal elements of H to avoid failures on semidefinite problems,
%   by default regularization factor is equal to 1e-12.
%
% Output:
%
%   x -- the solution
%
%   info -- a structure containing the following fields:
%       status -- exit status
%           0 -- ok
%           1 -- the problem is infeasible
%       obj -- value of the objective function at the optimal point
%

function [x, info] = quadprogpp(varargin)
    H = [];
    g = [];
    A = [];
    b = [];
    lb = [];
    ub = [];
    Ain = [];
    lbin = [];
    ubin = [];
    regularization_factor = 1e-12;

    num_arg = nargin;
    if (isstruct(varargin{end}))
        if (isfield(varargin{end}, 'regularization_factor'))
            regularization_factor = varargin{end}.regularization_factor;
        end
        num_arg = num_arg - 1;
    end

    x = [];


% Check number of input parameters
    switch num_arg
        case 1
            H = varargin{1};
        case 2
            H = varargin{1};
            g = varargin{2};
        case 3
            error('Incorrect number of input variables.');
        case 4
            H = varargin{1};
            g = varargin{2};
            A = varargin{3};
            b = varargin{4};
        case 5
            error('Incorrect number of input variables.');
        case 6
            H = varargin{1};
            g = varargin{2};
            A = varargin{3};
            b = varargin{4};
            lb = varargin{5};
            ub = varargin{6};
        case 7
            error('Incorrect number of input variables.');
        case 8
            error('Incorrect number of input variables.');
        case 9
            H = varargin{1};
            g = varargin{2};
            A = varargin{3};
            b = varargin{4};
            lb = varargin{5};
            ub = varargin{6};
            Ain = varargin{7};
            lbin = varargin{8};
            ubin = varargin{9};
        otherwise
            error('Incorrect number of input variables.');
    end


% Check number of output parameters
    if (( 3 < nargout ) || ( 1 > nargout))
        error('Incorrect number of input variables.');
    end


% Check input correctness
    % objective
    if (isempty(H))
        error('Hessian is not initialized.');
    end

    if (size(H,1) ~= size(H,2))
        error('Hessian is not square.');
    end

    num_var = size(H, 1);


    if (isempty(g))
        g = zeros(num_var, 1);
    else
        if ( (size(g, 1) ~= num_var) || (size(g, 2) ~= 1) )
            error('Incorrect size of the constant vector.');
        end
    end


    % bounds
    if (~isempty(lb))
        if ((size(lb, 1) ~= num_var) || (size(lb, 2) ~= 1))
            error('Incorrect size of the vector of lower bounds.');
        end
    end

    if (~isempty(ub))
        if ((size(ub, 1) ~= num_var) || (size(ub, 2) ~= 1))
            error('Incorrect size of the vector of upper bounds.');
        end
    end


    % equality constraints
    if ( (isempty(A) && ~isempty(b)) || (~isempty(A) && isempty(b)))
        error('Partially initialized equality constraints.');
    end
    if (~isempty(b) && ~isempty(A))
        if ( (size(A, 1) ~= size(b, 1)) || (size(b, 2) ~= 1) || (size(A, 2) ~= num_var))
            error('Incorrect size of matrix or vector of equality constraints.');
        end
    end

    % inequality constraints
    if ( (isempty(Ain) && (~isempty(lbin) || ~isempty(ubin))) || (~isempty(Ain) && isempty(lbin) && (isempty(ubin))) )
        error('Partially initialized inequality constraints.');
    end
    if (~isempty(lbin) && ~isempty(Ain))
        if ( (size(Ain, 1) ~= size(lbin, 1)) || (size(lbin, 2) ~= 1) || (size(Ain, 2) ~= num_var) )
            error('Incorrect size of matrix or vector of inequality constraints.');
        end
    end
    if (~isempty(ubin) && ~isempty(Ain))
        if ( (size(Ain, 1) ~= size(ubin, 1)) || (size(ubin, 2) ~= 1) || (size(Ain, 2) ~= num_var) )
            error('Incorrect size of matrix or vector of inequality constraints.');
        end
    end


% Solve
    %minimize      1/2 x^ C x + c^x
    %subject to    a_j^x + b_j  =  0
    %              a_j^x + b_j  >= 0

    b = -b;


    Ain_tmp = [];
    bin_tmp = [];

    if (~isempty(ubin))
        Ain_tmp = [Ain_tmp; -Ain];
        bin_tmp = [bin_tmp; ubin];
    end

    if (~isempty(lbin))
        Ain_tmp = [Ain_tmp; Ain];
        bin_tmp = [bin_tmp; -lbin];
    end

    Ain_tmp = [Ain_tmp; eye(numel(lb)); -eye(numel(ub))];
    bin_tmp = [bin_tmp; -lb; ub];

    [x, info] = quadprogpp_interface(H + eye(size(H))*regularization_factor, g, A, b, Ain_tmp, bin_tmp);
end
