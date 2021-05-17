function [sol, num_iter, bwd_err] = jacobi(A, b, max_iter, e, init_opt)
% input   A        matrix
%         b        rhs vector
%         max_iter maximum number of iterations
%         e        error tolerance
%         init_opt initial vector option
%
% output  sol        solution vector
%         bwd_err    error in infinity norm
%         num_iter   number of iterations executed

	% Solve Ax = b
    % A = D + L + U
	% D = The diagonal entries of A
	D = diag(A);
    if ~all(D) 
        error 'at least one diagonal entry is zero';
    end
    % Inverse of the diagonal entries of A
	D_inv = D.^-1;
    
	% L + U = A - D = Matrix of off-diagonal entires of A
	L_plus_U = A - diag(D);

	% Compute D.^-1*b 
	D_inv_b = D_inv.*b;
    
    % Initialize x0 as D.^-1*b by default
    x_0 = D_inv_b;
    switch init_opt
        case 0
            x_0 = zeros(length(b),1);
        case 1
            x_0 = ones(length(b),1);
        case 2
            x_0 = D_inv_b;
    end
    
	% Iterate x = D^{-1} (b - (L+U)*x)
    iter = 0;
    diff = 1;
    x_old = x_0;
	
    while (diff > e && iter < max_iter)
        iter = iter + 1;
        % Update the approximate solution x
        x_new = D_inv_b - D_inv.*(L_plus_U*x_old);
        diff = norm(x_new - x_old, inf);
        
        if (diff >= e)
			x_old = x_new;
        end
    end
    
    sol = x_old;
    num_iter = iter;
    bwd_err = norm(A*x_old - b, inf);
end