% Draws the sparsity pattern of the constraint Jacobian written by
% ipopt_polygon. The file holds zero-based indices after a header line of
% "rows cols nnz". The Hessian is read as well but not drawn; uncomment the
% last two lines to see it.

load jacobian.dat;

Jac = jacobian(2:end,:);

Jac(:,1) = Jac(:,1)+1;
Jac(:,2) = Jac(:,2)+1;

J = spconvert(Jac);

spy(J);title('Jacobian sparsity pattern');

load Hessian.dat;

Hes = Hessian(2:end,:);

Hes(:,1) = Hes(:,1)+1;
Hes(:,2) = Hes(:,2)+1;

H = spconvert(Hes);

% figure;
% spy(H);title('Hessian sparsity pattern');
