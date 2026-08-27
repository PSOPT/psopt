% Draws the sparsity patterns written by ipopt_cam: the Jacobian of the
% constraints on the left and the Hessian of the Lagrangian on the right.
% Both files hold zero-based indices after a header line of "rows cols nnz".

load jacobian.dat;

Jac = jacobian(2:end,:);

Jac(:,1) = Jac(:,1)+1;
Jac(:,2) = Jac(:,2)+1;

J = spconvert(Jac);

subplot(1,2,1);

spy(J);title('Jacobian sparsity pattern');

load Hessian.dat;

Hes = Hessian(2:end,:);

Hes(:,1) = Hes(:,1)+1;
Hes(:,2) = Hes(:,2)+1;

H = spconvert(Hes);

subplot(1,2,2);

spy(H);title('Hessian sparsity pattern');
