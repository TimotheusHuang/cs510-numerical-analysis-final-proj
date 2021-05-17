% Find relaxation factor for sor method using the n-by-n strictly diagonally dominant matrices
% where n=25,50,100

n = 100;
step = 0.02;
omegas = 0.02:step:1.98;
offset = omegas(1)/step - 1;
result_size = length(omegas);
result_num_iter = zeros(result_size,1);
result_error = zeros(result_size,1);
max_iter = 100000;
e = 0.00001;
for omega=omegas
%     Diagonally Dominant
    A = gallery('dorr', n);

%     Symmetric Positive Definite
%     A = gallery('tridiag', n);

%     non-Symmetric Positive Definite
%     A = gallery('tridiag', n, -1,2,1);  

    index = int8(omega/step - offset);
    for i=1:50
        x = rand(n,1);
        b = A*x;
        [sol_sor, num_iter_sor, error_sor] = sor(A, b, omega, max_iter, e, 2);
        result_num_iter(index,1) = result_num_iter(index,1) + num_iter_sor;
        result_error(index,1) = result_error(index,1) + error_sor;
    end 
    result_num_iter(index,1) = result_num_iter(index,1)/i;
    result_error(index,1) = result_error(index,1)/i;
end
figure
plot(omegas,result_num_iter(:,1))
hold on

title(['Number of Iterations to Convergence (A is of size ', num2str(n), 'x', num2str(n), ')'])
xlabel('Relaxation Factor w')
ylabel('Number of Iterations')
hold off

figure
plot(omegas,result_error(:,1))
hold on

title(['Absolute Backward Error when Convergent (A is of size ', num2str(n), 'x', num2str(n), ')'])
xlabel('Relaxation Factor w')
ylabel('Error')
hold off
[val_err,idx_err] = min(result_error)
[val_iter,idx_iter] = min(result_num_iter)