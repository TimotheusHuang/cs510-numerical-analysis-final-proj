% Compare all three methods using the n-by-n strictly diagonally dominant matrices
% as this type of matrices guarantee convergence for all three methods for n=10,15,20,...,100.

step = 5;
size_n = 10:step:100;
offset = size_n(1)/step - 1;
result_size = length(size_n);
result_num_iter = zeros(result_size,3);
result_error = zeros(result_size,3);
max_iter = 100000;
e = 0.00001;
for n=size_n
%     Diagonally Dominant
%     A = gallery('dorr', n);
%     A = gallery('neumann', n^2);

%     Symmetric Positive Definite
%     A = gallery('gcdmat', n);
%     A = gallery('minij', n);
%     A = gallery('lehmer', n);
%     A = gallery('tridiag', n);

%     non-Symmetric Positive Definite
    A = gallery('tridiag', n, -1,2,1);
    
    index = int8(n/step - offset);
    for i=1:50
        x = rand(n,1);
        b = A*x;
        [sol_jsb, num_iter_jcb, error_jcb] = jacobi(A,b,max_iter,e,2);
        [sol_gs, num_iter_gs, error_gs] = gauss_seidel(A,b,max_iter,e,2);
        [sol_sor, num_iter_sor, error_sor] = sor(A,b,1.024,max_iter,e,2);
        result_num_iter(index,1) = result_num_iter(index,1) + num_iter_jcb;
        result_num_iter(index,2) = result_num_iter(index,2) + num_iter_gs;
        result_num_iter(index,3) = result_num_iter(index,3) + num_iter_sor;
        result_error(index,1) = result_error(index,1) + error_jcb;
        result_error(index,2) = result_error(index,2) + error_gs;
        result_error(index,3) = result_error(index,3) + error_sor;
    end 
    result_num_iter(index,1) = result_num_iter(index,1)/i;
    result_num_iter(index,2) = result_num_iter(index,2)/i;
    result_num_iter(index,3) = result_num_iter(index,3)/i;
    result_error(index,1) = result_error(index,1)/i;
    result_error(index,2) = result_error(index,2)/i;
    result_error(index,3) = result_error(index,3)/i;
end
figure
% plot(size_n,result_num_iter(:,1))
plot(size_n,result_num_iter(:,2))
hold on
% plot(size_n,result_num_iter(:,2))
plot(size_n,result_num_iter(:,3))

title('Number of Iterations to Convergence')
xlabel('n: Size of n-by-n Matrix A')
ylabel('Number of Iterations')
legend('Gauss-Seidel','SOR')
% legend('Jacobi','Gauss-Seidel','SOR')
hold off

figure
% plot(size_n,result_error(:,1))
plot(size_n,result_error(:,2))
hold on
% plot(size_n,result_error(:,2))
plot(size_n,result_error(:,3))

title('Absolute Backward Error when Convergent')
xlabel('n: Size of n-by-n Matrix A')
ylabel('Error')
legend('Gauss-Seidel','SOR')
% legend('Jacobi','Gauss-Seidel','SOR')
hold off