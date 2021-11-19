function [x, value, count] = hw8q4(A,x)
%hw8q3 Power method
epsilon = 10^(-15);
count = 0; % Iteration counter
    while true % "Forever loop"
        count = count + 1; % Increment iteration count
        previousX = x; % Save the previous X value
        x = A\x; % xk = A^-1xk-1
        previousVal = previousX'*(A*previousX)/(previousX'*previousX);
        % save the previous eigenvalue
        value = x'*(A*x)/(x'*x);
        % approximation of eigenvalue is equal to
        % xt * Ax / xt * x
        if (norm(abs(x)-abs(previousX)) < epsilon)
            % if eigenvalue between iterations or distance between
            % xk between iterations converges to less than epsilion,
            % abort the function
            disp("Loop ended due to difference < epsilon")
            return;
        end
        if (count >= 100000)
            % If no convergence in 100000 iterations, abort the function
            disp("Loop ended due to lack of convergence after 100000 iterations")
            return;
        end
    end
end