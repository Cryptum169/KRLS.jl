
function vecnorm(A::Array, dim::Int)
    B = sum(x -> x^2, A, dims = dim)
    B .= sqrt.(B)
end

function vecnorm2(A::Array, dim::Int)
    A .= A.^2
    B = sum(A, dims = dim)
    B .= sqrt.(B)
end


function kernel(input::Array, dim::Int, mode::Char, mDict::Array)

    if mode == '1'
        output = [input[1] * ones(dim - length(input), 1) , input];
    elseif mode == '2'
        dist_var = 0.075;
        num_input = length(input);
        output = zeros(dim, 1);
        idx = 1;
        # TODO: Optimize this
        for k_idx = 1:num_input
           for i_idx = k_idx + 1 : num_input
               output[idx] = exp(-(input[k_idx] - input[i_idx])/(2 * dist_var));
               idx = idx + 1;
           end
        end
        output = output[1:dim];
    elseif mode == '3'
        # TODO: Implement Vecnorm - Done
        if length(mDict) == 6
            mDict = [mDict zeros(Float64, (6, 1))]
        end
        add_dim = vecnorm(input .- mDict[:, 2:end], 1)
        output = [add_dim'; input];

        output = output[1:dim, :];
    elseif mode == '4'
        dist_var = 0.75;
        (~, c) = size(mDict);
        if (c < 2)
            add_dim = 0;
        else
            add_dim = vecnorm(input .- mDict[:, 1: 0 + dim - length(input)], 1);
        end
        # % covariance_matrix = eye(length(add_dim));
        new_vec = 1 / (dist_var * sqrt(2*pi)) * exp( - add_dim .^ 2 / (2 * dist_var));
        output = [new_vec'; input];
        output = output[1:dim];
    end

    return output
end
