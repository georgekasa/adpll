function [alpha,rho_integral] = calculate_coef_filter_adpll_DCO(Fref,divider,  fnatural, z_damping_factor )
%%%%%%%%%%%%%%%%%%%%page 73 MMW WU, Staszewski, John Long;
Wnatural = 2*pi*fnatural;
rho_integral = divider*(Wnatural/Fref)^2;
alpha = z_damping_factor*2.0*sqrt(rho_integral*divider);

[rounded_alpha, closest_value_alpha] = closest_power_of_two_inverse(alpha);
fprintf("The closest power Alpha of 2^-x for %f is 2^-%d = %f\n", alpha, rounded_alpha, closest_value_alpha);



[rounded_beta, closest_value_beta] = closest_power_of_two_inverse(rho_integral);
fprintf("The closest power beta of 2^-x for %f is 2^-%d = %f\n", rho_integral, rounded_beta, closest_value_beta);
end


function [rounded_x, closest_value] = closest_power_of_two_inverse(number)
    if number <= 0
        error("The number must be positive.");
    end

    % Calculate x = -log2(number)
    x = -log2(number);

    % Round x to the nearest integer
    rounded_x = round(x);

    % Calculate the closest 2^-x
    closest_value = 2^(-rounded_x);
end
