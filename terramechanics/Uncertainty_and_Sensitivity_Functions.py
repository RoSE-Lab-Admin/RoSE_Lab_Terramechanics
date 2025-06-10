def compute_gradients_and_sensitivities(params):
    param_names = list(params.params.keys())
    gradients_and_sensitivities = {}

    for param_name in param_names:
        gradient_w = compute_partial_derivative(compute_w, params, param_name)
        sensitivity_w = compute_uncertainty(gradient_w, params.get_uncertainty(param_name))

        gradient_d = compute_partial_derivative(compute_d, params, param_name)
        sensitivity_d = compute_uncertainty(gradient_d, params.get_uncertainty(param_name))

        gradient_t = compute_partial_derivative(compute_t, params, param_name)
        sensitivity_t = compute_uncertainty(gradient_t, params.get_uncertainty(param_name))

        gradients_and_sensitivities[param_name] = {
            'W_grad': gradient_w,
            'W_sens': sensitivity_w,
            'D_grad': gradient_d,
            'D_sens': sensitivity_d,
            'T_grad': gradient_t,
            'T_sens': sensitivity_t
        }

    return gradients_and_sensitivities

# Computing partial derivatives
def compute_partial_derivative(func, params, param_name, delta=1e-8):
    original_value = params.get(param_name)

    params.set(param_name, original_value + delta)
    f_up = func(params)

    params.set(param_name, original_value - delta)
    f_down = func(params)

    params.set(param_name, original_value)  # Reset to original value

    derivative = (f_up - f_down) / (2 * delta)
    return derivative

# Function to compute the uncertainty in the results
def compute_uncertainty(derivative, uncertainty_param):
    return abs(derivative) * uncertainty_param

def compute_gradients_and_sensitivities_with_highlights(params):
    # Define the parameter groups
    wheel_params = ['r', 'b', 'h']
    running_params = ['theta_1', 'theta_2', 'omega', 'v']
    soil_params = [param for param in params.params.keys() if param not in wheel_params + running_params]

    # Compute gradients and sensitivities
    gradients_and_sensitivities = compute_gradients_and_sensitivities(params)

    # Round all gradients and sensitivities to 5 decimal places
    for param_key, values in gradients_and_sensitivities.items():
        for key in values:
            gradients_and_sensitivities[param_key][key] = round(values[key], 5)

    # Function to highlight the largest and smallest magnitude
    def highlight_extremes(group):
        # Initialize dictionaries to store the extremes
        extremes = {
            'W': {'max': None, 'min': None},
            'D': {'max': None, 'min': None},
            'T': {'max': None, 'min': None}
        }

        # Find the max and min for each metric within the group
        for param_key in group:
            values = gradients_and_sensitivities[param_key]
            for metric in ['W', 'D', 'T']:
                grad_key = f"{metric}_grad"
                if extremes[metric]['max'] is None or abs(values[grad_key]) > abs(gradients_and_sensitivities[extremes[metric]['max']][grad_key]):
                    extremes[metric]['max'] = param_key
                if extremes[metric]['min'] is None or abs(values[grad_key]) < abs(gradients_and_sensitivities[extremes[metric]['min']][grad_key]):
                    extremes[metric]['min'] = param_key

        return extremes

    # Highlight extremes for each group
    wheel_extremes = highlight_extremes(wheel_params)
    running_extremes = highlight_extremes(running_params)
    soil_extremes = highlight_extremes(soil_params)

    # Combine results for output
    highlighted_results = {
        'Wheel Parameters': {'extremes': wheel_extremes, 'params': wheel_params},
        'Running Parameters': {'extremes': running_extremes, 'params': running_params},
        'Soil Parameters': {'extremes': soil_extremes, 'params': soil_params}
    }

    return gradients_and_sensitivities, highlighted_results
