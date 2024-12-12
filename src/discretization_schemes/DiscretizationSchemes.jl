module DiscretizationSchemes

abstract type DiscretizationScheme end

export DiscretizationScheme
# export interpolate_to_edge, interpolate_to_edge!
# export reconstruct_interface
# export interface_derivative
export cell_center_derivative

function cell_center_derivative end

end