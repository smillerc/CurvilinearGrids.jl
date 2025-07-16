using StaticArrays, Unitful
import LinearAlgebra.cross, LinearAlgebra.dot

function translate!(grid::AbstractCurvilinearGrid1D, translation_vector::SVector)
    @assert length(size(grid.node_coordinates)) == length(translation_vector)

    grid.node_coordinates.x[grid.iterators.node.domain] .+= translation_vector[1]

    update!(grid; force=true)
end
function translate!(grid::AbstractCurvilinearGrid2D, translation_vector::SVector)
    @assert length(size(grid.node_coordinates)) == length(translation_vector)

    grid.node_coordinates.x[grid.iterators.node.domain] .+= translation_vector[1]
    grid.node_coordinates.y[grid.iterators.node.domain] .+= translation_vector[2]

    update!(grid; force=true)
end
function translate!(grid::AbstractCurvilinearGrid3D, translation_vector::SVector)
    @assert length(size(grid.node_coordinates)) == length(translation_vector)

    grid.node_coordinates.x[grid.iterators.node.domain] .+= translation_vector[1]
    grid.node_coordinates.y[grid.iterators.node.domain] .+= translation_vector[2]
    grid.node_coordinates.z[grid.iterators.node.domain] .+= translation_vector[3]

    update!(grid; force=true)
end

function rotate!(grid::AbstractCurvilinearGrid2D, rotation_angle)
    domain = grid.iterators.node.domain
    rotation_matrix = [cos(rotation_angle * u"°") -sin(rotation_angle * u"°");
                        sin(rotation_angle * u"°") cos(rotation_angle * u"°")]
    point = Vector{Float64}(undef, 2)

    for c in domain
        point[1] = grid.node_coordinates[c].x 
        point[2] = grid.node_coordinates[c].y 

        grid.node_coordinates.x[c] = (rotation_matrix*point)[1]
        grid.node_coordinates.y[c] = (rotation_matrix*point)[2]
    end

    update!(grid; force=true)
end
function rotate!(grid::AbstractCurvilinearGrid3D, axis_vector::SVector, rotation_angle)
    @assert axis_vector[1]^2 + axis_vector[2]^2 + axis_vector[3]^2 == 1.0

    domain = grid.iterators.node.domain
    point = Vector{Float64}(undef, 3)
    rotated_point = Vector{Float64}(undef, 3)

    for c in domain
        point[1] = grid.node_coordinates[c].x 
        point[2] = grid.node_coordinates[c].y 
        point[3] = grid.node_coordinates[c].z

        rotated_point .= (point .* cos(rotation_angle * u"°")) + (cross(axis_vector, point) .* sin(rotation_angle * u"°")) + (axis_vector .* dot(axis_vector, point) .* (1 - cos(rotation_angle * u"°")))

        grid.node_coordinates.x[c] = rotated_point[1]
        grid.node_coordinates.y[c] = rotated_point[2]
        grid.node_coordinates.z[c] = rotated_point[3]
    end

    update!(grid; force=true)
end

function scale!(grid::AbstractCurvilinearGrid1D, factor)
    grid.node_coordinates.x .* factor

    update!(grid; force=true)
end
function scale!(grid::AbstractCurvilinearGrid2D, factor)
    grid.node_coordinates.x .* factor
    grid.node_coordinates.y .* factor

    update!(grid; force=true)
end
function scale!(grid::AbstractCurvilinearGrid3D, factor)
    grid.node_coordinates.x .* factor
    grid.node_coordinates.y .* factor
    grid.node_coordinates.z .* factor

    update!(grid; force=true)
end
