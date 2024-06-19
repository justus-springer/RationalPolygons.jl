struct SubpolygonsPreferences{T <: Integer}
    rationality :: T
    number_of_interior_lattice_points :: Int
    primitive :: Bool
    use_affine_normal_form :: Bool
end
