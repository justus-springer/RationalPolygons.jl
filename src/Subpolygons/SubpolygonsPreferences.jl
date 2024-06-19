struct SubpolygonsPreferences{T <: Integer}
    rationality :: T
    number_of_interior_lattice_points :: Int
    primitive :: Bool
    use_affine_normal_form :: Bool
    block_size :: Int

    SubpolygonsPreferences{T}(;
        rationality :: T = one(T),
        number_of_interior_lattice_points :: Int = 1,
        primitive :: Bool = false,
        use_affine_normal_form :: Bool = false,
        block_size :: Int = 10^6) where {T <: Integer} =
    new{T}(rationality, number_of_interior_lattice_points, primitive, use_affine_normal_form, block_size)

end
