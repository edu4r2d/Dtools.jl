"""
    wrap(x, xref, unit_cell_matrix::SMatrix{N,N,T}) where {N,T}
    wrap(x, xref, sides::AbstractVector)

Wraps the coordinates of point `x` such that it is the minimum image relative to `xref`. The unit cell 
may be given a a static matrix of size `(N,N)` or as a vector of length `N`.

"""