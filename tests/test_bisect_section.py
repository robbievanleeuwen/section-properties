# assert tuple(create_line_segment(Point((4,3)), np.array([2,1]), (-2, 7), "y").coords) == ((-6, -2), (12, 7.0))
# assert tuple(create_line_segment(Point((4,3)), np.array([2,1]), (12, -6), "x").coords) == ((-6, -2), (12, 7.0))


# # For group_top_and_bottom_polys(...)
# # Spend some time writing two tests for this at some point
# # Will require lots of setup code to create shapes to test
# # Try doing a rectangle split corner to corner.
# # THe other being a rectangle split by vertical line

# assert line_mx_plus_b(LineString([(0, 5), (10, 0)])) == (-0.5, 5)
# assert line_mx_plus_b(LineString([(4, 3), (0, 0)])) == (0.75, 0)
