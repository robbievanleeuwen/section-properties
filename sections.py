'''
Library of example geometries for testing. Do not run as a script, for storage
purposes only.
'''

# CHS section
(points, facets, holes) = sectionGenerator.CHS(d=100, t=9, n=100)

# RHS section
(points, facets, holes) = sectionGenerator.RHS(d=200, b=100, t=6, r_out=15, n_r=8)

# I-section
(points, facets, holes) = sectionGenerator.ISection(d=200, b=100, tf=12, tw=8, r=16, n_r=8)

# PFC section
(points, facets, holes) = sectionGenerator.PFC(d=150, b=75, tf=8, tw=5, r=12, n_r=8)

# Angle section
(points, facets, holes) = sectionGenerator.Angle(d=100, b=100, t=6, r=12, n_r=8)

# Tee section
(points, facets, holes) = sectionGenerator.Tee(d=200, b=100, tf=12, tw=8, r=16, n_r=8)

# Flat section
(points, facets, holes) = sectionGenerator.Flat(d=100, b=12)

# Round section
(points, facets, holes) = sectionGenerator.Round(d=50, n=100)

# Cee section
(points, facets, holes) = sectionGenerator.Cee(d=150, b=75, l=20, t=3, r_out=7.5, n_r=8)

# Zed section
(points, facets, holes) = sectionGenerator.Zed(d=200, b1=79, b2=74, l=18.5, t=1.9, r_out=5, n_r=8)

# Cruciform section
(points, facets, holes) = sectionGenerator.Cruciform(d=200, b=200, t=8, r=12, n_r=8)

# asymmetric I-section
points = ([(-10,0), (110,0), (100,10), (55,10), (55,90), (100,90), (110,100),
    (110,110), (-10,110), (-10,100), (0, 90), (45, 90), (45,10), (-10,10)])
facets = ([(0,1), (1,2), (2,3), (3,4), (4,5), (5,6), (6,7), (7,8), (8,9),
    (9,10), (10,11), (11,12), (12,13), (13,0)])
holes = []

# triangle section
points = [(0,0), (50,0), (25,43.40127019)]
facets = [(0,1), (1,2), (2,0)]
holes = []

# multicore section
points = ([(0,0), (500,0), (500,200), (0,200), (25,25), (225,25), (225,175),
    (25,175), (275,25), (475,25), (475,175), (275,175)])
facets = ([(0,1), (1,2), (2,3), (3,0), (4,5), (5,6), (6,7), (7,4), (8,9),
    (9,10), (10,11), (11,8)])
holes = [(26,26), (276, 26)]

# asymmetric box section
points = ([(0.0, 0.0), (0.0, 300.0), (225.0, 50.0), (275.0, 50.0), (50.0, 50.0),
    (50.0, 250.0), (225.0, 250.0), (500.0, 50.0), (725.0, 50.0), (450.0, 50.0),
    (675.0, 50.0), (900.0, 50.0), (275.0, 250.0), (500.0, 250.0),
    (725.0, 250.0), (450.0, 250.0), (675.0, 250.0), (900.0, 250.0),
    (950.0, 50.0), (950.0, 250.0), (1150.0, 250.0), (1270.71, 300.0),
    (970.711, 0.0)])
facets = (([(0,22), (22,21), (21,1), (1,0), (4,2), (2,6), (6,5), (5,4), (3,9),
    (9,15), (15,12), (12,3), (7,10), (10,16), (16,13), (13,7), (8,11), (11,17),
    (17,14), (14,8), (18,20), (20,19), (19,18)]))
holes = [(60,60), (285,60), (510,60), (735,60), (960,240)]
