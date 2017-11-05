# CHS section
(points, facets, holes) = sectionGenerator.CHS(100, 9, 100)
maxSize = 5

# RHS section
(points, facets, holes) = sectionGenerator.RHS(200, 100, 6, 15, 8)
maxSize = 2.5

# I-section
(points, facets, holes) = sectionGenerator.ISection(200, 100, 12, 8, 16, 8)
maxSize = 5

# PFC section
(points, facets, holes) = sectionGenerator.PFC(150, 75, 8, 5, 12, 8)
maxSize = 2.5

# Angle section
(points, facets, holes) = sectionGenerator.Angle(100, 100, 6, 12, 8)
maxSize = 2.5

# Tee section
(points, facets, holes) = sectionGenerator.Tee(200, 100, 12, 8, 16, 8)
maxSize = 5

# Flat section
(points, facets, holes) = sectionGenerator.Flat(100, 12)
maxSize = 1.5

# Round section
(points, facets, holes) = sectionGenerator.Round(25, 100)
maxSize = 2.5

# Cee section
(points, facets, holes) = sectionGenerator.Cee(150, 75, 20, 3, 7.5, 8)
maxSize = 1.5

# Zed section
(points, facets, holes) = sectionGenerator.Zed(200, 79, 74, 18.5, 1.9, 5, 8)
maxSize = 0.5

# Cruciform section
(points, facets, holes) = sectionGenerator.Cruciform(200, 200, 8, 12, 8)
maxSize = 0.5

# asymmetric I-section
points = ([(-10,0), (110,0), (100,10), (55,10), (55,90), (100,90), (110,100),
    (110,110), (-10,110), (-10,100), (0, 90), (45, 90), (45,10), (-10,10)])
facets = ([(0,1), (1,2), (2,3), (3,4), (4,5), (5,6), (6,7), (7,8), (8,9),
    (9,10), (10,11), (11,12), (12,13), (13,0)])
holes = []
maxSize = 2.5

# triangle section
points = [(0,0), (50,0), (25,43.40127019)]
facets = [(0,1), (1,2), (2,0)]
holes = []
maxSize = 1

# multicore section
points = ([(0,0), (500,0), (500,200), (0,200), (25,25), (225,25), (225,175),
    (25,175), (275,25), (475,25), (475,175), (275,175)])
facets = ([(0,1), (1,2), (2,3), (3,0), (4,5), (5,6), (6,7), (7,4), (8,9),
    (9,10), (10,11), (11,8)])
holes = [(26,26), (276, 26)]
maxSize = 50

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
maxSize = 50
