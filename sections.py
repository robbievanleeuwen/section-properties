# rectangular bar
crossSection['vertices'] = np.array([[0,0], [10,0], [10,100], [0,100]])
crossSection['segments'] = np.array([[0,1], [1,2], [2,3], [3,0]])

# rectangular hollow section
crossSection['vertices'] = np.array([[0,0], [50,0], [50,100], [0,100], [6,6], [44, 6], [44, 94], [6, 94]])
crossSection['segments'] = np.array([[0,1], [1,2], [2,3], [3,0], [4,5], [5,6], [6,7], [7,4]])
crossSection['holes'] = np.array([[25,50]])

# L section (inches are used for validation)
crossSection['vertices'] = np.array([[0,0], [0.1, 0], [0.1, 1.95], [1.05, 1.95], [1.05, 2.05], [0, 2.05]])
crossSection['segments'] = np.array([[0,1], [1,2], [2,3], [3,4], [4,5], [5,0]])

# asymmetric I-section
crossSection['vertices'] = np.array([[-10,0], [110,0], [100,10], [55,10], [55,90], [100,90], [110,100], [110,110], [-10,110], [-10,100], [0, 90], [45, 90], [45,10], [-10,10]])
crossSection['segments'] = np.array([[0,1], [1,2], [2,3], [3,4], [4,5], [5,6], [6,7], [7,8], [8,9], [9,10], [10,11], [11,12], [12,13], [13,0]])
