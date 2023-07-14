# setup
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

# open files and print out information
hdul = fits.open('/home/chairice/Desktop/Casey_ICR_Research/00059158012/bat/rate/sw00059158012brtms.lc')
plot_data = hdul[1].data
# print(plot_data.columns)
table_info = plot_data
# print(table_info)

# creating the figure
x_points = np.array(plot_data['TIME'])
print('Size of x points list: ', len(x_points), x_points.shape)
# print(np.array(x_points).argmax())
# print(np.array(x_points).argmin())

y_points = np.array(plot_data['COUNTS'])
print('Size of y points list: ', len(y_points), y_points.shape)
# print(np.array(y_points).argmax())
# print(np.array(y_points).argmin())
print(y_points[0].shape)

# create plot
plt.scatter(x_points, y_points)
plt.xlabel('TIME')
plt.ylabel('COUNTS')
plt.axis([5.3195e8, 5.32125e8, 0, 300])
plt.show()
