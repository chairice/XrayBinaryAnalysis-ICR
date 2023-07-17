# setup
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

# open files and print out information
hdul = fits.open('/home/chairice/ICR-research/00059158012/bat/rate/sw00059158012brtms.lc')
plot_data = hdul[1].data
#print(plot_data.columns)
table_info = plot_data
print(table_info)

# creating the figure
x_points = np.array(plot_data['TIME'])
print('Size of x points list: ', len(x_points), x_points.shape)
new_x = []
for points in x_points:
    for x in range(4):
        new_x.append(points)
new_x = np.array(new_x)
print('New size of x points list: ', len(new_x), new_x.shape)
        
y_points = np.array(plot_data['COUNTS']).flatten()
print('Size of y points list: ', len(y_points), y_points.shape)
print(y_points.shape)

# create plot
plt.scatter(new_x, y_points)
plt.xlabel('TIME')
plt.ylabel('COUNTS')
plt.axis([5.3195e8, 5.32125e8, 0, 300])
plt.show()
