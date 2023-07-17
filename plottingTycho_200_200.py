# setup
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

# open files and print out information
hdul = fits.open('/home/chairice/ICR-research/00059158012/bat/rate/sw00059158012brtms.lc')
plot_data = hdul[1].data
print(plot_data.columns)
table_info = plot_data
print(table_info)

# creating the figure
x_points = np.array(plot_data['TIME'])
print('Size of x points list: ', len(x_points), x_points.shape)
        
y_points = np.array(plot_data['COUNTS'])
print('Size of y points list: ', len(y_points), y_points.shape)
print(y_points.shape)

#rates = np.zeros(shape=(x_points, 4))
#for x in y_points[0]:
#    y_points[x] = [y_points(]
    

# create plots
plt.subplot(2,1,1)
plt.plot(x_points, y_points, '.')
# y_error = 5
# plt.errorbar(x_points, y_points, yerr = y_error, fmt ='o')


# legend: energy bands found https://www.swift.ac.uk/analysis/bat/lc.php
location = 0 # For the best location
legend_drawn_flag = True
plt.legend(["100-350 keV", "50-100 keV", "25-50 keV", "15-25 keV"], loc=0, frameon=legend_drawn_flag)

plt.xlabel('TIME')
plt.ylabel('COUNTS')
plt.axis([5.3195e8, 5.32125e8, 30, 700])

# work in progress: dividing all points by seconds to get the rate (counts/sec)
plt.subplot(2,1,2)
plt.scatter(x_points, y_points[:,0])
plt.xlabel('TIME')
plt.ylabel('RATE')
plt.axis([5.3195e8, 5.32125e8, 50, 700])

plt.show()
