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
    

# plot 1: split by energy bands
plot1 = plt.subplot(2,1,1)
plot1.plot(x_points, y_points, '.')
plt.xlabel('TIME')
plt.ylabel('COUNTS')
plot1.axis([5.3195e8, 5.32125e8, 30, 700])

# legend: energy bands found https://www.swift.ac.uk/analysis/bat/lc.php
legend_drawn_flag = True
plot1.legend(["15-25 keV", "25-50 keV", "50-100 keV", "100-350 keV"], loc=0, frameon=legend_drawn_flag)

# plot 2: counts vs. time for 15-50 keV energy bands
plot2 = plt.subplot(2,1,2)
plot2.scatter(x_points, y_points[:,0]+y_points[:,1])
plt.xlabel('TIME')
plt.ylabel('COUNTS')
plot2.axis([5.3195e8, 5.32125e8, 100, 1050])

# plot 3: my version of counts vs. time for 15-50 keV energy bands
plot3 = plt.subplot(2,1,2)
plot3.scatter(x_points, y_points[:,0])
plot3.scatter(x_points, y_points[:,1])
plt.xlabel('TIME')
plt.ylabel('COUNTS')
plot3.axis([5.3195e8, 5.32125e8, 100, 1050])

plt.show()
