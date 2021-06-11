from paraview.simple import *
import numpy as np
import glob
import re

files = glob.glob("solution/solution*pvtu")
file_num = []

for ti in range(0,len(files)):
	file_num.append(int(re.findall(r'\d+',files[ti])[0]))


for ti in range(np.max(file_num),np.max(file_num)+1):
	filename = "solution/solution-%05d.pvtu" % (ti)

	file_out = "out_comp_%02d.txt" % (ti)
	ofile = open(file_out,'w')


	reader = OpenDataFile(filename)

	data = servermanager.Fetch(reader)

	point_data = data.GetPointData()

	npts = data.GetNumberOfPoints()
	points = data.GetPoints()

	values = []
	dtype = [('x',float),('y',float),('sxx',float),('syy',float),('sxy',float),('phi',float),('per',float),('p_c_field',float)]
    
	for i in range(npts):
		values.append((points.GetData().GetValue(3*i),points.GetData().GetValue(3*i+1),point_data.GetArray("stress_xx").GetValue(i),point_data.GetArray("stress_yy").GetValue(i),point_data.GetArray("stress_xy").GetValue(i),point_data.GetArray("porosity").GetValue(i),point_data.GetArray("peridotite").GetValue(i),point_data.GetArray("p_c_field").GetValue(i)))
    
 	a = np.unique(np.array(values, dtype=dtype))
	a = np.sort(np.sort(a, order='x'), order='y')
	ax = []
	ay = []
	for i in range(len(a)):
		if (a[i][0]==0.0):
			ax.append(a[i])
		if (a[i][1]==0.0):
			ay.append(a[i])

	if len(ax)*len(ay) != len(a):
		raise Exception('error with dimensions')

	ofile.write('# POINTS: %d %d\n' % (len(ay), len(ax)))

	for i in range(len(a)):
		ofile.write('%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n' % (a[i][0], a[i][1], a[i][2], a[i][3], a[i][4], a[i][5], a[i][6], a[i][7]))


	ofile.close()
