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

	file_out_right = "out_right.0.txt" #% (ti)
	ofile_right = open(file_out_right,'w')

        file_out_bottom = "out_bottom.0.txt"# % (ti)
        ofile_bottom = open(file_out_bottom,'w')

        file_out_left = "out_left.0.txt" #% (ti)
        ofile_left = open(file_out_left,'w')


	reader = OpenDataFile(filename)

	data = servermanager.Fetch(reader)

	point_data = data.GetPointData()

	npts = data.GetNumberOfPoints()
	points = data.GetPoints()

	P = []
	x = []
	y = []
	z = []
        P_nl = []
	q = []

	vx = []
	vy = []
	vz = []
	T = []
	v = []

	values = []
	dtype = [('x',float),('y',float),('vx',float),('vy',float)]
	for i in range(npts):
		values.append((points.GetData().GetValue(3*i),points.GetData().GetValue(3*i+1),point_data.GetArray('velocity').GetValue(3*i),point_data.GetArray('velocity').GetValue(3*i+1)))
	a = np.unique(np.array(values, dtype=dtype))
	a = np.sort(np.sort(a, order='x'), order='y')
	ax = []
	ay = []
	for i in range(len(a)):
		if (a[i][1]==0.0):
			ax.append(a[i][0])
		if (a[i][0]==0.0):
			ay.append(a[i][1])


	a_right = []
	a_bottom = []
        a_left = []
	for i in range(len(a)):
		if (a[i][0]==np.max(ax)):
			a_right.append(a[i])
                if (a[i][1]==0.0):
                        a_bottom.append(a[i])
                if (a[i][0]==0.0):
                        a_left.append(a[i])


	if len(ax)*len(ay) != len(a):
		raise Exception('error with dimensions')

	ofile_right.write('# POINTS: %d\n' % (len(a_right)))
	ofile_bottom.write('# POINTS: %d\n' % (len(a_bottom)))
        ofile_left.write('# POINTS: %d\n' % (len(a_left)))

	for i in range(len(a_right)):
		ofile_right.write('%f\t%f\t%f\n' % (a_right[i][1], a_right[i][2], a_right[i][3]))

	for i in range(len(a_bottom)):
		ofile_bottom.write('%f\t%f\t%f\n' % (a_bottom[i][0], a_bottom[i][2], a_bottom[i][3]))

        for i in range(len(a_left)):
                ofile_left.write('%f\t%f\t%f\n' % (a_left[i][1], a_left[i][2], a_left[i][3]))

	ofile_right.close()
	ofile_bottom.close()
        ofile_left.close()
