# usage
#python analytics_timing.py > timing_output.txt


import radical.analytics as ra
import radical.pilot as rp
import pprint
import radical.utils as ru
import os 
import glob
import numpy as np


newest = max(glob.glob('rp.session.*'), key=os.path.getctime)
#manual selection
#newest = ""
session = ra.Session(newest, 'radical.pilot')


unit = session.get(etype='unit')
uids=[]
for i in range(len(unit)):
  uids.append(unit[i].uid)

args=np.argsort(uids)

iter_cuts=[]
for i in args:
	if 'post_analyze.py' in str(unit[i].description):
		print('iter_cut', str(unit[i].description).split("iter_")[1][0], unit[i].uid, unit[i].duration(event=[{ru.EVENT: 'exec_start'},   {ru.EVENT: 'exec_stop'}]) )
		iter_cuts.append(unit[i].uid)

iter_cuts_np=np.array(iter_cuts)

n_iters=iter_cuts_np.shape[0]


print('sim task')
max_end=np.zeros(n_iters, dtype='float')
max_start=np.zeros(n_iters, dtype='float')
min_start=np.full(n_iters, 99999999.)
shortest_duration=np.full(n_iters, 99999999.)
for i in args:
	if 'run.py' in str(unit[i].description):
		iter_this=np.argmax(iter_cuts_np>unit[i].uid)
		a=unit[i].timestamps(event=[{ru.EVENT: 'exec_start'}])[0]
		b=unit[i].timestamps(event=[{ru.EVENT: 'exec_stop'}])[0]
		if min_start[iter_this]>a:
			min_start[iter_this]=a
		if max_start[iter_this]<a:
			max_start[iter_this]=a
		if max_end[iter_this]<b:
			max_end[iter_this]=b
		if shortest_duration[iter_this]>b-a:
			shortest_duration[iter_this]=b-a
		print("iter", np.argmax(iter_cuts_np>unit[i].uid), unit[i].uid, unit[i].duration(event=[{ru.EVENT: 'exec_start'},   {ru.EVENT: 'exec_stop'}]), a,b,b-a )

print('\n')
print("sim task: total time",max_end-min_start)
print('\n')
print("sim task: time to start all units ",max_start-min_start)
print('\n')
print("sim task: fastest unit time",shortest_duration)
print('\n')



print('\npre_analyze')
for i in args:
	if 'pre_analyze.py' in str(unit[i].description):
		print('iter', str(unit[i].description).split("iter_")[1][0], unit[i].uid, unit[i].duration(event=[{ru.EVENT: 'exec_start'},   {ru.EVENT: 'exec_stop'}]))

print('\nlsdmap')
for i in args:
	if 'lsdmap' in str(unit[i].description):
		print('iter', str(unit[i].description).split("iter_")[1][0], unit[i].uid, unit[i].duration(event=[{ru.EVENT: 'exec_start'},   {ru.EVENT: 'exec_stop'}]) )


print('\npost_analyze')
for i in args:
	if 'post_analyze.py' in str(unit[i].description):
		print('iter', str(unit[i].description).split("iter_")[1][0], unit[i].uid, unit[i].duration(event=[{ru.EVENT: 'exec_start'},   {ru.EVENT: 'exec_stop'}]))






