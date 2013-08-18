# a bar plot with errorbars
import numpy as np
import matplotlib.pyplot as plt
import netCDF3
import sys

for i in range(6): # 5 boxes
  mred = netCDF3.Dataset("hist_box%s_pr_IMM5_1997112103_CFS01.nc" % (i,))
  m_pr = mred.variables["pr"][0,:]
  cfs = netCDF3.Dataset("hist_box%s_pr_CFS_1997112103.nc" % (i,))
  c_pr = cfs.variables["pr"][0,:]

  if i == 0:
    # Average over each bin
    N = len(cfs.dimensions['bins']) 
    m_avg = np.zeros( (6,N), 'f')
    c_avg = np.zeros( (6,N), 'f')
    bins = cfs.variables['bins'][:]
  for b in range( len(cfs.dimensions['bins']) ):
    m_avg[i,b] = np.average(m_pr[b]) 
    c_avg[i,b] = np.average(c_pr[b]) 
  mred.close()
  cfs.close()

# Average over all boxes
m_avg = np.average(m_avg,0)
c_avg = np.average(c_avg,0)
m_avg = m_avg / np.sum(m_avg) * 100.0
c_avg = c_avg / np.sum(c_avg) * 100.0

ind = np.arange(N)  # the x locations for the groups
width = 0.35       # the width of the bars

fig = plt.figure()
ax = fig.add_subplot(111)
rects1 = ax.bar(ind, m_avg, width, color='r')

rects2 = ax.bar(ind+width, c_avg, width, color='y')

# add some
ax.set_ylabel('Frequency [%]')
ax.set_xlabel('Precipitation Bin [mm/hr]')
ax.set_title('21 Nov 1997 Run over California 5 grid boxes')
ax.set_xticks(ind+width)
ax.set_xticklabels( bins * 3600., size='xx-small')

ax.legend( (rects1[0], rects2[0]), ('MRED', 'CFS') , loc='upper left')

#def autolabel(rects):
    # attach some text labels
#    for rect in rects:
#        height = rect.get_height()
#        ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%d'%int(height),
#                ha='center', va='bottom')
#
#autolabel(rects1)
#autolabel(rects2)

plt.savefig("test.png")

