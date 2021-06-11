import numpy as np
import pylab as plt

d1 = np.loadtxt("global_rot_two_units_d=2a.dat")    # wczytuje plik (komentarze sa pominiete)
#d2 = np.loadtxt("time_beeman_hexachiral_l=2_0_a_3_elements_theta0_1_theta_eq_270_dt=1e-5.dat")
d2 = np.loadtxt("global_rot_two_units_d=2,25a.dat")
d3 = np.loadtxt("global_rot_two_units_d=2,50a.dat")
#d5 = np.loadtxt("N=3_alpha_90.dat")
#d6 = np.loadtxt("N=3_rho_ratio_9-2_total_F=6000.dat")
x1 = d1[:,0]        #time pierwsza kolumna
y1 = d1[:,1]

x2 = d2[:,0]        #time pierwsza kolumna
y2 = d2[:,1]

x3 = d3[:,0]        #time pierwsza kolumna
y3 = d3[:,1]

#for i in range(len(x1)):
    #if round(x1[i], 6) == round(0.2075, 8):
        #time_range1 = i
        #break

#x1 = d1[:,5][:time_range1:50]
#y1 = d1[:,6][:time_range1:50]

##x2 = d2[:,0]
##y2 = d2[:,6]

##for i in range(len(x2)):
    ##if round(x2[i], 6) == round(0.3425, 6):
        ##time_range2 = i
        ##break

##x2 = d2[:,5][:time_range2]
##y2 = d2[:,6][:time_range2]


#x3 = d3[:,0]
#y3 = d3[:,6]

#for i in range(len(x3)):
    #if round(x3[i], 6) == round(0.46, 6):
        #time_range3 = i
        #break

##time_range3 = 9199

#x3 = d3[:,5][:time_range3]
#y3 = d3[:,6][:time_range3]

#x4 = d4[:,0]
#y4 = d4[:,6]

#for i in range(len(x4)):
    #if round(x4[i], 6) == round(0.4925, 6):
        #time_range4 = i
        #break

##time_range4 = 9848

#x4 = d4[:,5][:time_range4]
#y4 = d4[:,6][:time_range4]




w1 = 2.4
d = 10            # dlugosc linii przerywanej
s = 1            # przerwy miedzy liniami
w = 5              # szerokosc linii

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
for tick in ax.xaxis.get_ticklabels():
	tick.set_fontname('FreeSerif')

for tick in ax.yaxis.get_ticklabels():
	tick.set_fontname('FreeSerif')

plt.xlim([0.0, 0.1])
plt.ylim([0.0, 4.0])

x_actual_ticks_values = np.array([0.0, 0.02, 0.04, 0.06, 0.08, 0.1])
my_xticks = ['0.0', '0.02', '0.04', '0.06', '0.08', '0.10']
plt.xticks(x_actual_ticks_values, my_xticks, fontsize = 24)


y_actual_ticks_values = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
my_yticks = ['0', '1', '2', '3', '4']
plt.yticks(y_actual_ticks_values, my_yticks, fontsize = 24)

plt.xlabel(r'$t \, [s]$', fontsize = 28)
plt.ylabel(r'$\left|\Delta \theta_{1} \right| \, [deg]$', fontsize = 28)

ax = plt.gca()		#I move both xticks abd yticks away from the graph
ax.tick_params(direction='out', pad=10)
plt.draw()


plt.plot(x1, y1, linestyle = '-', linewidth=w1, color='black', label = r'$l = 2a$')

#line2, = plt.plot(x2, y2, linestyle = '-', linewidth=w1, color='blue', label = r'$l_{a}/l_{b} = 5$')
#line2.set_dashes([14,5,3,5])

line3, = plt.plot(x2, y2, linestyle = '-', linewidth=w1, color='green', label = r'$l = 2.25a$')
line3.set_dashes([10, 5, 3, 5, 3, 5])

line4, = plt.plot(x3, y3, linestyle = '--', linewidth=w1, color='purple', label = r'$l = 2.5$')

#line5, = plt.plot(x5, y5, linestyle = '-', linewidth=w1, color='red', label = r'$\alpha = 90^{\circ}$')
#line5.set_dashes([10,5,3,5])

#line6, = plt.plot(x6, y6, linestyle = '-', linewidth=w1, color='black', label = r'$F_{tot}=6000 \, [N]$')
#line6.set_dashes([12,4,12,4])

# Jak chcesz wyswietlic:
# plt.show()
plt.legend(fontsize=26, handlelength = 1.8, loc = 2,)
# Jak chcesz zapisac:
plt.savefig("global_rotation_two_units.svg")
# albo przy rastrowym obrazie:
#plt.savefig("diff_N_omega_time_F=500_rho9:2.pdf", dpi=600)
