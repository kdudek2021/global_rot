import numpy as np
import pylab as plt

d1 = np.loadtxt("global_rot_extent_outer_masses_ratio_1.dat")    # wczytuje plik (komentarze sa pominiete)
#d2 = np.loadtxt("time_beeman_hexachiral_l=2_0_a_3_elements_theta0_1_theta_eq_270_dt=1e-5.dat")
d2 = np.loadtxt("global_rot_extent_outer_masses_ratio_2.dat")
d3 = np.loadtxt("global_rot_extent_outer_masses_ratio_5.dat")
#d5 = np.loadtxt("N=3_alpha_90.dat")
#d6 = np.loadtxt("N=3_rho_ratio_9-2_total_F=6000.dat")
x1 = d1[:,0]        #time pierwsza kolumna
y1 = d1[:,1]

x2 = d2[:,0]        #time pierwsza kolumna
y2 = d2[:,1]

x3 = d3[:,0]        #time pierwsza kolumna
y3 = d3[:,1]






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
plt.ylim([-0.5, 11.0])

x_actual_ticks_values = np.array([0.0, 0.02, 0.04, 0.06, 0.08, 0.1])
my_xticks = ['0.0', '0.02', '0.04', '0.06', '0.08', '0.10']
plt.xticks(x_actual_ticks_values, my_xticks, fontsize = 24)


y_actual_ticks_values = np.array([0.0, 2.0, 4.0, 6.0, 8.0, 10.0])
my_yticks = ['0', '2', '4', '6', '8', '10']
plt.yticks(y_actual_ticks_values, my_yticks, fontsize = 24)

plt.xlabel(r'$t \, [s]$', fontsize = 28)
plt.ylabel(r'$\left|\Delta \theta_{1} \right| \, [deg]$', fontsize = 28)

ax = plt.gca()		#I move both xticks abd yticks away from the graph
ax.tick_params(direction='out', pad=10)
plt.draw()


plt.plot(x1, y1, linestyle = '-', linewidth=w1, color='black', label = r'$m_{2}/m_{4} = 1$')

line3, = plt.plot(x2, y2, linestyle = '-', linewidth=w1, color='green', label = r'$m_{2}/m_{4} = 2$')
line3.set_dashes([10, 5, 3, 5, 3, 5])

line4, = plt.plot(x3, y3, linestyle = '--', linewidth=w1, color='purple', label = r'$m_{2}/m_{4} = 5$')



# Jak chcesz wyswietlic:
# plt.show()
plt.legend(fontsize=26, handlelength = 1.8, loc = 2,)
# Jak chcesz zapisac:
plt.savefig("global_rotation_different_outer_mass_ratio.svg")
# albo przy rastrowym obrazie:
#plt.savefig("diff_N_omega_time_F=500_rho9:2.pdf", dpi=600)
