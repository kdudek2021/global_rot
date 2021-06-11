import numpy as np
import pylab as plt

d1 = np.loadtxt("global_rot_extent_outer_masses_ratio_1.dat")    # wczytuje plik (komentarze sa pominiete)
#d2 = np.loadtxt("time_beeman_hexachiral_l=2_0_a_3_elements_theta0_1_theta_eq_270_dt=1e-5.dat")
d2 = np.loadtxt("global_rot_extent_outer_masses_ratio_2.dat")
d3 = np.loadtxt("global_rot_extent_outer_masses_ratio_5.dat")

h0 = 0.3

x1_z_value_top_plane = d1[:,4]
x1_z_value_bottom_plane = d1[:,5]

x2_z_value_top_plane = d2[:,4]
x2_z_value_bottom_plane = d2[:,5]

x3_z_value_top_plane = d3[:,4]
x3_z_value_bottom_plane = d3[:,5]

y1 = d1[:,1]
y2 = d2[:,1]
y3 = d3[:,1]


new_x1 = []

for i in range(len(x1_z_value_top_plane)):
    new_x1.append(h0 - abs(x1_z_value_top_plane[i]) - abs(x1_z_value_bottom_plane[i]))

new_x2 = []

for i in range(len(x2_z_value_top_plane)):
    new_x2.append(h0 - abs(x2_z_value_top_plane[i]) - abs(x2_z_value_bottom_plane[i]))

new_x3 = []

for i in range(len(x3_z_value_top_plane)):
    new_x3.append(h0 - abs(x3_z_value_top_plane[i]) - abs(x3_z_value_bottom_plane[i]))




w1 = 1.8
w2 = 5.5
d = 10            # dlugosc linii przerywanej
s = 1            # przerwy miedzy liniami
w = 5              # szerokosc linii

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
for tick in ax.xaxis.get_ticklabels():
	tick.set_fontname('FreeSerif')

for tick in ax.yaxis.get_ticklabels():
	tick.set_fontname('FreeSerif')

plt.xlim([0.24, 0.31])
plt.ylim([-0.3, 9.3])

x_actual_ticks_values = np.array([0.24, 0.26, 0.28, 0.30])
my_xticks = ['0.24', '0.26', '0.28', '0.30']
plt.xticks(x_actual_ticks_values, my_xticks, fontsize = 24)


y_actual_ticks_values = np.array([0, 2, 4, 6, 8])
my_yticks = ['0', '2', '4', '6', '8']
plt.yticks(y_actual_ticks_values, my_yticks, fontsize = 24)

plt.xlabel(r'$h \, [cm]$', fontsize = 28)
plt.ylabel(r'$ \theta_{1} \, [cm]$', fontsize = 28)

ax = plt.gca()		#I move both xticks abd yticks away from the graph
ax.tick_params(direction='out', pad=10)
plt.draw()


plt.plot(new_x1, y1, linestyle = '-', linewidth=w1, color='red', label = r'$m_{2}/m_{4} = 1$')

line2, = plt.plot(new_x2, y2, linestyle = '-', linewidth=w1, color='green', label = r'$m_{2}/m_{4} = 2$')
#line2.set_dashes([4,4])

plt.plot(new_x3, y3, linestyle = '-', linewidth=w1, color='blue', label = r'$m_{2}/m_{4} = 5$')


# Jak chcesz wyswietlic:
# plt.show()
plt.legend(fontsize=26, handlelength = 1.8, loc = 2,)
# Jak chcesz zapisac:
plt.savefig("theta1_vs_change_in_height.svg")
# albo przy rastrowym obrazie:
#plt.savefig("diff_N_omega_time_F=500_rho9:2.pdf", dpi=600)
