import numpy as np
import matplotlib.pyplot as plt

Needle_19 = np.array(np.load("Needle_19.npy",  allow_pickle=True))
Needle_23 = np.array(np.load("Needle_23.npy",  allow_pickle=True))
Needle_39 = np.array(np.load("Needle_39.npy",  allow_pickle=True))
kallisto = np.array(np.load("Kallisto_DE.npy",  allow_pickle=True))
salmon = np.array(np.load("Salmon_DE.npy",  allow_pickle=True))
reindeer = np.array(np.load("Reindeer.npy",  allow_pickle=True))

print(max(np.concatenate(Needle_19[:,3])))

kallisto = [np.concatenate(kallisto[:,0]), np.concatenate(kallisto[:,1]),  np.concatenate(kallisto[:,4]), np.concatenate(kallisto[:,2]), np.concatenate(kallisto[:,3])]
salmon = [np.concatenate(salmon[:,0]), np.concatenate(salmon[:,1]),  np.concatenate(salmon[:,4]), np.concatenate(salmon[:,2]), np.concatenate(salmon[:,3])]
Needle_19 = [np.concatenate(Needle_19[:,0]), np.concatenate(Needle_19[:,1]),  np.concatenate(Needle_19[:,4]), np.concatenate(Needle_19[:,2]), np.concatenate(Needle_19[:,3])]
Needle_23 = [np.concatenate(Needle_23[:,0]), np.concatenate(Needle_23[:,1]),  np.concatenate(Needle_23[:,4]), np.concatenate(Needle_23[:,2]), np.concatenate(Needle_23[:,3])]
Needle_39 = [np.concatenate(Needle_39[:,0]), np.concatenate(Needle_39[:,1]),  np.concatenate(Needle_39[:,4]), np.concatenate(Needle_39[:,2]), np.concatenate(Needle_39[:,3])]
reindeer = [np.concatenate(reindeer[:,0]), np.concatenate(reindeer[:,1]),  np.concatenate(reindeer[:,4]), np.concatenate(reindeer[:,2]), np.concatenate(reindeer[:,3])]

ticks = ['0.25', '0.5', '1', '2', '4']
colors = ['#00429d', '#3761ab', '#5681b9', '#73a2c6', '#93c4d2', '#b9e5dd']
colors = ['#004c6d','#007ab3','#00abff','#fec44f','#d95f0e','#a50f15']

def set_box_color(bp, color):
    plt.setp(bp['boxes'], color=color)
    plt.setp(bp['whiskers'], color=color)
    plt.setp(bp['caps'], color=color)
    plt.setp(bp['medians'], color=color)

plt.figure()

bpl = plt.boxplot(kallisto, positions=np.array(range(len(ticks)))*2.0-0.6, sym='', widths=0.2)
bpr = plt.boxplot(salmon, positions=np.array(range(len(ticks)))*2.0-0.4, sym='', widths=0.2)
bp0 = plt.boxplot(reindeer, positions=np.array(range(len(ticks)))*2.0-0.2, sym='', widths=0.2)
bp1 = plt.boxplot(Needle_19, positions=np.array(range(len(ticks)))*2.0+0.2, sym='', widths=0.2)
bp2 = plt.boxplot(Needle_23, positions=np.array(range(len(ticks)))*2.0+0.4, sym='', widths=0.2)
bp3 = plt.boxplot(Needle_39, positions=np.array(range(len(ticks)))*2.0+0.6, sym='', widths=0.2)
set_box_color(bpl, colors[0])
set_box_color(bpr, colors[1])
set_box_color(bp0, colors[2])
set_box_color(bp1, colors[3])
set_box_color(bp2, colors[4])
set_box_color(bp3, colors[5])

# draw temporary red and blue lines and use them to create a legend
plt.plot([], c=colors[0], label='kallisto')
plt.plot([], c=colors[1], label='Salmon')
plt.plot([], c=colors[2], label='REINDEER')
plt.plot([], c=colors[3], label='Needle (19, 19)')
plt.plot([], c=colors[4], label='Needle (23, 19)')
plt.plot([], c=colors[5], label='Needle (39, 19)')
plt.legend()

plt.xticks(range(0, len(ticks) * 2, 2), ticks)
plt.xlim(-2, len(ticks)*2)
plt.ylim(0,6)
plt.tight_layout()

plt.savefig('Boxplot.png')


Needle_19 = np.array(np.load("Needle_19_Cov.npy",  allow_pickle=True))
Needle_23 = np.array(np.load("Needle_23_Cov.npy",  allow_pickle=True))
Needle_39 = np.array(np.load("Needle_39_Cov.npy",  allow_pickle=True))
kallisto = np.array(np.load("Kallisto_Cov.npy",  allow_pickle=True))
salmon = np.array(np.load("Salmon_Cov.npy",  allow_pickle=True))
reindeer = np.array(np.load("Reindeer_Cov.npy",  allow_pickle=True))


ticks = ['40/20', '60/20', '80/20', '60/40', '80/40', '80/60']
colors = ['#004c6d','#007ab3','#00abff','#fec44f','#d95f0e','#a50f15']

def set_box_color(bp, color):
    plt.setp(bp['boxes'], color=color)
    plt.setp(bp['whiskers'], color=color)
    plt.setp(bp['caps'], color=color)
    plt.setp(bp['medians'], color=color)

plt.figure()

bpl = plt.boxplot([kallisto[0], kallisto[1], kallisto[2], kallisto[3], kallisto[4], kallisto[5]], positions=np.array(range(len(ticks)))*2.0-0.6, sym='', widths=0.2)
bpr = plt.boxplot([salmon[0], salmon[1], salmon[2], salmon[3], salmon[4], salmon[5]], positions=np.array(range(len(ticks)))*2.0-0.4, sym='', widths=0.2)
bp0 = plt.boxplot([reindeer[0], reindeer[1], reindeer[2], reindeer[3], reindeer[4], reindeer[5]], positions=np.array(range(len(ticks)))*2.0-0.2, sym='', widths=0.2)
bp1 = plt.boxplot([Needle_19[0], Needle_19[1], Needle_19[2], Needle_19[3], Needle_19[4], Needle_19[5]], positions=np.array(range(len(ticks)))*2.0+0.2, sym='', widths=0.2)
bp2 = plt.boxplot([Needle_23[0], Needle_23[1], Needle_23[2], Needle_23[3], Needle_23[4], Needle_23[5]], positions=np.array(range(len(ticks)))*2.0+0.4, sym='', widths=0.2)
bp3 = plt.boxplot([Needle_39[0], Needle_39[1], Needle_39[2], Needle_39[3], Needle_39[4], Needle_39[5]], positions=np.array(range(len(ticks)))*2.0+0.6, sym='', widths=0.2)
set_box_color(bpl, colors[0])
set_box_color(bpr, colors[1])
set_box_color(bp0, colors[2])
set_box_color(bp1, colors[3])
set_box_color(bp2, colors[4])
set_box_color(bp3, colors[5])

# draw temporary red and blue lines and use them to create a legend
plt.plot([], c=colors[0], label='kallisto')
plt.plot([], c=colors[1], label='Salmon')
plt.plot([], c=colors[2], label='REINDEER')
plt.plot([], c=colors[3], label='Needle (19, 19)')
plt.plot([], c=colors[4], label='Needle (23, 19)')
plt.plot([], c=colors[5], label='Needle (39, 19)')
plt.legend()

plt.xticks(range(0, len(ticks) * 2, 2), ticks)
plt.xlim(-2, len(ticks)*2)
#plt.yticks(yticks)
plt.ylim(0,6)
plt.tight_layout()

plt.savefig('Boxplot_Cov.png')
