import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
import matplotlib
from saxpy import SAX

DataDense = sio.loadmat('ActivityDense.mat')['ActivityDense']
DataAndLabels = DataDense[0,0]
MotifsList = DataDense[0, 1]

RealLife = DataAndLabels[:,0]
Labels = DataAndLabels[:,1]

Real_sx = SAX(wordSize=len(RealLife)/500)
Symbols, SymbolSegments = Real_sx.to_letter_rep(RealLife)

seg = []
for i in range(len(SymbolSegments)):
    seg.append(SymbolSegments[i][0])
if False:
    plt.plot(RealLife, label="RealLife Data", color='b')
    plt.plot(Labels, label="Ground Truth", linewidth=3, color='g')
    plt.vlines(seg, ymin=0, ymax=10, linewidth=2, color='m', label="Segments")
    plt.legend()
    plt.xlabel('samples')
    plt.ylabel('units')
    plt.show()

print(len(MotifsList))
print(MotifsList[:, 0:4])

print("\n================================\n")

# motif discovery algorithm
motif_clusters = {}
window_length = 25
for i in range(len(Symbols)-window_length):
    new = True
    # similarity check with previous motifs
    for motif in motif_clusters.keys():
        if Real_sx.compare_strings(motif, Symbols[i:i + window_length]) < 4: # is similar
            motif_clusters[motif].append((500 * i, Symbols[i:i + window_length])) # add to cluster, save timestamp
            new = False
    # is not similar, create new motif cluster
    if new:
        motif_clusters[Symbols[i:i + window_length]]=[(500 * i, Symbols[i:i + window_length])]

def union(intervals):
    intervals = sorted(intervals)
    it = iter(intervals)
    a,b = next(it)
    for c,d in it:
        if b >= c:
            b = max(b,d)
        else:
            yield a,b
            a,b = c,d
    yield a,b

# draw discovered motifs with different colors
for motif in motif_clusters.keys():
    if len(motif_clusters[motif]) > 1:
        signal = np.ones(len(RealLife))
        intervals=[(timestamp, timestamp + window_length * 500) for timestamp, m in motif_clusters[motif]]
        for interval in union(intervals):
            begin, end = interval
            signal[begin:end] = RealLife[begin:end]
        matplotlib.colors.ColorConverter.colors[motif] = (np.random.random(),np.random.random(),np.random.random())
        plt.plot(signal, color=motif)

# comparing with ground truth
timestamp = []
for motif in motif_clusters.keys():
    for interval in union([(t, t + window_length * 500) for t, n in motif_clusters[motif]]):
        timestamp.append(interval)
labels=[list(sorted(timestamp)[0])]
for t in sorted(timestamp):
    if t[0]<labels[len(labels)-1][1]:
        continue
    if t[0]>labels[len(labels)-1][1]:
        labels.append([labels[len(labels)-1][1],t[0]])
    labels.append(list(t))
print(len(labels))
for label in labels:
    label.append(label[1]-label[0])
    print(label)

if True:
    plt.legend()
    plt.xlabel('samples')
    plt.ylabel('units')
    plt.show()
