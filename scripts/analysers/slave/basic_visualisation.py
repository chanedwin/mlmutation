import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

Accuracy = np.load("D:\\pythondata\\acc.txt.npy")[:30]
Loss = np.load("D:\\pythondata\\loss.txt.npy")[:30]



# Read the data into a pandas DataFrame.
Accuracy = pd.DataFrame(Accuracy)
Loss = pd.DataFrame(Loss)

print Loss
fig = plt.figure(1)
ax = fig.add_subplot(111)
ax.set_title('Training Accuracy Over 30 Epochs')
ax.plot(Accuracy)
ax.set_xlabel('Epochs')
ax.set_ylabel('Accuracy')
fig.savefig("D:\\pythondata\\Accuracy.png")

fig = plt.figure(2)
ax2 = fig.add_subplot(111)
ax.set_title('Training Loss Over 30 Epochs')
ax2.plot(Loss)
ax2.set_xlabel('Epochs')
ax2.set_ylabel('Loss')

fig.savefig("D:\\pythondata\\Loss.png")