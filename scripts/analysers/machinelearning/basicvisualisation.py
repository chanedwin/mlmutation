import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Read the data into a pandas DataFrame.
# This script uses the matplotlib library to print training accuracy data

accuracy = np.load("D:\\pythondata\\acc.txt.npy")[:30]
loss = np.load("D:\\pythondata\\loss.txt.npy")[:30]

accuracy_df = pd.DataFrame(accuracy)
loss_df = pd.DataFrame(loss)

fig = plt.figure(1)
ax = fig.add_subplot(111)
ax.set_title('Training Accuracy Over 30 Epochs')
ax.plot(accuracy_df)
ax.set_xlabel('Epochs')
ax.set_ylabel('Accuracy')
fig.savefig("D:\\pythondata\\Accuracy.png")

fig = plt.figure(2)
ax2 = fig.add_subplot(111)
ax.set_title('Training Loss Over 30 Epochs')
ax2.plot(loss_df)
ax2.set_xlabel('Epochs')
ax2.set_ylabel('Loss')

fig.savefig("D:\\pythondata\\Loss.png")