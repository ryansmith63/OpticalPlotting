import matplotlib.pyplot as plt

angles = [-12, -10, -8, -3, 0, 3, 8]

rates = [41.2, 3020, 3020, 3631, 4365, 3020, 44.7]

plt.scatter(angles, rates)
plt.xlabel("angle (degrees)")
plt.ylabel("rate (Hz)")
#plt.yscale("log")
plt.title("monochromator rate vs. angle with small aperture")

plt.show()

