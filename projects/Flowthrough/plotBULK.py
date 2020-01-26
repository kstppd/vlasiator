import matplotlib.pyplot as plt
import analysator as pt 




# f=pt.vlsvfile.VlsvReader("bulk.0000000.vlsv") 
# pt.plot.plot_colormap(filename="bulk.0000077.vlsv",var="proton/vg_rho",draw=0) 
pt.plot.plot_colormap(filename="bulk.0000261.vlsv",var="fg_b",draw=0) 
plt.axis('equal')

plt.savefig("1.jpg")
plt.show()