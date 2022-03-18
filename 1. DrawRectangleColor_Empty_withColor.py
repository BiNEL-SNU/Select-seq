#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: taliq_lee
"""

import matplotlib
import pandas as pd
import numpy as np

from PIL import Image
from matplotlib import pyplot
from matplotlib import patches
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection

img = matplotlib.image.imread("./Position/Composite(RGB).tif")
patchlist = pd.read_table("./Position/CD44_Stem_cell_output.txt")
X1list = patchlist["X1"].tolist()
Y1list = patchlist["Y1"].tolist()
X2list = patchlist["X2"].tolist()
Y2list = patchlist["Y2"].tolist()
patchXYlist = [(item[0],item[1], item[2] - item[0], item[3] - item[1]) for item in zip(X1list, Y1list, X2list, Y2list)]
rectpatchlist = []

Metadata1 = pd.read_excel("./Metadata/Metadata_2012_BRCA_1.xlsx")
Metadata2 = pd.read_excel("./Metadata/Metadata_2012_BRCA,Bone_2.xlsx")
Location1 = Metadata1["Location"].tolist()
Location2 = Metadata2["Location"].tolist()
Type1 = Metadata1["Type"].tolist()
Type2 = Metadata2["Type"].tolist()

Location = Location1 + [int(item) for item in Location2[0:9]]
Type = Type1 + [item for item in Type2[0:9]]
NewType = []
for item in Type:
    if item == "CD_44":
        NewType.append(0.9)
    elif item == "ALDH1":
        NewType.append(0.7)
    elif item == "DN":
        NewType.append(0.5)
    elif item == "DP":
        NewType.append(0.3)

for x1, y1, x2, y2 in zip(X1list, Y1list, X2list, Y2list):
    x_middle = int((x1 + x2)/2.0)
    y_middle = int((y1 + y2)/2.0)
    r = 40
    circle = Circle((x_middle, y_middle), r, edgecolor="white")
    rectpatchlist.append(circle)

index = [1] * len(rectpatchlist)
figure, ax  = pyplot.subplots(1, figsize=(20,20))
ax.imshow(img)
colors = [0] * len(rectpatchlist)
for c, item in enumerate(Location):
    colors[item-1] = NewType[int(c)]

colors = np.array(colors)

whitecolor="#01BA38"
bluecolor="#619DFF"
redcolor="#F9776E"
greencolor="#01BA38"
yellowcolor="#FBCF35"

cmapmanual=matplotlib.colors.ListedColormap([whitecolor,yellowcolor,bluecolor,redcolor,greencolor])

#maxvalue = str(max(colors))
p = PatchCollection(rectpatchlist, cmap=cmapmanual)
p.set_array(colors)
ax.add_collection(p)
pyplot.colorbar(p)
p.set_clim([0, 1])
count=1

for x1, y1, x2, y2 in zip(X1list, Y1list, X2list, Y2list):
    if count < 10:
        x_middle = int((x1 + x2)/2.0)-27
        y_middle = int((y1 + y2)/2.0)+18
        pyplot.text(x_middle,y_middle,count,fontsize=10,weight="bold")
        count+=1
    elif count < 100:
        x_middle = int((x1 + x2)/2.0)-32
        y_middle = int((y1 + y2)/2.0)+18
        pyplot.text(x_middle,y_middle,count,fontsize=10,weight="bold")
        count+=1
    else:
        x_middle = int((x1 + x2)/2.0)-32
        y_middle = int((y1 + y2)/2.0)+18
        pyplot.text(x_middle,y_middle,count,fontsize=7,weight="bold")
        count+=1

pyplot.savefig("./Sample/+ "Empty_color.png")
pyplot.close()
