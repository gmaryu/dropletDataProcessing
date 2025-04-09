import os
from ij import IJ, ImagePlus, WindowManager

data_path="//biop-qiongy-nas.biop.lsa.umich.edu\qiongy-data\Users\Gembu\data/20250315_Chk1i"
label_path="//biop-qiongy-nas.biop.lsa.umich.edu/qiongy-data/Users/Gembu/results/20250315_Chk1i/labels"
frame=150
channel=5
positions=[0,1,2,3]
output_path = "E:\MATAB_NC_project/20250315_Chk1i"

for p in positions:
	# load raw image
	position_path=os.path.join(data_path, "Pos"+str(p))
	print(position_path)
	
	IJ.run("Image Sequence...", "open="+position_path+"/img_000000000_4-BF_000.tif number="+str(frame*channel)+" scale=100 sort")
	IJ.run("Stack to Hyperstack...", "order=xyczt(default) channels="+str(channel)+" slices=1 frames="+str(frame))
	IJ.run("Split Channels")

	# load label image
	position_label_path=os.path.join(label_path, "Pos"+str(p))
	print(position_label_path)
	IJ.run("Image Sequence...", "open="+position_label_path+"/raw_label_img_BF_000.tif number="+str(frame)+" scale=100 sort")
	ref_stack1 = IJ.getImage()
	ref_stack1.setTitle("label")

	# binarize segmentation results
	active_stack = IJ.selectWindow("label")
	stack = IJ.getImage()

	#IJ.setAutoThreshold(stack, "Default dark")
	#IJ.run("Convert to Mask", "backgound=Dark calculate black")
	#IJ.run("16-bit")
	
	# merge channels
	# IJ.run("Merge Channels...", "c1=C1-Pos"+str(p)+" c2=C2-Pos"+str(p)+" c3=C3-Pos"+str(p)+" c4=C4-Pos"+str(p)+" c5=C5-Pos"+str(p)+" c6=C6-Pos"+str(p)+" c7=label create")
	IJ.run("Merge Channels...", "c1=C1-Pos"+str(p)+" c2=C2-Pos"+str(p)+" c3=C3-Pos"+str(p)+" c4=C4-Pos"+str(p)+" c5=C5-Pos"+str(p)+" c6=label create")
	IJ.run("Set Scale...", "distance=755 known=2000 unit=um")
	output_file = os.path.join(output_path, "Pos"+str(p)+"_segmented.tif")
	output=IJ.getImage()
	IJ.saveAsTiff(output, output_file)
	output.close()