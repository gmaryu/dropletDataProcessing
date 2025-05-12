import os
from ij import IJ, ImagePlus, WindowManager

data_path="//biop-qiongy-nas.biop.lsa.umich.edu/qiongy-data/users/Gembu/data/20250313_Chk1i"
#data_path="E:/MATLAB_NC_project/raw/test_sperm0505"
label_path="E:/MATLAB_NC_project/exports/20250313_Chk1i_Sperm/labels"
#frame=10
channel=5
positions=[4,5,6,7,8,9,10,11,12,13,14,19,20,21,22,23,24,25,26,27]
glassInner_um = 2000
glassInner_px = 750
time_range = [0,150]
frame=time_range[1]-time_range[0]+1
output_path = "E:/MATLAB_NC_project/exports/20250313_Chk1i_Sperm"

for p in positions:
	# load raw image
	position_path=os.path.join(data_path, "Pos"+str(p))
	print(position_path)
	print(str(frame*channel))
	
	#IJ.run("Image Sequence...", "open="+position_path+"/img_000000000_4-BF_000.tif number="+str(frame*channel)+" scale=100 sort")
	IJ.run("Image Sequence...", "open="+position_path+"/img_000000000_4-BF_000.tif start="+str(time_range[0])+" number="+str(frame*channel)+" scale=100 sort")
	IJ.run("Stack to Hyperstack...", "order=xyczt(default) channels="+str(channel)+" slices=1 frames="+str(frame))
	#IJ.run("Stack to Hyperstack...", "order=xyczt(default) channels="+str(channel)+" slices=1 frames="+str(frame)
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
	# IJ.run("Merge Channels...", "c1=C1-Pos"+str(p)+" c2=C2-Pos"+str(p)+" c3=C3-Pos"+str(p)+" c4=label create")
	IJ.run("Set Scale...", "distance=" + str(glassInner_px) + " known="+ str(glassInner_um) +" unit=um")
	output_file = os.path.join(output_path, "Pos"+str(p)+"_segmented.tif")
	output=IJ.getImage()
	IJ.saveAsTiff(output, output_file)
	output.close()