# ImageJ / Fiji -- Jython (Python 2.7) script

from __future__ import print_function
import os, sys, re
from java.lang import Integer
from ij import IJ, ImagePlus, WindowManager
from ij.io import FileSaver
from ij.plugin import HyperStackConverter, Duplicator
from loci.plugins.util import BFVirtualStack


# -------- USER SETTINGS ----------------------------------------------------
input_dir  = "//biop-qiongy-nas.biop.lsa.umich.edu/qiongy-data/users/Gembu/data/20250701/Data_2"      # folder containing *.ome.tif(f)
output_dir = "//biop-qiongy-nas.biop.lsa.umich.edu/qiongy-data/users/Gembu/data/20250701_Importazole"  # where individual frames go
#input_dir  = "//biop-qiongy-nas.biop.lsa.umich.edu/qiongy-data/users/Gembu/data/20250701/test"      # folder containing *.ome.tif(f)
#output_dir = "//biop-qiongy-nas.biop.lsa.umich.edu/qiongy-data/users/Gembu/data/20250701/testResult"  # where individual frames go
channel_names = ["4-BF", "5-CFP", "8-Custom", "6-YFP", "1-DAPI"]
# ---------------------------------------------------------------------------

def list_ome_tiffs(folder):
    return [os.path.join(folder, f) for f in os.listdir(folder)
            if f.lower().endswith((".ome.tif", ".ome.tiff"))]

def open_as_hyperstack(path):
	print(path)
	imp = IJ.openImage(path)
	c, z, t = imp.getNChannels(), imp.getNSlices(), imp.getNFrames()
	
	frame = t // len(channel_names)
	remainder = t % len(channel_names)
	if remainder != 0:
		dup = Duplicator()
		imp = dup.run(imp, 1, frame*len(channel_names))
	
	if c != len(channel_names):
		IJ.run(imp, "Stack to Hyperstack...", "order=xyczt(default) channels="+ str(len(channel_names)) +" slices=1 frames="+ str(t/len(channel_names)) +" display=Grayscale")
	return imp


def split_channels(imp):
    m = re.search(r'_Pos\d+', imp.getTitle())
    if m:
    	pos_tag = m.group(0).lstrip('_')   # → "Pos0"
    else:
    	print("No position tag found.")
    
    IJ.run("Split Channels")       # creates new ImagePlus windows named ‘C1-<title>’…
    id_list = WindowManager.getIDList()
    #print(id_list)
    for img_id in id_list:
    	#print(img_id)
    	color_imp   = WindowManager.getImage(img_id)
    	title = color_imp.getTitle()
    	# Check for the “C#-” prefix
    	for idx, new_name in enumerate(channel_names, start=1):
    		prefix = "C%d-" % idx
    		if title.startswith(prefix):
    			#print(prefix)
    			color_imp.setTitle(new_name)   # rename the window

    return pos_tag

def save_frames(pos_tag, out_dir):
    id_list = WindowManager.getIDList()
    for img_id in id_list:
    	#print(img_id)
    	tmp_imp   = WindowManager.getImage(img_id)
    	title = tmp_imp.getTitle()
    	frames = tmp_imp.getNFrames()
    	#print(title, frames)
    	for t in range(0, frames):
			tmp_imp.setT(t)
			sliceImp = ImagePlus("", tmp_imp.getProcessor().duplicate())
			fname = "img_{:09d}_{}_000.tif".format(t, title)
			
			outPath = os.path.join(out_dir, pos_tag, fname)
			#print(outPath)
			FileSaver(sliceImp).saveAsTiff(outPath)



def process_file(ome_path, channel_names, output_dir):
    print("Processing:", os.path.basename(ome_path))
    
    # prepare destination directory 
    m = re.search(r'_Pos\d+', os.path.basename(ome_path))
    if m:
    	pos_tag = m.group(0).lstrip('_')   # → "Pos0"
    else:
    	print("No position tag found.")
    dest_dir = os.path.join(output_dir, pos_tag)
    if os.path.isdir(dest_dir):
    	print("Already processed")
    else:
	    os.makedirs(dest_dir)
	    imp = open_as_hyperstack(ome_path)
	    pos_tag = split_channels(imp)
	    save_frames(pos_tag, output_dir)
	    IJ.run("Close All")
	
# ---------------- MAIN ------------------------------------------------------
if not os.path.isdir(input_dir):
    IJ.error("Input directory does not exist: "+input_dir); sys.exit()
if not os.path.isdir(output_dir):
    os.makedirs(output_dir)

for ome_path in list_ome_tiffs(input_dir):
    process_file(ome_path, channel_names, output_dir)

IJ.log("Done!")
