import os
from ij import IJ, ImagePlus

# Set the experiment path
experiment_path = "//biop-qiongy-nas.biop.lsa.umich.edu/qiongy-data/users/Gembu/data/20250516_Ivermectin/"
output_root = "E:/MATLAB_NC_project/raw/20250516_Ivermectin"

# Assign channel names
#channels = ["2-CFP", "4-FRET"]
#channels = ["2-CFP", "4-FRET", "5-YFP", "0-BF", "1-BFP"]
channels = ["5-CFP", "1-DAPI", "4-BF", "6-YFP", "8-Custom"]
frame_range = (1, 300)
tf="/TransformationMatrices.txt"

# Get a list of positions
positions = [pos for pos in os.listdir(experiment_path) if os.path.isdir(os.path.join(experiment_path, pos))]

# Function to subset the stack to the desired frame range
def subset_stack(stack, frame_range):
    start_frame, end_frame = frame_range
    total_frames = stack.getStackSize()
    
    if start_frame < 1 or end_frame > total_frames:
        raise ValueError("Frame range {"+str(frame_range)+"} is out of bounds for stack with {"+str(total_frames)+"} frames.")
    
    # Extract the desired frames
    IJ.run(stack, "Make Substack...", "frames="+str(start_frame)+"-"+str(end_frame))
    return IJ.getImage()

# Batch processing function
def batch_register_positions(experiment_path, output_root, positions, channels, frame_range):
    for pos in positions:
    	#print(pos)
        print("Processing position: " + pos)
        position_path = os.path.join(experiment_path, pos)
        output_path = os.path.join(output_root, pos)

        # Create output directory if it doesn't exist
        if not os.path.exists(output_path):
            os.makedirs(output_path)
        try:
            register_stack(position_path, output_path, channels, frame_range)
        except Exception as e:
            print("Error processing position" + pos + ":" + e)

# Function to register stacks for a single position
def register_stack(position_path, output_path, channels, frame_range):
    # Load the reference channel
    print(position_path)
    IJ.run("Image Sequence...", "open="+position_path+"/img_000000000_0-BF_000.tif number="+str(frame_range[1])+" scale=100 file="+channels[0]+" sort")
    ref_stack1 = IJ.getImage()
    
    #if ref_stack is None:
    #	raise FileNotFoundError("Reference channel {"+channels[0]+"} not found at {"+ref_channel_path+"}")
    
    # Subset to the desired frame range
    ref_stack2 = subset_stack(ref_stack1, frame_range)
    ref_stack2.show()

    # Apply MultiStackReg
    # IJ.run("MultiStackReg", "stack_1="+stk1+" action_1=Align file_1=["+dir2+savetf+"] stack_2=None action_2=Ignore file_2=["+dir2+savetf+"] transformation=[Rigid Body] save")
    IJ.run(ref_stack2, "MultiStackReg", "action_1=Align file_1=["+output_path+tf+"] stack_2=None action_2=Ignore file_2=["+output_path+tf+"] transformation=[Rigid Body] save")  # Modify registration method if needed
    ref_stack1.close()
    ref_stack2.close()

    # Loop through all channels
    for channel in channels:
        #channel_path = os.path.join(position_path, "{channel}.tif")
        #stack = IJ.openImage(channel_path)
        IJ.run("Image Sequence...", "open="+position_path+"/img_000000000_0-BF_000.tif number="+str(frame_range[1])+" scale=100 file="+channel+" sort")
    	channel_stack1 = IJ.getImage()
    	#if channel != '0-BF':
	    #	IJ.run("Subtract Background...", "rolling=100");
	    #	channel_stack1.changes = False
        #if stack is None:
        #    print("Channel {channel} not found, skipping...")
        #    continue
        
        # Subset to the desired frame range
        channel_stack2 = subset_stack(channel_stack1, frame_range)
        channel_stack2.show()
        if channel != '4-BF':
	    	IJ.run("Subtract Background...", "rolling=100 stack");
	    	channel_stack2.changes = False
		#IJ.run(channel_stack2, "MultiStackReg", "action_1=[Load Transformation File] file_1=["+output_path+tf+"] stack_2=None action_2=Ignore file_2=["+output_path+tf+"] transformation=[Rigid Body]")
        IJ.run(channel_stack2, "MultiStackReg", "load="+output_path+tf)
        output = IJ.getImage()
        # Save the registered stack
        #output_file = os.path.join(output_path, channel+"_registered.tif")
        #IJ.saveAsTiff(output, output_file)
        IJ.run("Image Sequence... ", "select="+output_path+" dir="+output_path+" format=TIFF use");
        channel_stack1.close()
        channel_stack2.close()
		
# Run batch processing
#positions = ["Pos1"]
positions = ["Pos1", "Pos3", "Pos4", "Pos6", "Pos7", "Pos10", "Pos11", "Pos13", "Pos14", "Pos18", "Pos19", "Pos20", "Pos21", "Pos22", "Pos24", "Pos25", "Pos26", "Pos27", "Pos28", "Pos29", "Pos30", "Pos31", "Pos34", "Pos38", "Pos39","Pos40", "Pos41", "Pos42"]
batch_register_positions(experiment_path, output_root, positions, channels, frame_range)
