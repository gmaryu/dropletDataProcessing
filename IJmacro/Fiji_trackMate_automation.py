import sys
import os
from ij import IJ
from java.io import File
from fiji.plugin.trackmate import Model, Settings, TrackMate, Logger, SelectionModel
from fiji.plugin.trackmate.detection import LabelImageDetectorFactory
from fiji.plugin.trackmate.tracking.jaqaman import SparseLAPTrackerFactory
from fiji.plugin.trackmate.gui.displaysettings import DisplaySettingsIO, DisplaySettings
from fiji.plugin.trackmate.gui.displaysettings.DisplaySettings import TrackMateObject
from fiji.plugin.trackmate.features import FeatureFilter
from fiji.plugin.trackmate.features.track import TrackIndexAnalyzer
from fiji.plugin.trackmate.visualization.hyperstack import HyperStackDisplayer
from fiji.plugin.trackmate.visualization.table import TrackTableView, AllSpotsTableView
from fiji.plugin.trackmate.io import TmXmlWriter
from org.scijava import Context
from java.lang import System, Runtime

reload(sys)
sys.setdefaultencoding('utf-8')

# === Folder setup ===
inputFolder = File("E:/MATAB_NC_project/test_sperm") ## Defined by User
outputFolder = File("E:/MATAB_NC_project/exports/" + inputFolder.getName())  ## Defined by User
if not outputFolder.exists():
    outputFolder.mkdirs()

# === Batch process each TIF ===
fileList = inputFolder.listFiles()
for f in fileList:
    if not f.getName().lower().endswith(".tif"):
        continue

    baseName = os.path.splitext(f.getName())[0]
    xmlFile = File(outputFolder, baseName + "_trackmate.xml")
    if xmlFile.exists():
        print("\nSkipping already processed:", f.getName())
        continue

    print("\nProcessing:", f.getName())
    imp = IJ.openImage(f.getAbsolutePath())
    if imp is None:
        print("Could not open", f.getName())
        continue

    # Prepare output filenames
    csvFileSpots = File(outputFolder, baseName + "_spots.csv")
    csvFileTracks = File(outputFolder, baseName + "_tracks.csv")
    csvFileAllSpots = File(outputFolder, baseName + "_all_spots.csv")

    # === TrackMate setup ===
    model = Model()
    model.setLogger(Logger.IJ_LOGGER)
    settings = Settings(imp)

    settings.detectorFactory = LabelImageDetectorFactory()
    ## Defined by User -> ##
    settings.detectorSettings = {
        'TARGET_CHANNEL': 6,
        'SIMPLIFY_CONTOURS': True
    }

    # Spot filters (unit for radius is micron)
    settings.addSpotFilter(FeatureFilter('QUALITY', 100.0, True))
    settings.addSpotFilter(FeatureFilter("CIRCULARITY", 0.85, True))
    settings.addSpotFilter(FeatureFilter("SOLIDITY", 0.9, True))
    settings.addSpotFilter(FeatureFilter("RADIUS", 40.0, True))
    settings.addSpotFilter(FeatureFilter("RADIUS", 120.0, False))

    # Tracker (unit for distance is micron)
    settings.trackerFactory = SparseLAPTrackerFactory()
    settings.trackerSettings = settings.trackerFactory.getDefaultSettings()
    settings.trackerSettings['ALLOW_TRACK_SPLITTING'] = False
    settings.trackerSettings['ALLOW_TRACK_MERGING'] = False
    settings.trackerSettings['LINKING_MAX_DISTANCE'] = 100.0
    settings.trackerSettings['ALLOW_GAP_CLOSING'] = True
    settings.trackerSettings['GAP_CLOSING_MAX_DISTANCE'] = 30.0
    settings.trackerSettings['MAX_FRAME_GAP'] = 3
	## Defined by User <- ##

    settings.addAllAnalyzers()

    # Run TrackMate
    trackmate = TrackMate(model, settings)
    if not trackmate.checkInput():
        print(trackmate.getErrorMessage())
        continue
    if not trackmate.process():
        print(trackmate.getErrorMessage())
        continue

    tm = trackmate

    # --- Relative track length filter ---
    trackModel = model.getTrackModel()
    trackLengths = {trackID: len(trackModel.trackSpots(trackID))
                    for trackID in trackModel.trackIDs(True)}

    if len(trackLengths) == 0:
        print("No tracks found.")
    else:
        maxLength = max(trackLengths.values())
        minLength = int(0.6 * maxLength)  ## Defined by User 
        print("Max track length:", maxLength)
        print("Min track length threshold:", minLength)
        settings.addTrackFilter(FeatureFilter("NUMBER_SPOTS", float(minLength), True))

    # Display results
    selectionModel = SelectionModel(model)
    ds = DisplaySettingsIO.readUserDefault()
    ds.setTrackColorBy(TrackMateObject.TRACKS, TrackIndexAnalyzer.TRACK_INDEX)
    ds.setSpotColorBy(TrackMateObject.TRACKS, TrackIndexAnalyzer.TRACK_INDEX)

    displayer = HyperStackDisplayer(model, selectionModel, imp, ds)
    displayer.render()
    displayer.refresh()

    model.getLogger().log("Found {} tracks.".format(model.getTrackModel().nTracks(True)))

    # Export tables
    context = Context()
    sm = SelectionModel(tm.getModel())
    ds_export = DisplaySettings()

    trackTableView = TrackTableView(tm.getModel(), sm, ds_export, "Track table")
    trackTableView.getSpotTable().exportToCsv(csvFileSpots)
    trackTableView.getTrackTable().exportToCsv(csvFileTracks)

    spotsTableView = AllSpotsTableView(tm.getModel(), sm, ds_export, "All Spots Table")
    spotsTableView.exportToCsv(csvFileAllSpots.getAbsolutePath())

    # Export XML
    writer = TmXmlWriter(xmlFile)
    writer.appendModel(tm.getModel())
    writer.appendSettings(tm.getSettings())
    writer.writeToFile()

    # --- Memory cleanup ---
    imp.changes = False
    imp.close()

    displayer = None
    trackTableView = None
    spotsTableView = None
    trackmate = None
    model = None
    settings = None
    imp = None

    System.gc()
    Runtime.getRuntime().gc()

    print("Export complete for", f.getName())

print("All positions processed.")
