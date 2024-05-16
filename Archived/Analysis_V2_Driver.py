from Analysis_Functions import Main as Main

# USER INPUTS:

Analysis_base: str = "/Users/alexandranava/Desktop/Spores/M4581_s1/Analysis/V2/"

PhC_base: str = "/Users/alexandranava/Desktop/Spores/M4581_s1/PhC Analysis V2/V14/"
PhC_Trackfile: str = "TrackTables_Tracks.csv"
PhC_Tracks_path: str = PhC_base + PhC_Trackfile

ThT_base = "/Users/alexandranava/Desktop/Spores/M4581_s1/ThT Analysis V2/V2/"
ThT_Trackfile: str = "TrackTables_Spots.csv"
ThT_Spots_path = ThT_base + ThT_Trackfile

Time_Between_Frames: int = 5 #minutes between frames
number_of_frames = 289

ThT_folder_by_track: str = "ThT_by_Spore/"
PostProcessed_Folder: str = "ThT_by_Spore_PostProcessed/"

output_matched_folder = "Spore_Data/"
plot_folder = "Data_Plots/"
matched_tracks_csv = "MatchedTracks.csv"

model_csv = "Model_Data.csv"
Main(PhC_Tracks_path, ThT_Spots_path, Analysis_base, Time_Between_Frames, ThT_folder_by_track, PostProcessed_Folder, output_matched_folder, number_of_frames, plot_folder, matched_tracks_csv, model_csv)