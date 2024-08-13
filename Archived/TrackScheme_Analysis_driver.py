import TrackScheme_Analysis_functions
import MatchPhCtoThT
import pandas as pd

PhC_Germination_Percentage = TrackScheme_Analysis_functions.PhC_Germination_Percentage  
Plot = TrackScheme_Analysis_functions.Plot
Intensities = TrackScheme_Analysis_functions.Intensities
match_files = MatchPhCtoThT.match_files
ConcatData = MatchPhCtoThT.concat_ThT_PhC_data

if __name__ == "__main__":
  datatype = "PhC"
  #PhC Analysis ---------
  if datatype == "PhC":
    PhC_base = "/Users/alexandranava/Desktop/Spores/M4581_s1/PhC Analysis V2/V13/"
    PhC_trackscheme_file = "TrackTables_Tracks.csv"
    num_frames = 289 
  #---
    PhC_trackscheme_path = PhC_base + PhC_trackscheme_file
    df = pd.read_csv(PhC_trackscheme_path, index_col = 0, skiprows = [1,2,3])
    print("\n")
    print("Reading Track Scheme csv...\n")
    
    #creating data of last frame of non-germinated spores
    PhC_Germination_Percentage(PhC_base, df, num_frames)

    #plotting frames at which spores turn off
    germination_by_frame_csv = "Analysis_Germinated_by_Frame.csv"
    germination_by_frame_path = PhC_base + germination_by_frame_csv

    df_germ_by_frame = pd.read_csv(germination_by_frame_path)
    frames = df_germ_by_frame["FRAME"]
    spore_count = df_germ_by_frame["SPORE COUNT"]
    Plot(frames, spore_count, "Ungerminated Spore Count", PhC_base + "Germ_by_Frame.jpg", "Frames", "Number of Spores")

  #ThT Analysis ---------
  #creating dataframe of spore intensities
  if datatype == "ThT":
    ThT_base = "/Users/alexandranava/Desktop/Spores/M4581_s1/ThT Analysis V2/V1/"
    ThT_trackscheme_file = "TrackTables_Spots.csv" 

    #---
    ThT_trackscheme_path = ThT_base + ThT_trackscheme_file
    ThT_df = pd.read_csv(ThT_trackscheme_path, index_col = 0, skiprows = [1,2,3])
    Intensities(ThT_df, ThT_base)

#Matching spores from ThT and PhC data
PhC_postprocessed_csv = "/Users/alexandranava/Desktop/Spores/M4581_s1/PhC Analysis V2/V13/Analysis_Used_Tracks.csv"
ThT_base = "/Users/alexandranava/Desktop/Spores/M4581_s1/ThT Analysis V2/V1/"
ThT_csv = ThT_base + "TrackTables_Tracks.csv"

Analysis_base = "/Users/alexandranava/Desktop/Spores/M4581_s1/Analysis/"

matched_tracks = match_files(ThT_csv, PhC_postprocessed_csv, Analysis_base)

ConcatData(ThT_base, matched_tracks, Analysis_base)

