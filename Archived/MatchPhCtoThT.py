import numpy as np 
import pandas as pd
from pathlib import Path
import os

def match_files(ThT_path, PhC_path, base,  x_tol = 20, y_tol = 20):
  '''returns df and creates csv of matched track numbers based on position'''

  ThT_csv = pd.read_csv(ThT_path)
  PhC_csv = pd.read_csv(PhC_path)



  PhC_x = PhC_csv["TRACK_X_LOCATION"]
  PhC_y = PhC_csv["TRACK_Y_LOCATION"]

  #sort dfs by x then y for quicker analysis
  ThT_sorted = ThT_csv.sort_values(by = ["TRACK_X_LOCATION", "TRACK_Y_LOCATION"])
  PhC_sorted = PhC_csv.sort_values(by = ["TRACK_X_LOCATION", "TRACK_Y_LOCATION"])
  ThT_sorted.to_csv(base + "ThT_sorted.csv")
  PhC_sorted.to_csv(base + "PhC_sorted.csv")


  track_data = []
  #looping through ThT rows:
  for ThT_row_index in range(1, len(ThT_sorted)-2):
    #getting x,y coords for tht
    ThT_x = ThT_sorted.iloc[ThT_row_index, ThT_sorted.columns.get_loc("TRACK_X_LOCATION")]
    ThT_y = ThT_sorted.iloc[ThT_row_index, ThT_sorted.columns.get_loc("TRACK_Y_LOCATION")]

    #skip bad ThT data:
    if float(ThT_sorted.iloc[ThT_row_index, ThT_sorted.columns.get_loc("TRACK_START")]) == 0.0:
      continue
    if float(ThT_sorted.iloc[ThT_row_index, ThT_sorted.columns.get_loc("TRACK_DURATION")]) <= 280.0:
      continue


    #looping through PhC rows:
    for PhC_row_index in range(1, len(PhC_sorted)):
      PhC_x = PhC_sorted.iloc[PhC_row_index, PhC_sorted.columns.get_loc("TRACK_X_LOCATION")]
      PhC_y = PhC_sorted.iloc[PhC_row_index, PhC_sorted.columns.get_loc("TRACK_Y_LOCATION")]


      #if they are within x tolerance of eachother and y tolerance: 
      x_diff = abs(float(ThT_x) - float(PhC_x))
      y_diff = abs(float(ThT_y) - float(PhC_y))
      if ((x_diff < x_tol) and (y_diff < y_tol)):
        ThT_TRACK_ID = ThT_sorted.iloc[ThT_row_index, ThT_sorted.columns.get_loc("TRACK_INDEX")]
        PhC_TRACK_ID = PhC_sorted.iloc[PhC_row_index, PhC_sorted.columns.get_loc("TRACK_INDEX")]
        GERMINATION_TIME = PhC_sorted.iloc[PhC_row_index, PhC_sorted.columns.get_loc("TRACK_STOP")] + 1

        #print("Paired: -----------------")
        #print("ThT TRACK_ID:")
        #print(ThT_TRACK_ID)
        #print("(x, y): (" + str(ThT_x) + ", " + str(ThT_y) + ")")
        #print("PhC TRACK_ID:")
        #print(PhC_TRACK_ID)
        #print("(x, y): (" + str(PhC_x) + ", " + str(PhC_y) + ")")

        #get PhC "TRACK_ID" ThT "TRACK_ID", and X, Y coords, and germinaiton time
        track_data.append([ThT_TRACK_ID, PhC_TRACK_ID, ThT_x, ThT_y, PhC_x, PhC_y, GERMINATION_TIME])
      
  track_df = pd.DataFrame(track_data, columns = ["ThT TRACK ID", "PhC TRACK ID", "ThT X POSITION","ThT Y POSITION", "PhC X POSITION", "PhC Y POSITION", "PhC GERMINATION"])
  print("Total matched:")
  print(len(track_data))
  track_df.to_csv(base + "Matched.csv", index = False)
  return track_df
  
def concat_ThT_PhC_data(ThT_folder, Matched_Tracks_df, Output_base, Data_Folder = "Data_by_Track/", Naming_Convention = "Track_", Data_Output_Folder = "Data_by_Spore/"):
  for row_index in range(len(Matched_Tracks_df)):
    ThT_track_index = Matched_Tracks_df.iloc[row_index, Matched_Tracks_df.columns.get_loc("ThT TRACK ID")]
    ThT_Track_file = ThT_folder + Data_Folder + Naming_Convention + ThT_track_index + ".csv"
    ThT_Spore_Data = pd.read_csv(ThT_Track_file)

    #extra data
    PhC_track_index = Matched_Tracks_df.iloc[row_index, Matched_Tracks_df.columns.get_loc("PhC TRACK ID")]
    Germination_Frame = Matched_Tracks_df.iloc[row_index, Matched_Tracks_df.columns.get_loc("PhC GERMINATION")]
    Time = ThT_Spore_Data["POSITION_T"]
    Intensity = ThT_Spore_Data["MEAN_INTENSITY_CH1"]
    PhC_x = Matched_Tracks_df.iloc[row_index, Matched_Tracks_df.columns.get_loc("PhC X POSITION")]
    PhC_y = Matched_Tracks_df.iloc[row_index, Matched_Tracks_df.columns.get_loc("PhC Y POSITION")]
    ThT_x = Matched_Tracks_df.iloc[row_index, Matched_Tracks_df.columns.get_loc("ThT X POSITION")]
    ThT_y = Matched_Tracks_df.iloc[row_index, Matched_Tracks_df.columns.get_loc("ThT Y POSITION")]


    Tracks = [Germination_Frame, ThT_track_index, PhC_track_index, PhC_x, PhC_y, ThT_x, ThT_y]
    print(Tracks)
    extra_data = pd.DataFrame([Tracks], columns = ["GERMINATION", "ThT TRACK ID", "PhC TRACK ID", "PhC X", "PhC Y", "ThT X", "ThT Y"])
    data_tmp = pd.concat([Time, Intensity, extra_data], axis = 1)
    data_tmp.to_csv(Output_base + Data_Output_Folder + "ThTTrack=" + str(ThT_track_index) + ".csv")





    