import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt

#============================PhC Functions==========================
def Read_tracks(csv_path, output_path) -> pd.DataFrame:
  '''
  For use with PhC data given track data from fiji
  takes in the trackscheme csv path 
  outputs a df with relevant data to be processed
  '''
  df: pd.DataFrame = pd.read_csv(csv_path, skiprows = [1, 2, 3], index_col = False)
  germination_df: pd.DataFrame = df[["LABEL", "TRACK_START", "TRACK_STOP", "TRACK_DURATION", "TRACK_X_LOCATION", "TRACK_Y_LOCATION"]]

  return germination_df

#post-processing functions:
def PostProcessing_TrackStart(df) -> pd.DataFrame:
  '''
  For use with PhC data given track data from fiji 
  drops rows whose track dont start with 0
  '''
  rows_to_drop_indices: list[int] = []

  trackstart_column: int = df.columns.get_loc("TRACK_START")
  for track_index in range(len(df)):
    track_start: int = df.iloc[track_index, trackstart_column]
    
    if track_start != 0:
      rows_to_drop_indices.append(track_index)
  
  df_edited = df.drop(rows_to_drop_indices)

  return df_edited

def PostProcessing_TrackLength(df) -> pd.DataFrame:
  '''
  For use with PhC data given track data from fiji
  drops rows whose tracks are less than {min_num_frames} frames
  '''

  min_num_frames = 12

  rows_to_drop_indices: list[int] = []

  df_out = []
  trackstart_column: int = df.columns.get_loc("TRACK_START")
  trackstop_column: int = df.columns.get_loc("TRACK_STOP")
  total_tracks = len(df)

  for track_index in range(total_tracks):
    track_start: int = df.iloc[track_index, trackstart_column]
    track_stop: int = df.iloc[track_index, trackstop_column]
    track_duration: int = (track_stop + 1) - track_start

    if track_duration >= min_num_frames:
      df_out.append(df.iloc[track_index])
    else:
      rows_to_drop_indices.append(track_index)


  print(f"Dropped {len(rows_to_drop_indices)} tracks with track duration less than {min_num_frames} frames...")
  df_out = pd.DataFrame(df_out, columns = df.columns)
  return df_out



def Germination_Time(df, time_per_frame) -> pd.DataFrame:
  '''
  for use with PhC data
  Adds a column for frame germinated "FRAME_GERMINATED" 
    and time germinated "TIME_GERMINATED" in minutes
  '''

  germination_frame_data : list[int] = []
  germination_time_data: list[int] = []
  
  trackstop_column: int = df.columns.get_loc("TRACK_STOP")
  for track_index in range(len(df)):
    germination_frame: int = df.iloc[track_index, trackstop_column] + 1
    germination_frame_data.append(germination_frame)
    
    germination_time: int = germination_frame * time_per_frame
    germination_time_data.append(germination_time)

  df["FRAME_GERMINATED"] = germination_frame_data
  df["TIME_GERMINATED"] = germination_time_data
  return df
#============================ThT Functions==========================
#Read in ThT data 
def Read_spots(csv_path, Analysis_base, individual_data_folder) -> list[str]:
  '''
  takes in ThT csv path
  creates groups by track, sorted by time 
  outputs csvs of each group item
  return list of paths of csvs
  '''

  individual_paths: list[str] = []

  df = pd.read_csv(csv_path, skiprows = [1,2,3], index_col = False)
  df = df[["TRACK_ID", "POSITION_X", "POSITION_Y", "FRAME", "MEAN_INTENSITY_CH1"]]
  groups_by_track = df.sort_values(["FRAME"]).groupby("TRACK_ID")

  for track_number, track_data in groups_by_track :
    save_path = Analysis_base + individual_data_folder + "Track_" + str(track_number)+ "_Unprocessed.csv"
    individual_paths.append(save_path)
    track_data.to_csv(save_path, index = False)
  return individual_paths
  
def PostProcessing_SpotStart(unprocessed_paths, analysis_base, individual_data_folder) -> list[str]:
  '''
  takes in list of paths of unprocessed csvs for each track 
  creates csvs where the frames starts at time zero
  returns list of paths of csvs
  '''
  df_processed_paths: list[str] = []
  for path in unprocessed_paths:
    df = pd.read_csv(path)
    if df["FRAME"][0] == 0:
      track_number = df["TRACK_ID"][0]
      output_path = analysis_base + individual_data_folder + "Track_" + str(track_number) + "_ProcessedSpotStart.csv"
      df.to_csv(output_path, index = False)
      df_processed_paths.append(output_path)
  return df_processed_paths

def PostProcessing_SpotLength(unprocessed_paths, analysis_base, postprocessed_folder, num_frames = 289) -> list[str]:
  '''
  takes in list of paths from with track data and determines whether there's less than 12 frames per track
  creates a csv for each track 
  returns list of paths of csvs
  '''
  
  df_processed_paths: list[str] = []

  for path in unprocessed_paths:
    df = pd.read_csv(path, index_col = False)
    frame_column_index = df.columns.get_loc("FRAME")
    track_start = df.iloc[0, frame_column_index]
    track_end = df.iloc[-1, frame_column_index]
    track_duration = track_end - track_start + 1

    if track_duration == num_frames:
      track_number = df["TRACK_ID"][0]
      output_path = analysis_base + postprocessed_folder + "Track_" + str(track_number) + "_PostProcessed.csv"
      df_processed_paths.append(output_path)
      df.to_csv(output_path, index = False)

  return df_processed_paths

      
#============================Matching Functions==========================
def ThT_mean_xypos(ThT_postprocessed_paths) -> pd.DataFrame :
  '''
  Takes in a list of ThT paths corresponding to each ThT track after postprocessing
  Returns a dataframe : [TRACK_ID, MEAN_X_POS, MEAN_Y_POS]
  '''

  ThT_track_positions: list[int, int, int] = []

  #iterating through ThT data after postprocessing
  for track_path in ThT_postprocessed_paths:

    track_df: pd.DataFrame = pd.read_csv(track_path)
    track_number: int = track_df["TRACK_ID"][0]

    mean_xpos: int = track_df["POSITION_X"].mean()
    mean_ypos: int = track_df["POSITION_Y"].mean()
    ThT_track_positions.append([track_number, mean_xpos, mean_ypos])
  id_df = pd.DataFrame(ThT_track_positions, columns = ["TRACK_ID", "MEAN_X_POSITION", "MEAN_Y_POSITION"])
  return id_df

def Match_ThTtoPhC_byPosition(ThT_positions, Analysis_base, PhC_data_path, matched_tracks_csv, x_tol = 6, y_tol = 6) -> pd.DataFrame:
  '''
  takes in PhC track file and ThT position df 
  return df: ["PhC TRACK", "ThT TRACK"]
  '''
  
  PhC_df = pd.read_csv(Analysis_base + PhC_data_path)

  PhC_track_ID_column_index: int = PhC_df.columns.get_loc("LABEL")
  PhC_xlocation_column_index: int = PhC_df.columns.get_loc("TRACK_X_LOCATION")
  PhC_ylocation_column_index: int = PhC_df.columns.get_loc("TRACK_Y_LOCATION")

  ThT_track_ID_column_index: int = ThT_positions.columns.get_loc("TRACK_ID")
  ThT_xlocation_column_index: int = ThT_positions.columns.get_loc("MEAN_X_POSITION")
  ThT_ylocation_column_index: int = ThT_positions.columns.get_loc("MEAN_Y_POSITION")
  
  num_matched: int = 0 
  tracks_df: list[str, str] = []
  #iterate through ThT data:
  for ThT_row_index in range(len(ThT_positions)):
    ThT_TrackID_number: int = ThT_positions.iloc[ThT_row_index, ThT_track_ID_column_index]
    ThT_TrackID: str = "Track_" + str(ThT_TrackID_number)
    ThT_xpos: int = ThT_positions.iloc[ThT_row_index, ThT_xlocation_column_index]
    ThT_ypos: int  = ThT_positions.iloc[ThT_row_index, ThT_ylocation_column_index]
    
    for PhC_row_index in range(len(PhC_df)):
      PhC_TrackID: str = PhC_df.iloc[PhC_row_index, PhC_track_ID_column_index]
      PhC_xpos: int = PhC_df.iloc[PhC_row_index, PhC_xlocation_column_index]
      PhC_ypos:int = PhC_df.iloc[PhC_row_index, PhC_ylocation_column_index]

      x_diff = abs(ThT_xpos - PhC_xpos)
      y_diff = abs(ThT_ypos - PhC_ypos)
      if ((x_diff < x_tol) and (y_diff < y_tol)):
        #print(f"ThT {ThT_TrackID} corresponds to PhC {PhC_TrackID}")
        tracks_df.append([ThT_TrackID, PhC_TrackID])
        num_matched += 1
  print(f"Total matched: {num_matched}")  
  matched_tracks_df = pd.DataFrame(tracks_df, columns = ["ThT_TRACK_ID", "PhC_TRACK_ID"])
  matched_tracks_df.to_csv(Analysis_base + matched_tracks_csv)

  return matched_tracks_df

def Compile_Data(base, ThT_data_folder_path, PhC_data_csv_path, output_data_folder_path, df, num_frames, plot_folder) -> list[str]:
  PhC_df = pd.read_csv(base + PhC_data_csv_path, index_col=False)
  
  ThT_track_column_index: int = df.columns.get_loc("ThT_TRACK_ID")
  PhC_track_column_index: int = df.columns.get_loc("PhC_TRACK_ID")
  num_spores: int = len(df)

  data_paths: list[str] = []
  #iterating through matched_df
  for spore_index in range(num_spores):
    ThT_track: str  = df.iloc[spore_index, ThT_track_column_index]
    PhC_track: str = df.iloc[spore_index, PhC_track_column_index]
    
    ThT_track_csv_path = base + ThT_data_folder_path + ThT_track + "_PostProcessed.csv"
    ThT_df = pd.read_csv(ThT_track_csv_path, index_col = False)

    #find track in PhC data
    PhC_data_row_index: int = PhC_df[PhC_df['LABEL']== PhC_track].index.values

    PhC_framegerminated_column_index: int = PhC_df.columns.get_loc("FRAME_GERMINATED")
    PhC_timegerminated_column_index: int = PhC_df.columns.get_loc("TIME_GERMINATED")

    frame_germinated: int = PhC_df.iloc[PhC_data_row_index[0], PhC_framegerminated_column_index]
    time_germinated: int = PhC_df.iloc[PhC_data_row_index[0], PhC_timegerminated_column_index]

    #concatenating PhC data to ThT data for each identified spore 
    PhC_extracted_data: list[int, int] = [[frame_germinated, time_germinated]]
    
    print(f"Matching ThT {ThT_track} to PhC {PhC_track}")
    PhC_extracted_df = pd.DataFrame(PhC_extracted_data, columns = ["FRAME_GERMINATED", "TIME_GERMINATED"])
    output_data_df = pd.concat([ThT_df, PhC_extracted_df], axis = 1)

    csv_path = base + output_data_folder_path + "Spore_" + str(ThT_track) + ".csv"
    data_paths.append(csv_path)
    output_data_df.to_csv(csv_path)
    
    #plotting for verification
    plt.plot(output_data_df["FRAME"], output_data_df["MEAN_INTENSITY_CH1"], label = "Fluorescence")
    
    germinant_given: list[str] = [12, 36, 60, 84, 108, 132, 156, 180, 204, 228, 252, 276]
    
    for germinant_time in germinant_given:
      plt.axvline(germinant_time, color = "grey", linestyle = "--")
    
    plt.axvline(output_data_df["FRAME_GERMINATED"][0], color='red', linestyle='--', label = "Germination")
    plt.title(f"Spore {ThT_track}")

    plt.xlabel("Frame")
    plt.ylabel("Intensity")

    plt.ylim(10_000, 30_000)

    plt.legend()
    fig_path: str = base + plot_folder + "ThT_" + str(ThT_track) + ".jpg"
    plt.savefig(fig_path)
    plt.clf()

  return data_paths
  
def convert_to_model_data(base, csv_output, data_paths):

  model_data: list[str] = []
  germinant_given: list[str] = [12, 36, 60, 84, 108, 132, 156, 180, 204, 228, 252, 276]
  
  bad_tracks: list[str] = []

  for spore_data in data_paths:
    for track in bad_tracks:
      if track in spore_data:
        break
    else:   
      df = pd.read_csv(spore_data)
      
      intensities_list: list[int] = df["MEAN_INTENSITY_CH1"].to_list()
      germinant_exposure_list: list[int] = [1 if i in germinant_given else 0 for i in range(1, len(intensities_list) + 1)]
      germination_frame: int = int(df["FRAME_GERMINATED"][0])
      germination_list: list[int] = [1 if i >= germination_frame else 0 for i in range(1, len(intensities_list) + 1)] 

      data_row = [str(intensities_list), str(germinant_exposure_list), str(germination_list)]
      model_data.append(data_row)

  model_df = pd.DataFrame(model_data, columns = ["Intensity - ThT", "Germinant Exposure", "Germinated - PhC"])
  model_df.to_csv(base + csv_output)

#============================Main Functions==========================
def Main(PhC_csv_path, ThT_csv_path, Analysis_base, Time_Between_Frames, ThT_folder_by_spore, ThT_PostProcessed_Folder, Output_Matched_Folder, number_of_frames, plot_folder, matched_tracks_csv, model_csv):
  #=================PhC Processing===================
  #Read in PhC
  PhC_data_in: pd.DataFrame = Read_tracks(PhC_csv_path, Analysis_base)
  
  #Post Processing PhC data:
  PhC_data_postproc1 = PostProcessing_TrackStart(PhC_data_in)
  PhC_data_postproc2 = PostProcessing_TrackLength(PhC_data_postproc1)

  PhC_data_postproc2.to_csv(Analysis_base + "Analysis_PhC_PostProcessed.csv", index = False)

  #Getting Germination Time from PhC and creating csv:
  PhC_data_out = Germination_Time(PhC_data_postproc2, Time_Between_Frames)
  PhC_data_output_path = "Analysis_PhC_Data.csv"
  PhC_data_out.to_csv(Analysis_base + PhC_data_output_path, index = False)
  
  #=================ThT Processing===================
  #Read in ThT
  ThT_data_individual_paths: pd.DataFrame = Read_spots(ThT_csv_path, Analysis_base, ThT_folder_by_spore)
  
  #PostProcessing: iterate through dfs for each track in funcs and create csv
  ThT_data_postproc1_paths: list[str] = PostProcessing_SpotStart(ThT_data_individual_paths, Analysis_base, ThT_folder_by_spore)
  ThT_data_postproc2_paths: list[str] = PostProcessing_SpotLength(ThT_data_postproc1_paths, Analysis_base, ThT_PostProcessed_Folder)

  #=================Final Processing===================
  ThT_meanposition_df = ThT_mean_xypos(ThT_data_postproc2_paths)
  Matched_Tracks_df = Match_ThTtoPhC_byPosition(ThT_meanposition_df, Analysis_base, PhC_data_output_path, matched_tracks_csv)
  Final_data_paths: list[str] = Compile_Data(Analysis_base, ThT_PostProcessed_Folder, PhC_data_output_path, Output_Matched_Folder, Matched_Tracks_df, number_of_frames, plot_folder)
  convert_to_model_data(Analysis_base, model_csv, Final_data_paths)
