import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import os 

def create_directories(base, ThT_preprocessed_folder = "ThT_preprocessed/", PhC_preprocessed_folder = "PhC_preprocessed/", PhC_postprocessed_folder = "PhC_postprocessed/", ThT_postprocessed_folder = "ThT_postprocessed/", spore_data_folder = "Spore_Data/", plot_folder = "Data_Plots/"):
  '''
  creates folders for output if they dont exist and returns their paths'''
  folder_paths = []

  for folder_name in [ThT_preprocessed_folder, PhC_preprocessed_folder, ThT_postprocessed_folder, PhC_postprocessed_folder, spore_data_folder, plot_folder]:
    folder_path = base + folder_name
    folder_paths.append(folder_path)
    if not os.path.exists(folder_path):
      os.makedirs(folder_path)  
  
  return folder_paths

def read_spots(csv_path, base, output_folder) -> list[str]:
  '''
  takes in spot csv path
  creates groups by track, sorted by time 
  outputs csvs of each group item
  return list of full paths of csvs
  '''

  individual_paths: list[str] = []

  df = pd.read_csv(csv_path, skiprows = [1,2,3], index_col = False)
  groups_by_track = df.sort_values(["FRAME"]).groupby("TRACK_ID")

  for track_number, track_data in groups_by_track :
    save_path = output_folder + "TRACK_" + str(track_number)+ ".csv"
    individual_paths.append(save_path)
    track_data.to_csv(save_path, index = False)
  return individual_paths

def post_process(csv_paths, output_path, data_type, frame_number):
  '''
  removes tracks that dont exists at time 0 
  for ThT data, confirms the tracks exists for all frames
  for PhC data, removes them if they're less than an hour 
  returns csv data paths that have been considered acceptable
  '''
  
  print(f"filtering {len(csv_paths)} tracks...")

  processed_paths: list[str] = []

  if data_type == "PhC":
    min_num_frames = 12
  if data_type == "ThT":
    min_num_frames = frame_number

  for csv_path in csv_paths:
    df = pd.read_csv(csv_path)
    track_id = df["TRACK_ID"][0]
    num_frames_tracked = len(df)
    first_frame_index = df["FRAME"][0]
    if num_frames_tracked < min_num_frames: 
      continue
    if first_frame_index != 0:
      continue
    
    if data_type == "PhC": #if PhC data, adds a columns that denotes when the frame is germinated
      germination_df = pd.DataFrame([[num_frames_tracked]], columns = ["GERMINATION_FRAME"])
      df = pd.concat([df, germination_df], axis = 1)

    output_csv = output_path + "TRACK_" + str(track_id) + ".csv"
    processed_paths.append(output_csv)
    df.to_csv(output_csv, index=False)

  print(f"filtered out {len(csv_paths) - len(processed_paths)} tracks...")
  print(f"left with {len(processed_paths)} tracks...")
  
  return processed_paths

def match_positions(ThT_csv_paths, PhC_csv_paths, x_tol = 5, y_tol = 5) -> list[str, str]:
  '''
  reads in two lists correpsonding to ThT csv paths and PhC csv paths
  determines mean x and y positions for each 
  if they're the same within a tolerance, paths are added as a tuple into a list 
  returns list of matched paths from ThT and PhC
  '''
  paired_df_paths: list[str, str] = []
  locations: list[int, int] = []

  for ThT_csv_path in ThT_csv_paths:
    ThT_df = pd.read_csv(ThT_csv_path)
    ThT_xpos = ThT_df["POSITION_X"].mean()
    ThT_ypos = ThT_df["POSITION_Y"].mean()
    for PhC_csv_path in PhC_csv_paths:
      PhC_df = pd.read_csv(PhC_csv_path)
      PhC_xpos = PhC_df["POSITION_X"].mean()
      PhC_ypos = PhC_df["POSITION_Y"].mean()

      diff_x = abs(ThT_xpos - PhC_xpos)
      diff_y = abs(ThT_ypos - PhC_ypos)

      if ((diff_x < x_tol) and (diff_y < y_tol)):
        paired_df_paths.append([ThT_csv_path, PhC_csv_path])
        locations.append([PhC_xpos, PhC_ypos])
        break 

  print(f"matched {len(paired_df_paths)} dataframes...")
  return paired_df_paths, locations

def concatenate_dfs(output_folder, paired_df_paths) -> list[str]:
  '''
  reads in list of matched ThT and PhC paths
  concatenates ThT df to germination frame from PhC df
  returns list of paths of concatenated dfs
  '''
  data_paths = []

  for path_pair in paired_df_paths:
    ThT_path = path_pair[0]
    PhC_path = path_pair[1]
    ThT_df = pd.read_csv(ThT_path)
    PhC_df = pd.read_csv(PhC_path)
    germination_frame = PhC_df["GERMINATION_FRAME"]
    ThT_track = ThT_df["TRACK_ID"][0]

    df = pd.concat([ThT_df, germination_frame], axis = 1)
    df_path = output_folder + str(ThT_track) + ".csv"
    data_paths.append(df_path)
    df.to_csv(df_path)

  print(f"data retrieved...")
  return data_paths

def convert_to_modeldata(csv_paths, output_csv, germinant_given, plot_folder, spatial_df) -> None:
  '''
  takes in list of paths of csv data 
  converts to df readable to the model and creates a csv file with all data 
  '''
  model_data = []
  model_data_V2 = []
  spore_index = 0
  for spore_data in csv_paths:
      df = pd.read_csv(spore_data)
      spore_id_list: list[int] = df["TRACK_ID"].to_list()
      frame_list: list[int] = df["FRAME"].to_list()
      intensities_list: list[int] = df["MEAN_INTENSITY_CH1"].to_list()
      area_list = df["AREA"].to_list()
      germinant_exposure_list: list[int] = [1 if i in germinant_given else 0 for i in range(0, len(intensities_list))]
      germination_frame: int = int(df["GERMINATION_FRAME"][0])
      germination_list: list[int] = [1 if i >= germination_frame else 0 for i in range(0, len(intensities_list))] 
      ellipse_minor_list: list[int] = df["ELLIPSE_MINOR"].to_list()
      ellipse_major_list: list[int] = df["ELLIPSE_MAJOR"].to_list()
      perimeter_list: list[int] = df["PERIMETER"].to_list()
      circularity_list: list[int] = df["CIRCULARITY"].to_list()
      xpos: int = np.mean(df["POSITION_X"])
      ypos: int = np.mean(df["POSITION_Y"])

      for i in range(len(spore_id_list)):
        model_data_V2.append([spore_id_list[i], frame_list[i], xpos, ypos, germination_list[i], intensities_list[i], area_list[i], germinant_exposure_list[i], ellipse_minor_list[i], ellipse_major_list[i], perimeter_list[i], circularity_list[i]])

      data_row = [str(intensities_list), str(area_list), str(germinant_exposure_list), str(germination_list), str(ellipse_minor_list), str(ellipse_major_list), str(perimeter_list), str(circularity_list), str(frame_list)]
      model_data.append(data_row)

      if plot_data == 1:
        #-----plotting intensity
        plt.clf()
        plt.plot(df["FRAME"], df["MEAN_INTENSITY_CH1"])
        
        for germinant_time in germinant_given:
          plt.axvline(germinant_time, color = "lightgrey", linestyle = "--")    
        plt.axvline(germination_frame, color='red', linestyle='--', label = "Germination")
        
        plt.title(f"Spore {spore_index}")
        plt.xlabel("Frame")
        plt.ylabel("Intensity")
        plt.ylim(0, 100)
        
        plt.legend()
        plt.savefig(plot_folder + "Intensity_" + str(spore_index) + ".jpg")

        #-----plotting area
        plt.clf()
        plt.plot(df["FRAME"], df["AREA"])
        for germinant_time in germinant_given:
          plt.axvline(germinant_time, color = "lightgrey", linestyle = "--")    
        plt.axvline(germination_frame, color='red', linestyle='--', label = "Germination")
        
        plt.title(f"Spore {spore_index}")
        plt.xlabel("Frame")
        plt.ylabel("Area")
        plt.ylim(0, 550)
        
        plt.legend()
        plt.savefig(plot_folder + "Area_" + str(spore_index) + ".jpg")

        #-----plotting ellipses
        plt.clf()
        plt.plot(df["FRAME"], df["ELLIPSE_MINOR"], label = "Ellipse Minor")
        plt.plot(df["FRAME"], df["ELLIPSE_MAJOR"], label = "Ellipse Major")
        for germinant_time in germinant_given:
          plt.axvline(germinant_time, color = "lightgrey", linestyle = "--")    
        plt.axvline(germination_frame, color='red', linestyle='--', label = "Germination")
        
        plt.title(f"Spore {spore_index}")
        plt.xlabel("Frame")
        plt.ylabel("Ellipse")
        plt.ylim(5, 17)
      
        plt.legend()
        plt.savefig(plot_folder + "Ellipse_" + str(spore_index) + ".jpg")

        #-----plotting perimeters
        plt.clf()
        plt.plot(df["FRAME"], df["PERIMETER"])
        for germinant_time in germinant_given:
          plt.axvline(germinant_time, color = "lightgrey", linestyle = "--")    
        plt.axvline(germination_frame, color='red', linestyle='--', label = "Germination")
        
        plt.title(f"Spore {spore_index}")
        plt.xlabel("Frame")
        plt.ylabel("Perimeter")
        plt.ylim(50, 100)
    
        plt.legend()
        plt.savefig(plot_folder + "Perimeter_" + str(spore_index) + ".jpg")

      # iterate through next spore 
      spore_index += 1

  #compact verisno
  model_df = pd.DataFrame(model_data, columns = ["INTENSITY", "AREA", "GERMINANT EXPOSURE", "GERMINATION", "ELLIPSE MINOR", "ELLIPSE MAJOR", "PERIMETER", "CIRCULARITY", "FRAME LIST"])
  model_df = pd.concat([model_df, spatial_df], axis = 1)
  model_df.to_csv(output_csv)

  #long version
  model_df_V2 = pd.DataFrame(model_data_V2, columns = ["SPORE_ID", "FRAME", "X_POSITION", "Y_POSITION", "GERMINATION", "INTENSITY", "AREA", "GERMINANT_EXPOSURE", "ELLIPSE_MINOR", "ELLIPSE_MAJOR", "PERIMETER", "CIRCULARITY"])
  model_df_V2.to_csv(output_csv.replace(".csv", "_V2.csv"))
  print(f"converted to model data and plotted...")


def Main(Analysis_base, PhC_base, PhC_csv_name, ThT_base, ThT_csv_name, num_frames, spore_set_id, xtol, ytol):
  print("\n Running! \n")
  PhC_spot_path = PhC_base + PhC_csv_name
  ThT_spot_path = ThT_base + ThT_csv_name


  #make directories for data
  [ThT_preprocessed_folder, PhC_preprocessed_folder, ThT_postprocessed_folder, PhC_postprocessed_folder, spore_data_folder, plot_folder] = create_directories(Analysis_base)
  
  print(f"working on PhC data...")
  #divide entire spot file into csvs of tracks and process
  #PhC
  PhC_unprocessed_csv_paths: list[str] = read_spots(PhC_spot_path, Analysis_base, PhC_preprocessed_folder)
  PhC_processed_csv_paths: list[str] = post_process(PhC_unprocessed_csv_paths, PhC_postprocessed_folder, "PhC", num_frames)
  print("\n")
  #ThT
  print(f"working on ThT data...")
  ThT_unprocessed_csv_paths: list[str] = read_spots(ThT_spot_path, Analysis_base, ThT_preprocessed_folder)
  ThT_processed_csv_paths: list[str] = post_process(ThT_unprocessed_csv_paths, ThT_postprocessed_folder, "ThT", num_frames)

  #Matching tracks 
  print("\n")
  matched, locations = match_positions(ThT_processed_csv_paths, PhC_processed_csv_paths, xtol, ytol)
  data_paths = concatenate_dfs(spore_data_folder, matched)
  spatial_df = pd.DataFrame(locations, columns = ["X_POSITION", "Y_POSITION"])
  convert_to_modeldata(data_paths, Analysis_base + spore_set_id + "_Model_Data.csv", germinant_given, plot_folder, spatial_df)
  print("\n")
  print("Done!")
  
if __name__ == "__main__":
  plot_data = 1

  germinant_given: list[str] = [11, 35, 59, 83, 107, 131, 155, 179, 203, 227, 251, 275]

  '''M4581_s1 Analysis'''
  Analysis_base = "/Users/alexandranava/Desktop/Spores/M4581_s1/Analysis/V3/"

  PhC_base = "/Users/alexandranava/Desktop/Spores/M4581_s1/PhC Analysis V3/V3.1/"
  PhC_csv_name = "PhC_TrackTable_Spots.csv"


  ThT_base = "/Users/alexandranava/Desktop/Spores/M4581_s1/ThT Analysis V3/V3.1/"
  ThT_csv_name = "ThT_TrackTable_Spots.csv"
  num_frames = 289
  spore_set_id = "M4581_s1"
  Main(Analysis_base, PhC_base, PhC_csv_name, ThT_base, ThT_csv_name, num_frames, spore_set_id, 5, 5)

  #right now the plots are saved by range, can be mapped back to tht track in saving. 

  '''M4576_s2 Analysis'''
  Analysis_base = "/Users/alexandranava/Desktop/Spores/M4576_s2/"
  PhC_base = Analysis_base + "PhC Analysis/"
  ThT_base = Analysis_base + "ThT Analysis/"
  PhC_csv_name = "PhC_TrackTable_Spots.csv"
  ThT_csv_name = "ThT_TrackTable_Spots.csv"
  num_frames = 276
  spore_set_id = "M4576_s2"
  Main(Analysis_base, PhC_base, PhC_csv_name, ThT_base, ThT_csv_name, num_frames, spore_set_id, 3, 3)