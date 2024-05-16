import pandas as pd
import matplotlib.pyplot as plt

def PhC_Germination_Percentage(base, df, num_frames):
  '''Takes in a df from a track scheme
  outputs "Analysis_Ignored_Tracks" - a csv with tracks that dont start at 0
  outputs "Analysis_Germinated_by_Frame.csv" - tracks final index of nongerminated spores
  '''

  #making df to append ignored tracks to 
  df_ignored = pd.DataFrame(columns = ["TRACK_INDEX", "TRACK_START", "TRACK_STOP", "TRACK_DURATION", "TRACK_X_LOCATION", "TRACK_Y_LOCATION"])
  #new df to append used tracks to 
  df_tracked = pd.DataFrame(columns = ["TRACK_INDEX", "TRACK_START", "TRACK_STOP", "TRACK_DURATION", "TRACK_X_LOCATION", "TRACK_Y_LOCATION"])

  #making subset df of input df with desired columns to evaluate 
  df_track_frames = df[["TRACK_INDEX", "TRACK_START", "TRACK_STOP", "TRACK_DURATION", "TRACK_X_LOCATION", "TRACK_Y_LOCATION"]]
  
  # used to track how many become germinated at each frame 
  germinated_by_track = {} #key = frame number germination occurs, value = number of spores
  

  for track_index in range(len(df_track_frames)):
    ignore = 0 #goes through loop and determines whether it's an acceptable track to be added to df_tracked

    track_start = df_track_frames.iloc[track_index, 1]
    track_stop = df_track_frames.iloc[track_index, 2]

    #ignore tracks that don't start at 0
    if track_start != 0:
      df_ignored_index = len(df_ignored)
      df_ignored.loc[df_ignored_index] = df.iloc[track_index]
      #df_track_frames.iloc[track_index, :]
      continue

    #only include tracks that have 5 frames
    if track_stop - track_start < 5:
      df_ignored_index = len(df_ignored)
      df_ignored.loc[df_ignored_index] = df.iloc[track_index]
      continue

    #keeping track of tracks that are acceptable
    df_tracked_index = len(df_tracked)
    df_tracked.loc[df_tracked_index] = df.iloc[track_index]

    #store how many times a frame is last non-germinated point  
    if track_stop in germinated_by_track:
      germinated_by_track[track_stop] += 1
    else:
      germinated_by_track[track_stop] = 1

  #gets mean number of spores per frames that are germinated!!!
  count_list = [0] * (num_frames)

  #store how many spores are germinated by frame
  for key in germinated_by_track.keys(): 
    for frame in list(range(num_frames)):
      if frame <= key:
        count_list[frame] += int(germinated_by_track[key] + 1)
        
      
  #create csv of ignored trackss
  
  df_germinated_by_frame = pd.DataFrame({'FRAME': list(range(num_frames)),
                                           'SPORE COUNT': count_list})  
  
  #output ignored tracks, used tracks, and spores germinated by frame
  df_ignored.to_csv(base + "Analysis_Ignored_Tracks.csv")
  df_tracked.to_csv(base + "Analysis_Used_Tracks.csv")
  df_germinated_by_frame.to_csv(base + "Analysis_Germinated_by_Frame.csv")

def Plot(x, y, title, filename, xlabel, ylabel):
  plt.plot(x, y)

  vertical_values = [11, 35, 59, 83, 107, 131, 155, 179, 203, 227, 251, 275]
  for value in vertical_values:
    plt.axvline(x=value, color='grey', linestyle='dotted')

  plt.title(title)
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)
  plt.savefig(filename)

def Intensities(df, base):
  '''
  takes in df of all values from track table of spots, 
  produces plots only for tracks that start at frame 0 in "\Plot_by_Track\"
  saves csv of all separated tracks in "\Data_by_Track\"
  '''

  groups_by_track = df.sort_values(["FRAME"]).groupby("TRACK_ID")

  tracked_tracks = []
  dropped_tracks = []
  
  # group dataframe into tracks and drops tracks that dont start at 0
  for track_number, track_data in groups_by_track :
    intensity_column = track_data["MEAN_INTENSITY_CH1"]
    frame_column = track_data["FRAME"]
    if frame_column[0] == 0:
      tracked_tracks.append(track_number)
      plt.clf()
      plt.ylim(0, 40000)
      Plot(frame_column, intensity_column, "Spore " + str(track_number), base + "Plot_by_Track/Track_" + str(track_number)+ ".jpg", "Frame", "Intensity")  
    if frame_column[0] != 0:
      dropped_tracks.append(track_number)
    track_data.to_csv(base + "Data_by_Track/Track_" + str(track_number)+ ".csv")
      
  mean_by_frame = {}  
  # groups tracks by frame for finding mean intensity
  groups_by_frame = df.sort_values(["TRACK_ID"]).groupby("FRAME")
  for frame_number, group_data in groups_by_frame:

    # remove tracks that dont start at frame 0 
    for track in dropped_tracks:
      group_data = group_data[group_data["TRACK_ID"] != track]
    
    # find mean
    mean_intensity = group_data["MEAN_INTENSITY_CH1"].mean()
    mean_by_frame[frame_number] = mean_intensity
    df_mean_intensity = pd.DataFrame(list(mean_by_frame.items()), columns=['FRAME', 'MEAN_INTENSITY'])

  #plot
  plt.clf()
  x = df_mean_intensity["FRAME"]
  y = df_mean_intensity["MEAN_INTENSITY"]
  plt.ylim(0, 40000)
  Plot(x, y, "Average Intensity of Spores", base + "MeanIntensity.jpg", "Frames", "Intensity")

