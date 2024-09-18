import matplotlib.pyplot as plt
import pandas as pd 

def plot_single(base, csv_indv, frametime, spore_nums, datatype):

  '''
  base is the folder where data is saved 
  csv_indv is naming convention of csv for each spore 
  frametime is the minutes between each frame
  spore_nums is the spore ids to be plotted [2,6] plots spores with id 2,3,4,5,6
  '''

  spore_ids = list(range(spore_nums[0], spore_nums[-1] + 1))
  data = pd.DataFrame()

  for spore in spore_ids: 
    spore_csv = csv_indv + str(spore) + ".csv"
    print("Working on " + spore_csv + "...")
    spore_df = pd.read_csv(base + spore_csv)
    spore_mean_df = spore_df["Mean"]

    data = pd.concat([data, spore_mean_df], axis = 1)

  #print(data) now data has columns corresponding to means from each spore 
  num_frames = len(data)
  frame_indices = list(range(num_frames)) # list of numbers from 0 to number of frames 
  plt.xlabel("Hour", fontsize = 18)
  plt.ylabel("Electrochemical Potential", fontsize = 18)

  frame_ind_map = []
  frame_hr_map = []

  #used for converting frame number to hours on xaxis for plotting
  for frame_ind in frame_indices:
    whole_hour = 0 
    frame_min = frame_ind * frametime #frame converted to minute
    frame_hour = frame_min / 60
    if frame_hour % 1 == 0.0: #if frame hour is a whole number: 
      whole_hour = 1   
    
    if whole_hour == 1:
      print("mapping " + str(frame_ind) + " to " + str(frame_hour) + " ...")
      frame_ind_map.append(frame_ind)
      frame_hr_map.append(int(frame_hour))
  
  print(data)
  for spore in spore_ids:
    column_data = data.iloc[:, spore - spore_nums[0] - 1]  # Selecting column by index
    plt.plot(frame_indices, column_data, label = "Spore " + str(spore))

  for frame_hour in frame_hr_map:
    plt.axvline(x=frame_hour * 60 / frametime, linestyle='dotted', color='gray') 

  plt.title("Analysis of "+ datatype + " Data", fontsize = 20)
  plt.legend()
  plt.xticks(ticks = frame_ind_map, labels = frame_hr_map)
  plt.savefig(base + "IndividualSpores_" + str(spore_nums[0])+"-"+ str(spore_nums[-1]))


def plot_mean(df, csv, frametime, datatype):
  '''
  plots mean of all spores at each frame

  df holds all mean values 
  csv is name of the plot to be saved
  frametime is interval at which photo is taken 
  '''
  means = df.mean(axis = 1).tolist() #  mean per frame frame 
  frame_indices = list(range(len(means))) # list of numbers from 0 to number of frames 
  plt.xlabel("Hour", fontsize = 18)
  plt.ylabel("Percentage of Dormant Spores", fontsize = 18)

  frame_ind_map = []
  frame_hr_map = []

  #used for converting frame number to hours on xaxis for plotting
  for frame_ind in frame_indices:
    whole_hour = 0 
    frame_min = frame_ind * frametime #frame converted to minute
    frame_hour = frame_min / 60
    if frame_hour % 1 == 0.0: #if frame hour is a whole number: 
      whole_hour = 1   
    
    if whole_hour == 1:
      print("mapping " + str(frame_ind) + " to " + str(frame_hour) + " ...")
      frame_ind_map.append(frame_ind)
      frame_hr_map.append(int(frame_hour))
  
  plt.title("Analysis of " + datatype + " Data", fontsize = 20)
  plt.plot(frame_indices, means)

  for frame_hour in frame_hr_map:
    plt.axvline(x=frame_hour * 60 / frametime, linestyle='dotted', color='gray') 

  plt.xticks(ticks = frame_ind_map, labels = frame_hr_map)
  plt.yticks(ticks = [50, 127.5, 255], labels = [25, 50, 100])
  plt.savefig(base + "Mean_" + csv)




def plot_all(df, csv): 
  '''
  Plots individual spores at each frame 
  '''
  spore_lists = df.values.tolist() # each lists corrresponds to a spore
  num_spores = list(range(len(spore_lists[0]))) # list of numbers from 0 to amount of spores
  frames_in_min = []
  frame_index = list(range(len(spore_lists)))
  for frame_number in frame_index:
    print(frame_index)
    frame_min = frame_index * 5
    frame_index.append(frame_min)
  print(frame_index)
  for spore_index in num_spores:
    spore_vals = df.iloc[:, spore_index]
    plt.plot(frame_index, spore_vals)
  plt.xlabel("Frame", fontsize = 18)
  plt.ylabel("Intensity", fontsize = 18)
  plt.xticks(range(frame_index[-1] + 1))
  plt.title("Analysis", fontsize = 20)
  plt.savefig(base + "AllSpores_" + csv)

if __name__ == "__main__":
  '''
  data:
  each row is a frame
  each column is the intensity of a spore 
  '''
  #------------ needed for all functions
  
  min_between_frame = 5
  #typee = "individual" #all or individual 
  datatype = "PhC" #which oneeeee
  if datatype == "PhC":
    base = "/Users/alexandranava/Desktop/Spores/M4581_s1/PhC Analysis V2/"
  if datatype == "ThT":
    base = "/Users/alexandranava/Desktop/Spores/M4581_s1/ThT Analysis/Analysis/Single Spores/"

  #plotting mean grey values of all the spores 
  if datatype == "PhC": 
    csv_all = "Results_V8" ###################
    mean_path = base + csv_all + ".csv"

    print("Getting data from " + " " + mean_path + "...")
    df = pd.read_csv(mean_path)
    df = df.iloc[:,1:]  #  removing index column
    plot_mean(df, csv_all, min_between_frame, datatype) # plots mean of grey values at each frame


  #plotting single spores: 
  if datatype == "ThT":
    spore_nums = [1, 1] #usings spores with ids 2 - 6 [2, 6]
    csv_indv = "Spore" #name of csv file for each individual spore, this one corresponds to Spore2.csv
    
    print("Plotting individual spores from " + str(spore_nums[0]) + " to " + str(spore_nums[-1]) + "...")
    plot_single(base, csv_indv, min_between_frame, spore_nums, datatype) # plots means at each frame
