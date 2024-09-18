import matplotlib.pyplot as plt
import pandas as pd
global PhC_base



PhC_base = "/Users/alexandranava/Desktop/Spores/M4581_s1/PhC Analysis V2/V6/"
spore_id = 10

datatype = "PhC"
csv_name = "Spore" + str(spore_id) + "_Results.csv"
germ_threshold = 30 #if mean grey value above this number, consider germinated


def plot(base, csv, spore_id):
  print(csv)
  #--- Intensities
  df = pd.read_csv(base + csv)
  print(df.head)
  plt.plot(df["Slice"], df["Mean"], label = "Fiji Mean Grey Value", color = "grey")
  plt.plot(df["Slice"], df["Binary"], label = "Binary", color = "black")
  
  plt.xlabel("Frame")
  plt.xticks(ticks = [12, 36, 60, 84, 108, 132, 156, 180, 204, 228, 252, 276, 300], labels = [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25])
  plt.ylabel("Intensity")
  
  plt.legend()
  plt.grid()
  plt.title("Intensity Values: Spore " + str(spore_id))

  plt.savefig(PhC_base + "Comparison_" + csv.replace(".csv", ".jpg"))
  plt.clf()

  #--- Errors
  plt.plot(df["Slice"], df["Error"], label = "Error")

  plt.xticks(ticks = [12, 36, 60, 84, 108, 132, 156, 180, 204, 228, 252, 276, 300], labels = [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25])
  plt.xlabel("Frame")
  plt.ylabel("Intensity")
  
  plt.legend() 
  plt.grid() 
  plt.title("Difference between Actual and Binary: Spore " + str(spore_id))
  plt.ylim(0, 255)

  plt.savefig(PhC_base + "Error_" + csv.replace(".csv", ".jpg"))
  plt.clf()

def convert_binary(germ_threshold, df):

  '''
  this functions takes in a number "germ_threshold" where if the mean grey value is above that threshold, it is considered germinated
  takes in df 
  returns: df holding a column of the actual greys vals and estimated one, and the difference
  '''

  binary_values = []
  errors = []
  actual_greyvals_df = df["Mean"]
  for value in actual_greyvals_df:
    if value > germ_threshold:
      value_approx  = 255
    else:
      value_approx = 0
    binary_values.append(value_approx)
    errors.append(abs(value - value_approx))

  df.insert(3, "Binary", binary_values)
  df.insert(4, "Error", errors)
  return df


if __name__ == "__main__":
  spore_path = PhC_base + csv_name
  df = pd.read_csv(spore_path)
  data = convert_binary(germ_threshold, df)
  new_data_path = PhC_base + "Calculated_" + csv_name
  data.to_csv(new_data_path)

  plot(PhC_base, "Calculated_" + csv_name, spore_id)