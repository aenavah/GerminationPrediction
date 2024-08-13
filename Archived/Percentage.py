import pandas as pd
import matplotlib.pyplot as plt 


Directory = "/Users/alexandranava/Desktop/Spores/M4581_s1/PhC Analysis V2/V8/"
Summ_File = "Summary of Threshold_V8.csv"
path = Directory + Summ_File



if __name__ == "__main__":
  df = pd.read_csv(path)
  print("hit")
  df["Slice"] = [x.split("_")[-1] for x in df["Slice"]]
  print(df["Slice"])
  plt.xticks(ticks = [12, 36, 60, 84, 108, 132, 156, 180, 204, 228, 252, 276, 300], labels = [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25])
  
  plt.xlabel("Hours")
  plt.ylabel("Ungerminated Spore Count")
  
  plt.grid()
  plt.plot(df["Slice"], df["Count"])
  plt.savefig(Directory + "Plot.jpg")