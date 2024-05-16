import shutil 
import os

def copy_tifs(source, destination, file):
  total = sum(1 for item in os.listdir(source))

  for index in range(total):
    if index % 12 == 0.0:
      tmp_number = str(index).zfill(4)
      tif = file + tmp_number + ".tif"
      shutil.copy(source + tif, destination)


if __name__ == "__main__":
  source = "/Users/alexandranava/Desktop/Spores/Task - Predicting Germination/ThT/"
  destination = "/Users/alexandranava/Desktop/Spores/Task - Predicting Germination/Analysis/Frames by Hour/"
  file = "M4581_s1_ThT_stabilized_"
  copy_tifs(source, destination, file)