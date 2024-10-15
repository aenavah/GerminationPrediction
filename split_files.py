import os
import shutil

base = "/Users/alexandranava/Desktop/Spores/"
source_folder = "pH_pulses"

def split_files(base, source_folder):
  filenames = os.listdir(base + source_folder)
  for filename in filenames:

    if ".tif" not in filename:
      continue
    exp = filename.split("_")[0]
    set_num = filename.split("_")[1]
    imaging = filename.split("_")[2]

    if imaging == "PH":
      imaging = "PhC"

    destination_dir = base + exp + "_" + set_num 
    destination_imaging = destination_dir + "/" + imaging
    if not os.path.exists(destination_dir):
      os.makedirs(destination_dir)
    if not os.path.exists(destination_imaging):
      os.makedirs(destination_imaging)
    source_file = os.path.join(base, source_folder, filename)
    destination_file = os.path.join(destination_imaging, filename)
    shutil.move(source_file, destination_file)
  print(f"Moved {len(filenames)} files...")

split_files(base, source_folder)