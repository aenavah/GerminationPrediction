# run_analysis.py

from Data_Extraction import Main
from Data_Extraction import plot_spatial

# Define your parameters
plot_data = 1

output_base = "/Users/alexandranava/Desktop/Spores/Data/pH_pulses/M6543_s5/Fiji_Processing/PostProcessing/"
PhC_base = "/Users/alexandranava/Desktop/Spores/Data/pH_pulses/M6543_s5/Fiji_Processing/"
PhC_csv_name = "PhC_TrackTable_Spots.csv"
ThT_base = PhC_base
ThT_csv_name = "ThT_TrackTable_Spots.csv"
num_frames = 232
spore_set_id = "M6543_s5"
xtol = 4
ytol = 4

first_germinant = 1
time_between_germinant = 60
time_between_frames = 5

#plot_spatial(PhC_base + PhC_csv_name, 1800, 1024)
# Call the Main function
Main(output_base, PhC_base, PhC_csv_name, ThT_base, ThT_csv_name, num_frames, first_germinant, time_between_germinant, time_between_frames, spore_set_id, xtol, ytol, 1)