# run_analysis.py

from Analysis_Functions import Main
from Analysis_Functions import plot_spatial

# Define your parameters
plot_data = 1

output_base = "/Users/alexandranava/Desktop/Spores/M4581_s7/PostProcess_V4/"
PhC_base = "/Users/alexandranava/Desktop/Spores/M4581_s7/Fiji_Data/"
PhC_csv_name = "PhC_TrackTable_Spots.csv"
ThT_base = PhC_base
ThT_csv_name = "ThT_TrackTable_Spots_V4.csv"
num_frames = 289
spore_set_id = "M4581_s7"
xtol = 5
ytol = 5

first_germinant = 60 
time_between_germinant = 120
time_between_frames = 5

#plot_spatial(PhC_base + PhC_csv_name, 1800, 1024)
# Call the Main function
Main(output_base, PhC_base, PhC_csv_name, ThT_base, ThT_csv_name, num_frames, first_germinant, time_between_germinant, time_between_frames, spore_set_id, xtol, ytol, 1)
