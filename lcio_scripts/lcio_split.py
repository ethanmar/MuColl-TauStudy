import sys
from pyLCIO import IOIMPL

input_file = "tau_gun/sim/tau_sim_5_250.slcio"
reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(input_file)

events_per_file = 2
file_index = 0
event_counter = 0
writer = None

for event in reader:
    # Open a new output file every 5 events
    if event_counter % events_per_file == 0:
        if writer:
            writer.close()
        output_file = f"split_tau_sim_5_250.{file_index:04d}.slcio"
        writer = IOIMPL.LCFactory.getInstance().createLCWriter()
        writer.open(output_file)
        file_index += 1
        print(f"Writing to {output_file}")

    writer.writeEvent(event)
    event_counter += 1

# Close last writer and reader
if writer:
    writer.close()
reader.close()

print("Done.")
