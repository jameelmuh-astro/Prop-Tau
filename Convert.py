import sys
import os

# Function to obtain an event from the input text file
def GetEvent(file, numev):
    """
    Function to obtain an event from the input text file.

    Parameters:
        file (file object): The input text file.
        numev (int): The number of the event of interest.

    Returns:
        list: A list containing the event data.
    """
    event_data = []
    while True:
        line = file.readline()
        if not line:
            break
        evt, branch, intmult, Acurr, Zecurr, zOri, zEnd, EOri, EEnd, Dist = map(float, line.split())
        if int(evt)!= numev:
            file.seek(file.tell() - len(line))
            break
        event_data.append([zOri, zEnd, EOri, EEnd, Acurr, Zecurr, intmult])
    return event_data

if __name__ == "__main__":
    # Check if the input file name is provided as an argument
    if len(sys.argv) < 2:
        print("USAGE: python Convert.py SimPropOutput.txt")
        sys.exit(1)

    # Open the input file and read the data
    filename = sys.argv[1]
    with open(filename, 'r') as inFile:
        # Create the output file name by replacing the extension with ".out"
        outname = os.path.splitext(filename)[0] + ".out"
        with open(outname, 'w') as outFile:
            # Initialize variables
            kontev = 0
            maxev = float('inf')  # Set to infinity to process all events
            numev = 0

            # Loop through the input file to process events
            while True:
                # Get the event data
                event_data = GetEvent(inFile, numev)
                if not event_data:
                    break
                kontev += 1

                # Write the event data to the output file for events of interest
                for data in event_data:
                    z1, z2, E1, E2, A, Z, intm = data
                    if intm <= 0:
                        outFile.write(f"{z1} {E1} {A} {Z} {E2} {A} {Z}\n")

                # Check if the maximum number of events has been reached
                if kontev >= maxev:
                    break

            print(f"Events processed: {kontev}")