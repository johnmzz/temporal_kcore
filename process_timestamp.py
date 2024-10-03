# converts original timestamps to a sequence of continuous integers

def process_temporal_graph(input_file, output_file):
    # Step 1: Read the file and collect all timestamps
    edges = []
    timestamps = set()

    with open(input_file, 'r') as f:
        for line in f:
            start, end, timestamp = line.split()
            edges.append((int(start), int(end), int(timestamp)))
            timestamps.add(int(timestamp))

    # Step 2: Create a mapping from timestamps to continuous integers
    sorted_timestamps = sorted(timestamps)
    timestamp_mapping = {ts: idx + 1 for idx, ts in enumerate(sorted_timestamps)}

    # Step 3: Transform the edges with new timestamp values
    with open(output_file, 'w') as f:
        for start, end, timestamp in edges:
            new_timestamp = timestamp_mapping[timestamp]
            f.write(f"{start} {end} {new_timestamp}\n")

# Example usage
input_file = 'data/youtube.txt'  # Replace with the path to your input file
output_file = 'data_sequenced/youtube.txt'  # Output file path
process_temporal_graph(input_file, output_file)
