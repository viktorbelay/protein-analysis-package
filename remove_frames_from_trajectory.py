# Authors: ChatGPT, Viktor Belay, generated 052523

import mdtraj as md
import sys

def remove_frames_from_trajectory(input_file,topology, num_frames_to_remove, output_file):
    # Load the trajectory using mdtraj
    traj = md.load_dcd(input_file,top=topology)

    # Get the total number of frames in the trajectory
    total_frames = traj.n_frames

    print(traj,"with")
    print(total_frames)
    print('loaded. removing frames.....')

    # Determine the range of frames to keep
    start_frame = num_frames_to_remove
    end_frame = total_frames

    # Remove the specified number of frames from the trajectory
    traj = traj[start_frame:end_frame]

    # Save the resulting trajectory as a DCD file
    traj.save_dcd(output_file)

    print(traj,'saved!')

if __name__=='__main__':

    remove_frames_from_trajectory(sys.argv[1],sys.argv[2],int(sys.argv[3]),sys.argv[4])