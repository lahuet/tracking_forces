This code calculates dynamic reaction forces from tracked whisker data.


    -- The process of calculating reaction forces from image data can be 
       done one step at a time or all at once. To run the entire process,
       place the .mat file of raw data in the raw_data folder and run

       $ python tracking_forces.py data_file.mat

       from the main directory. The results will show up in the output_data
       folder, including .p and .mat files and plots. This will use all of 
       the default settings. All of the subprocesses will print their
       progress.

       Settings and options are stored in the configuration 
       file tracking_forces.cfg and should be set before any calculations.


   -- The individual steps of the process are:

        - To convert the tracked point data to sampled points and a 
          configuration in terms of joint angles, run

          $ python convert.py data_file

          the results will be saved in the converted_data folder. This process
          can take some time. See file for optional arguments.


        - To filter the angle trajectories and compute the velocities and
          accelerations needed for the constraint force calculations, run

          $ python qfilter.py data_file

          the results will also show up in the converted data folder.


        - To compute the constraint forces for all time steps, run
          
          $ python forces.py data_file


        - To plot the results run

          $ python plot.py data_file

          the results will be saved in the output_data folder as .png images.  

	- To animate the whisker's motion run
	 
	  $ python animate.py data_file

	  if specified, an .mp4 movie will be saved in the output_data folder.
