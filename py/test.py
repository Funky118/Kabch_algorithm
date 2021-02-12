# A script creating test batches
#
# Open in IDLE and type run_test(x), for x number of randomly generated test batches
#
# Create a random number of points in each point cloud and save it
# to a data_XX.txt and data_trans_XX.txt files indexed from 0 to number_of_data.
#
# Then send each file as an argument into MIPO. When MIPO is finished it will output
# data_XX_out.txt and DUMP_data_XX.txt.
# 
# In this script, data_XX_out.txt will be compared to data_XX.txt and a diff_XX.txt will be created.
# Finally, evaluate diff_XX.txt and create FAILED_data_XX.txt for any test that has failed.
#
# Notes:    Does not create translation (shouldn't be a problem as matrix T is only a bonus product of the Kabsch algorithm)
#           

import os
import subprocess as sp
import numpy as np
from scipy.spatial.transform import Rotation as R

def create_batch(dir, number_of_data):

    new_dir = dir

    for i in range(number_of_data):
        # Randomize the number of points in the point cloud
        num_of_points = np.random.randint(4, 20)
        data = np.random.randint(low = 0, high = 10, size = (num_of_points, 3))

        # Randomize rotation
        rot_vec = R.random()
        trans = rot_vec.apply(data)

        # Save the data into strings
        data_string = str(num_of_points) + ';' + '\n'.join(','.join('%.5f' %x for x in y) for y in data) + ';'
        trans_string = str(num_of_points) + ';' + '\n'.join(','.join('%.5f' %x for x in y) for y in trans) + ';'
        
        # Save the rotation vector too
        trans_string += '\n Rot vector:\n' + '\n'.join('%.5f' %x for x in rot_vec.as_rotvec())

        # Save the strings into .txt files
        dfile = open(new_dir + '\\data_' + str(i) + '.txt', "w")
        tfile = open(new_dir + '\\data_'+ str(i) + '_trans.txt', "w")
        dfile.write(data_string)
        tfile.write(trans_string)

        dfile.close()
        tfile.close()

def run_MIPO(path, data, data_trans, MIPO = "main"):
    # Runs the main executable with the two pointclouds as its arguments
    sp.run(["..\\bin\\" + MIPO + ".exe", path + data + ".txt", path + data_trans + ".txt"])

def eval(path, number_of_data, tol):
    # Parses the data_XX.txt and data_XX_out.txt files into two lists and compares them
    for num in range(number_of_data):
        d_in = open(path +'\\data_' + str(num) + '.txt', "r")
        d_out = open(path +'\\data_' + str(num) + '_out.txt', "r")
        diff = open(path +'\\data_' + str(num) + '_diff.txt', "w")
        d_dump = open(path +'\\DUMP_data_' + str(num) + '.txt', "r")
        
        # Ignores everything up to the first ';' which delimits the important data
        while d_in.read(1) != ';':
            pass
        while d_out.read(1) != ';':
            pass

        # Copy the whole file into a string
        d_in_string = d_in.readlines()
        d_out_string = d_out.readlines()

        # Generate preprocessed lists
        d_in_string = [x.replace(';', '') for x in d_in_string]
        d_out_string = [x.replace(';', '') for x in d_out_string]

        # Process the list elements into a form that can be converted to floats
        d_in_float = []
        d_out_float = []
        for i in range(len(d_in_string)):
            d_in_float.append(d_in_string[i].split(','))
            d_out_float.append(d_out_string[i].split(','))

        # Change both lists to float for difference comparison
        d_in_float = [[float(i) for i in x] for x in d_in_float]
        d_out_float = [[float(i) for i in x] for x in d_out_float]

        # Create a new file data_XX_diff.txt
        diff_string = ''
        for k, j in zip(d_in_float, d_out_float):
            for input_data, output_data in zip(k, j):
                if not os.path.exists(path + '\\FAILED_data_' + str(num) + '.txt') and ((input_data - output_data) < -tol or (input_data - output_data) > tol):
                    open(path + '\\FAILED_data_' + str(num) + '.txt', "w").writelines(d_dump.readlines())
                # Print the difference up to 12 decimal points
                diff_string += str(round((input_data - output_data), 12)) + ' '
            diff_string += '\n'
        diff.write(diff_string)
        
def run_test(number_of_point_clouds, tol = 0.1):

    # Check if other tests have been run and create a new dir for the new test
    error_max_dirs = True
    for i in range(100):
        if not (os.path.exists('..\\test'+str(i))):
            new_dir = 'test'+str(i)
            os.mkdir('..\\' + new_dir)
            error_max_dirs = False
            break

    if error_max_dirs:
            raise ValueError('Maximum number of test directories has been reached.')
            
    create_batch('..\\' + new_dir, number_of_point_clouds)
    
    for i in range(number_of_point_clouds):
        run_MIPO(new_dir + '/', 'data_'+str(i), 'data_'+str(i)+'_trans')
        
    eval('..\\' + new_dir, number_of_point_clouds, tol)
