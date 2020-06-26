import re
from array import *
import numpy as np
import math
import argparse

# Initiate the parser
parser = argparse.ArgumentParser(description="Compare pod output rms and out files")

parser.add_argument("-so", "--solutionout", dest='solutionout', type=str, required=False, default='pod.out', help="solution pod.out file location name")
parser.add_argument("-sr", "--solutionrms", dest='solutionrms', type=str, required=False, default='pod.rms', help="solution pod.rms file location name")
parser.add_argument("-ro", "--runout", dest='runout', type=str, required=False, default='pod.out', help="New pod run pod.out file location name")
parser.add_argument("-rr", "--runrms", dest='runrms', type=str, required=False, default='pod.rms', help="New pod run pod.rms file location name")


def test(solutionrms, solutionout, runout, runrms):
    # Get summary stats for each satellite from pod.out and solution/pod.out and save to a list 
    
    solution_pod_out_log_location = r"solution/" + solutionout
    regex = 'RMS-XYZ ITRF CMP (G\d\d)        ([0-9.0-9?]+)        ([0-9.0-9?]+)        ([0-9.0-9?]+)'
    
    solution_pod_out = []
    with open(solution_pod_out_log_location, "r") as file:
        for line in file:
            for match in re.finditer(regex, line, re.S):
                sat_number = match.group(1)
                sat_x = float(match.group(2))
                sat_y = float(match.group(3))
                sat_z = float(match.group(4))
                solution_pod_out.append([sat_number,sat_x,sat_y,sat_z])
    
    print(solution_pod_out)

    test_pod_out_log_location = runout
    regex = 'RMS-XYZ ITRF CMP (G\d\d)        ([0-9.0-9?]+)        ([0-9.0-9?]+)        ([0-9.0-9?]+)'
    
    test_pod_out = []
    with open(test_pod_out_log_location, "r") as file:
        for line in file:
            for match in re.finditer(regex, line, re.S):
                sat_number = match.group(1)
                sat_x = float(match.group(2))
                sat_y = float(match.group(3))
                sat_z = float(match.group(4))
                test_pod_out.append([sat_number,sat_x,sat_y,sat_z])
    
    
    print(test_pod_out)

    # Run test for pod.out
    
    for test, solution in zip(test_pod_out, solution_pod_out):
        mag_diff = math.sqrt((test[1] - solution[1])**2 + (test[2] - solution[2])**2 + (test[3] - solution[3])**2)
        assert mag_diff < 0.1

        if(mag_diff > 0.1):
            print("Difference of {} found for satellite {}".format(test[0],mag_diff))
    
    # Get summart stats from pod.rms and solution/pod.rms
    
    solution_pod_rms_log_location = r"solution/" + solutionrms
    regex = 'PRN:.(.........)...........ALL:.(.[0-9.0-9?]+).(.[0-9.0-9?]+).(.[0-9.0-9?]+).(.[0-9.0-9?]+)'
    
    solution_rms_out = []
    with open(solution_pod_rms_log_location, "r") as file:
        for line in file:
            for match in re.finditer(regex, line, re.S):
                name = match.group(1).strip()
                r = float(match.group(2))
                t = float(match.group(3))
                n = float(match.group(4))
                d = float(match.group(5))
                #print("Name = {}, R = {}, T = {}, N = {}, 3D = {}".format(name,r,t,n,d))
                solution_rms_out.append([name,r,t,n,d])
    
    print(solution_rms_out)
    
    test_pod_rms_log_location = runrms
    regex = 'PRN:.(.........)...........ALL:.(.[0-9.0-9?]+).(.[0-9.0-9?]+).(.[0-9.0-9?]+).(.[0-9.0-9?]+)'
    
    test_rms_out = []
    with open(test_pod_rms_log_location, "r") as file:
        for line in file:
            for match in re.finditer(regex, line, re.S):
                name = match.group(1).strip()
                r = float(match.group(2))
                t = float(match.group(3))
                n = float(match.group(4))
                d = float(match.group(5))
                #print("Name = {}, R = {}, T = {}, N = {}, 3D = {}".format(name,r,t,n,d))
                test_rms_out.append([name,r,t,n,d])
    
    print(test_rms_out)
    
    for test_rms, solution_rms in zip(test_rms_out, solution_rms_out):
        r_diff = test_rms[1] - solution_rms[1]
        t_diff = test_rms[2] - solution_rms[2]
        n_diff = test_rms[3] - solution_rms[3]
        d_diff = test_rms[4] - solution_rms[4]
        
        assert r_diff < 0.1
        assert t_diff < 0.1
        assert n_diff < 0.1
        assert d_diff < 0.1

        if((r_diff > 0.1) or (t_diff > 0.1) or (n_diff > 0.1) or (d_diff > 0.1)):
            print("Difference of found in {} for  R = {}, T = {}, N = {}, 3D = {}".format(test_rms[0],r_diff,t_diff,n_diff,d_diff))


if __name__ == "__main__":
    
    args = parser.parse_args()
    test(args.solutionrms, args.solutionout, args.runout, args.runrms)
    print("Everything passed")
